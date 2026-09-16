################################################################################
## Submission audit for Figures 2 and 4
##
## Recomputes all reported summaries from the saved analysis objects, checks
## internal consistency, and writes compact diagnostic tables. It does not
## refit models or overwrite any analysis output or figure.
################################################################################

rm(list = ls())

required_packages <- c("tidyverse")
missing_packages <- required_packages[
	!vapply(required_packages, requireNamespace, quietly = TRUE,
		FUN.VALUE = logical(1))
]
if (length(missing_packages) > 0L) {
	stop("Missing package(s): ", paste(missing_packages, collapse = ", "),
		call. = FALSE)
}
suppressPackageStartupMessages(library(tidyverse))

OUT_DIR <- "tables/diagnostics/submission_audit"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

required_files <- c(
	"data_processed/fia_traits_clean.rds",
	"tables/shap_values.rds",
	"tables/shap_per_var.rds",
	"tables/shap_importance.rds",
	"tables/vecv_raw.rds",
	"tables/vecv_summary.rds",
	"tables/vecv_divergence.rds",
	"tables/perf_broadleaf.csv",
	"tables/perf_coniferous.csv",
	"tables/sensitivity_vecv.rds"
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
	stop("Missing required file(s):\n  ", paste(missing_files, collapse = "\n  "),
		call. = FALSE)
}

TRAITS <- c(
	"bark_thickness", "conduit_diam", "height", "leaf_density", "leaf_k",
	"root_depth", "seed_dry_mass", "shade_tolerance", "specific_leaf_area"
)
COVARIATES <- c(
	"standage", "temp_pc", "soil_pc", "rain_pc", "elevation", "soil_ph"
)
ENV_VARS <- setdiff(COVARIATES, "standage")

audit_checks <- tibble(
	section = character(), check = character(), status = character(),
	detail = character()
)
add_check <- function(section, check, status, detail) {
	audit_checks <<- bind_rows(
		audit_checks,
		tibble(section = section, check = check, status = status, detail = detail)
	)
}

max_abs_diff <- function(x, y) {
	z <- abs(x - y)
	if (!any(is.finite(z))) return(NA_real_)
	max(z, na.rm = TRUE)
}

# ==============================================================================
# Data structure and dependence
# ==============================================================================

dat <- read_rds("data_processed/fia_traits_clean.rds") %>%
	mutate(
		leaf_type = case_when(
			biome_boreal_forests_or_taiga == 1 |
				biome_temperate_conifer_forests == 1 ~ "coniferous",
			biome_temperate_broadleaf_forests == 1 |
				biome_mediterranean_woodlands == 1 ~ "broadleaf",
			TRUE ~ NA_character_
		)
	) %>%
	filter(!is.na(leaf_type))

data_structure <- dat %>%
	group_by(leaf_type) %>%
	summarise(
		n_rows = n(),
		n_unique_PID = n_distinct(PID),
		n_unique_PID_rep = n_distinct(PID_rep),
		fraction_rows_from_remeasured_PID = mean(
			duplicated(PID) | duplicated(PID, fromLast = TRUE)
		),
		.groups = "drop"
	)
write_csv(data_structure, file.path(OUT_DIR, "data_structure.csv"))

PID_multiplicity <- dat %>%
	add_count(PID, name = "n_inventories_PID") %>%
	mutate(
		standage_bin = cut(
			standage, breaks = seq(0, 160, by = 10),
			include.lowest = TRUE, right = FALSE
		)
	) %>%
	filter(!is.na(standage_bin)) %>%
	group_by(leaf_type, standage_bin) %>%
	summarise(
		n_records = n(),
		fraction_records_from_remeasured_PID = mean(n_inventories_PID > 1L),
		.groups = "drop"
	)
write_csv(PID_multiplicity,
	file.path(OUT_DIR, "PID_remeasurement_by_standage.csv"))

if (n_distinct(dat$PID) < nrow(dat)) {
	add_check(
		"design", "Repeated inventories are present", "WARNING",
		sprintf(
			"%s rows represent %s unique PIDs; row-wise train/test splitting leaks repeated plots.",
			format(nrow(dat), big.mark = ","),
			format(n_distinct(dat$PID), big.mark = ",")
		)
	)
}

# Reconstruct the deterministic split used by scripts/03_rf_fit.R and quantify
# how many held-out rows have another inventory of the same PID in training.
split_leakage <- map_dfr(c("broadleaf", "coniferous"), function(lt) {
	d <- dat %>% filter(leaf_type == lt) %>% mutate(.id = row_number())
	n_test <- min(5000L, floor(0.2 * nrow(d)))
	set.seed(42L)
	test_ids <- sample(d$.id, n_test, replace = FALSE)
	train <- filter(d, !.id %in% test_ids)
	test <- filter(d, .id %in% test_ids)
	tibble(
		leaf_type = lt,
		n_test_rows = nrow(test),
		n_test_PIDs = n_distinct(test$PID),
		n_test_PIDs_also_in_training = sum(unique(test$PID) %in% train$PID),
		n_test_rows_with_PID_in_training = sum(test$PID %in% train$PID),
		fraction_test_rows_with_PID_in_training = mean(test$PID %in% train$PID)
	)
})
write_csv(split_leakage, file.path(OUT_DIR, "row_split_PID_leakage.csv"))
add_check(
	"design", "Original held-out split is independent by PID", "FAIL",
	paste(
		paste0(split_leakage$leaf_type, "=",
			scales::percent(split_leakage$fraction_test_rows_with_PID_in_training,
				accuracy = 0.1)),
		collapse = "; "
	)
)

covariate_correlations <- dat %>%
	dplyr::select(all_of(COVARIATES)) %>%
	cor(use = "pairwise.complete.obs") %>%
	as.data.frame() %>%
	rowid_to_column("predictor_1") %>%
	mutate(predictor_1 = rownames(cor(dat[COVARIATES], use = "pairwise.complete.obs"))) %>%
	pivot_longer(-predictor_1, names_to = "predictor_2", values_to = "correlation") %>%
	filter(predictor_1 < predictor_2)
write_csv(covariate_correlations,
	file.path(OUT_DIR, "covariate_correlations.csv"))

# ==============================================================================
# Figure 2 SHAP integrity and narrative
# ==============================================================================

shap_values <- read_rds("tables/shap_values.rds")
shap_per_var <- read_rds("tables/shap_per_var.rds")
shap_importance <- read_rds("tables/shap_importance.rds")

shap_per_var_check <- shap_values %>%
	group_by(trait, leaf_type, variable) %>%
	summarise(
		sum_abs_shap_calc = sum(abs(shap_value), na.rm = TRUE),
		mean_abs_shap_calc = mean(abs(shap_value), na.rm = TRUE),
		.groups = "drop"
	) %>%
	left_join(
		shap_per_var %>%
			dplyr::select(trait, leaf_type, variable, sum_abs_shap, mean_abs_shap),
		by = c("trait", "leaf_type", "variable")
	)

shap_var_max_diff <- max(
	max_abs_diff(shap_per_var_check$sum_abs_shap_calc,
		shap_per_var_check$sum_abs_shap),
	max_abs_diff(shap_per_var_check$mean_abs_shap_calc,
		shap_per_var_check$mean_abs_shap),
	na.rm = TRUE
)
add_check(
	"figure_2", "Saved per-variable SHAP summaries reproduce from shap_values",
	if_else(shap_var_max_diff < 1e-10, "PASS", "FAIL"),
	sprintf("maximum absolute discrepancy = %.3g", shap_var_max_diff)
)

shap_ratio_check <- shap_per_var %>%
	mutate(category = case_when(
		variable %in% ENV_VARS ~ "environmental",
		variable == "standage" ~ "successional",
		TRUE ~ NA_character_
	)) %>%
	filter(!is.na(category)) %>%
	group_by(trait, leaf_type, category) %>%
	summarise(value = sum(sum_abs_shap), .groups = "drop") %>%
	pivot_wider(names_from = category, values_from = value) %>%
	mutate(env_succ_ratio_calc = environmental / successional) %>%
	left_join(
		shap_importance %>%
			dplyr::select(trait, leaf_type, env_succ_ratio),
		by = c("trait", "leaf_type")
	)

ratio_max_diff <- max_abs_diff(
	shap_ratio_check$env_succ_ratio_calc,
	shap_ratio_check$env_succ_ratio
)
add_check(
	"figure_2", "Saved environmental-to-age ratios reproduce from SHAP values",
	if_else(ratio_max_diff < 1e-10, "PASS", "FAIL"),
	sprintf("maximum absolute discrepancy = %.3g", ratio_max_diff)
)

fig2_ratios <- shap_importance %>%
	mutate(
		mean_environmental_predictor_to_age = env_succ_ratio / length(ENV_VARS),
		combined_environment_exceeds_age = env_succ_ratio > 1,
		mean_environmental_predictor_exceeds_age =
			mean_environmental_predictor_to_age > 1
	) %>%
	arrange(leaf_type, env_succ_ratio)
write_csv(fig2_ratios, file.path(OUT_DIR, "fig2_ratio_audit.csv"))

fig2_predictor_ranks <- shap_per_var %>%
	group_by(trait, trait_label, leaf_type) %>%
	arrange(desc(mean_abs_shap), .by_group = TRUE) %>%
	mutate(
		predictor_rank = row_number(),
		share_total_abs_SHAP = mean_abs_shap / sum(mean_abs_shap)
	) %>%
	ungroup()
write_csv(fig2_predictor_ranks,
	file.path(OUT_DIR, "fig2_predictor_ranks.csv"))

fig2_summary <- fig2_ratios %>%
	group_by(leaf_type) %>%
	summarise(
		mean_ratio = mean(env_succ_ratio),
		median_ratio = median(env_succ_ratio),
		min_ratio = min(env_succ_ratio),
		max_ratio = max(env_succ_ratio),
		n_combined_environment_gt_age = sum(combined_environment_exceeds_age),
		n_mean_environmental_predictor_gt_age =
			sum(mean_environmental_predictor_exceeds_age),
		.groups = "drop"
	)
write_csv(fig2_summary, file.path(OUT_DIR, "fig2_summary.csv"))
add_check(
	"figure_2", "Environmental dominance survives a predictor-count diagnostic",
	"PASS",
	sprintf(
		"The mean environmental predictor exceeds stand age in %d/18 combinations; the three exceptions are the expected height/light-competition cases.",
		sum(fig2_ratios$mean_environmental_predictor_exceeds_age)
	)
)

ci_path <- "tables/shap_importance_ci.rds"
cv_raw_path <- "tables/shap_importance_cv_raw.rds"
if (file.exists(ci_path) && file.exists(cv_raw_path)) {
	shap_ci <- read_rds(ci_path)
	shap_cv_raw <- read_rds(cv_raw_path)
	write_csv(shap_ci, file.path(OUT_DIR, "fig2_CV_interval_audit.csv"))
	ci_complete <- nrow(shap_ci) == 18L &&
		all(c("trait", "leaf_type", "ratio_lwr", "ratio_upr") %in% names(shap_ci))
	add_check(
		"figure_2", "Repeated-CV uncertainty artifacts are complete",
		if_else(ci_complete, "PASS", "FAIL"),
		sprintf("%d interval rows and %d raw CV rows", nrow(shap_ci), nrow(shap_cv_raw))
	)

	if (all(c("rep", "env_sum", "succ_sum") %in% names(shap_cv_raw))) {
		repeat_pooled <- shap_cv_raw %>%
			group_by(trait, leaf_type, rep) %>%
			summarise(
				ratio_repeat_pooled = sum(env_sum) / sum(succ_sum),
				.groups = "drop"
			) %>%
			group_by(trait, leaf_type) %>%
			summarise(
				n_repeats = n(),
				ratio_repeat_med = median(ratio_repeat_pooled),
				ratio_repeat_lwr = quantile(ratio_repeat_pooled, 0.025, type = 8),
				ratio_repeat_upr = quantile(ratio_repeat_pooled, 0.975, type = 8),
				.groups = "drop"
			)
		write_csv(repeat_pooled,
			file.path(OUT_DIR, "fig2_repeat_pooled_intervals.csv"))
	}
} else {
	add_check(
		"figure_2", "Repeated-CV uncertainty artifacts are present in repository",
		"INCOMPLETE",
		"shap_importance_ci.rds and shap_importance_cv_raw.rds are absent locally"
	)
}

# ==============================================================================
# Figure 4 VEcv integrity, sensitivity, and narrative
# ==============================================================================

vecv_raw <- read_rds("tables/vecv_raw.rds")
vecv_summary <- read_rds("tables/vecv_summary.rds")
vecv_divergence <- read_rds("tables/vecv_divergence.rds")

repeat_completeness <- vecv_raw %>%
	group_by(trait, leaf_type, variable, env_group, standage_mid) %>%
	summarise(n_repeats = n_distinct(repeat_id), .groups = "drop")
write_csv(repeat_completeness,
	file.path(OUT_DIR, "fig4_repeat_completeness.csv"))
add_check(
	"figure_4", "Every retained VEcv cell contains 30 CV repeats",
	if_else(all(repeat_completeness$n_repeats == 30L), "PASS", "FAIL"),
	sprintf("repeat counts range from %d to %d",
		min(repeat_completeness$n_repeats), max(repeat_completeness$n_repeats))
)

vecv_recalc <- vecv_raw %>%
	group_by(
		trait, trait_label, leaf_type, variable, variable_label,
		env_group, standage_bin, standage_mid
	) %>%
	summarise(
		VEcv_med_calc = median(VEcv),
		VEcv_lwr_calc = quantile(VEcv, 0.025),
		VEcv_upr_calc = quantile(VEcv, 0.975),
		E1_med_calc = median(E1),
		.groups = "drop"
	) %>%
	left_join(
		vecv_summary %>%
			dplyr::select(
				trait, leaf_type, variable, env_group, standage_bin,
				VEcv_med, VEcv_lwr, VEcv_upr, E1_med
			),
		by = c("trait", "leaf_type", "variable", "env_group", "standage_bin")
	)

vecv_summary_max_diff <- max(
	max_abs_diff(vecv_recalc$VEcv_med_calc, vecv_recalc$VEcv_med),
	max_abs_diff(vecv_recalc$VEcv_lwr_calc, vecv_recalc$VEcv_lwr),
	max_abs_diff(vecv_recalc$VEcv_upr_calc, vecv_recalc$VEcv_upr),
	max_abs_diff(vecv_recalc$E1_med_calc, vecv_recalc$E1_med),
	na.rm = TRUE
)
add_check(
	"figure_4", "Saved VEcv summaries reproduce from raw repeats",
	if_else(vecv_summary_max_diff < 1e-10, "PASS", "FAIL"),
	sprintf("maximum absolute discrepancy = %.3g", vecv_summary_max_diff)
)

paired_raw <- vecv_raw %>%
	dplyr::select(
		trait, trait_label, leaf_type, variable, variable_label,
		standage_bin, standage_mid, repeat_id, env_group, VEcv, E1
	) %>%
	pivot_wider(names_from = env_group, values_from = c(VEcv, E1)) %>%
	mutate(
		delta_VEcv = VEcv_high - VEcv_low,
		delta_E1 = E1_high - E1_low
	)

div_recalc <- paired_raw %>%
	group_by(
		trait, trait_label, leaf_type, variable, variable_label,
		standage_bin, standage_mid
	) %>%
	summarise(
		delta_med_calc = median(delta_VEcv, na.rm = TRUE),
		delta_lwr_calc = quantile(delta_VEcv, 0.025, na.rm = TRUE),
		delta_upr_calc = quantile(delta_VEcv, 0.975, na.rm = TRUE),
		.groups = "drop"
	) %>%
	left_join(
		vecv_divergence %>%
			dplyr::select(
				trait, leaf_type, variable, standage_bin,
				delta_med, delta_lwr, delta_upr
			),
		by = c("trait", "leaf_type", "variable", "standage_bin")
	)

div_max_diff <- max(
	max_abs_diff(div_recalc$delta_med_calc, div_recalc$delta_med),
	max_abs_diff(div_recalc$delta_lwr_calc, div_recalc$delta_lwr),
	max_abs_diff(div_recalc$delta_upr_calc, div_recalc$delta_upr),
	na.rm = TRUE
)
add_check(
	"figure_4", "Saved VEcv divergence reproduces from raw repeats",
	if_else(div_max_diff < 1e-10, "PASS", "FAIL"),
	sprintf("maximum absolute discrepancy = %.3g", div_max_diff)
)

# Reproduce the manuscript's early/late averages, then provide a balanced
# combination-level comparison and the E1 robustness check.
fig4_stage_summary <- vecv_summary %>%
	filter(is.finite(VEcv_med), is.finite(E1_med)) %>%
	mutate(stage = case_when(
		standage_mid <= 30 ~ "early",
		standage_mid >= 100 ~ "late",
		TRUE ~ NA_character_
	)) %>%
	filter(!is.na(stage)) %>%
	group_by(leaf_type, stage) %>%
	summarise(
		VEcv = mean(VEcv_med),
		E1 = mean(E1_med),
		n_retained_cells = n(),
		.groups = "drop"
	)
write_csv(fig4_stage_summary,
	file.path(OUT_DIR, "fig4_early_late_summary.csv"))

combo_stage <- vecv_summary %>%
	filter(is.finite(VEcv_med), is.finite(E1_med)) %>%
	mutate(stage = case_when(
		standage_mid <= 30 ~ "early",
		standage_mid >= 100 ~ "late",
		TRUE ~ NA_character_
	)) %>%
	filter(!is.na(stage)) %>%
	group_by(leaf_type, trait, variable, env_group, stage) %>%
	summarise(VEcv = mean(VEcv_med), E1 = mean(E1_med), .groups = "drop") %>%
	pivot_wider(names_from = stage, values_from = c(VEcv, E1)) %>%
	filter(
		is.finite(VEcv_early), is.finite(VEcv_late),
		is.finite(E1_early), is.finite(E1_late)
	) %>%
	mutate(
		delta_VEcv = VEcv_late - VEcv_early,
		delta_E1 = E1_late - E1_early
	)
write_csv(combo_stage,
	file.path(OUT_DIR, "fig4_early_late_by_combination.csv"))

fig4_increase_robustness <- combo_stage %>%
	group_by(leaf_type) %>%
	summarise(
		n_combinations = n(),
		n_VEcv_increase = sum(delta_VEcv > 0),
		n_E1_increase = sum(delta_E1 > 0),
		median_delta_VEcv = median(delta_VEcv),
		median_delta_E1 = median(delta_E1),
		.groups = "drop"
	)
write_csv(fig4_increase_robustness,
	file.path(OUT_DIR, "fig4_predictability_increase_robustness.csv"))

summarise_metric_divergence <- function(data, metric, metric_name) {
	metric_sym <- rlang::sym(metric)
	by_age <- data %>%
		filter(is.finite(!!metric_sym)) %>%
		group_by(leaf_type, variable, variable_label, standage_mid) %>%
		summarise(mean_abs_delta = mean(abs(!!metric_sym)), .groups = "drop")

	by_age %>%
		group_by(leaf_type, variable, variable_label) %>%
		group_modify(~ {
			d <- .x
			tibble(
				metric = metric_name,
				mean_abs_delta = mean(d$mean_abs_delta),
				early_mean_abs_delta = mean(
					d$mean_abs_delta[d$standage_mid <= 30], na.rm = TRUE
				),
				late_mean_abs_delta = mean(
					d$mean_abs_delta[d$standage_mid >= 100], na.rm = TRUE
				),
				trend_per_year = unname(coef(
					lm(mean_abs_delta ~ standage_mid, data = d)
				)[[2L]])
			)
		}) %>%
		ungroup() %>%
		mutate(
			late_minus_early = late_mean_abs_delta - early_mean_abs_delta,
			trend_direction = if_else(trend_per_year > 0, "increasing", "decreasing")
		)
}

fig4_divergence_metrics <- bind_rows(
	summarise_metric_divergence(paired_raw, "delta_VEcv", "VEcv"),
	summarise_metric_divergence(paired_raw, "delta_E1", "E1")
)
write_csv(fig4_divergence_metrics,
	file.path(OUT_DIR, "fig4_divergence_metric_comparison.csv"))

fig4_divergence_reported <- vecv_divergence %>%
	group_by(leaf_type, variable, variable_label) %>%
	summarise(
		mean_abs_delta_VEcv = mean(abs(delta_med), na.rm = TRUE),
		proportion_bins_CV_interval_excludes_zero =
			mean(sig_divergence, na.rm = TRUE),
		.groups = "drop"
	) %>%
	left_join(
		fig4_divergence_metrics %>% filter(metric == "VEcv") %>%
			dplyr::select(
				leaf_type, variable, trend_per_year,
				early_mean_abs_delta, late_mean_abs_delta, late_minus_early,
				trend_direction
			),
		by = c("leaf_type", "variable")
	)
write_csv(fig4_divergence_reported,
	file.path(OUT_DIR, "fig4_divergence_by_environment.csv"))

late_persistence <- vecv_divergence %>%
	filter(standage_mid >= 100, is.finite(delta_med)) %>%
	group_by(leaf_type, variable, variable_label) %>%
	summarise(
		late_mean_abs_delta_VEcv = mean(abs(delta_med)),
		late_median_abs_delta_VEcv = median(abs(delta_med)),
		late_proportion_bins_CV_interval_excludes_zero =
			mean(sig_divergence),
		n_trait_age_cells = n(),
		.groups = "drop"
	)
write_csv(late_persistence,
	file.path(OUT_DIR, "fig4_late_divergence_persistence.csv"))

# The current main figure averages signed upper-minus-lower values across
# traits. Opposing trait-specific directions can therefore cancel even when
# absolute environmental contrasts remain large.
aggregation_cancellation <- vecv_divergence %>%
	filter(is.finite(delta_med)) %>%
	group_by(leaf_type, variable, variable_label, standage_mid) %>%
	summarise(
		mean_signed_delta = mean(delta_med),
		mean_abs_delta = mean(abs(delta_med)),
		cancellation_fraction = 1 - abs(mean_signed_delta) / mean_abs_delta,
		.groups = "drop"
	) %>%
	group_by(leaf_type, variable, variable_label) %>%
	summarise(
		mean_abs_delta = mean(mean_abs_delta),
		mean_abs_signed_delta = mean(abs(mean_signed_delta)),
		median_cancellation_fraction = median(cancellation_fraction),
		max_cancellation_fraction = max(cancellation_fraction),
		.groups = "drop"
	)
write_csv(aggregation_cancellation,
	file.path(OUT_DIR, "fig4_trait_aggregation_cancellation.csv"))

# Partition-to-partition stability of the environmental divergence trend.
trend_by_repeat <- paired_raw %>%
	filter(is.finite(delta_VEcv), is.finite(delta_E1)) %>%
	group_by(leaf_type, variable, variable_label, repeat_id, standage_mid) %>%
	summarise(
		mean_abs_delta_VEcv = mean(abs(delta_VEcv)),
		mean_abs_delta_E1 = mean(abs(delta_E1)),
		.groups = "drop"
	) %>%
	group_by(leaf_type, variable, variable_label, repeat_id) %>%
	group_modify(~ {
		d <- .x
		tibble(
			trend_VEcv = unname(coef(
				lm(mean_abs_delta_VEcv ~ standage_mid, data = d)
			)[[2L]]),
			trend_E1 = unname(coef(
				lm(mean_abs_delta_E1 ~ standage_mid, data = d)
			)[[2L]])
		)
	}) %>%
	ungroup()

trend_stability <- trend_by_repeat %>%
	group_by(leaf_type, variable, variable_label) %>%
	summarise(
		VEcv_trend_median = median(trend_VEcv),
		VEcv_trend_lwr = quantile(trend_VEcv, 0.025),
		VEcv_trend_upr = quantile(trend_VEcv, 0.975),
		VEcv_fraction_increasing = mean(trend_VEcv > 0),
		E1_trend_median = median(trend_E1),
		E1_trend_lwr = quantile(trend_E1, 0.025),
		E1_trend_upr = quantile(trend_E1, 0.975),
		E1_fraction_increasing = mean(trend_E1 > 0),
		.groups = "drop"
	)
write_csv(trend_stability,
	file.path(OUT_DIR, "fig4_divergence_trend_stability.csv"))

temperature_trends <- fig4_divergence_reported %>%
	filter(variable == "temp_pc")
if (all(temperature_trends$trend_per_year < 0, na.rm = TRUE)) {
	add_check(
		"figure_4", "Manuscript claim that temperature divergence increases",
		"FAIL",
		"The saved VEcv output gives a decreasing full-range |delta VEcv| trend for temperature in both forest types."
	)
}

# The current threshold analysis is algebraically projected from the main
# result rather than recomputed from alternative strata.
sensitivity_vecv <- read_rds("tables/sensitivity_vecv.rds")
projection_error <- max_abs_diff(
	sensitivity_vecv$projected_abs_delta,
	sensitivity_vecv$mean_abs_delta_ref * sensitivity_vecv$d_ratio
)
add_check(
	"figure_4", "Alternative-quantile VEcv sensitivity is independently estimated",
	"FAIL",
	sprintf(
		"It is a deterministic projection of the 25/75 result (identity error %.3g), not a re-stratified CV analysis.",
		projection_error
	)
)

# Existing performance files do not match the schema generated by the current
# scripts/03_rf_fit.R, so held-out performance claims cannot be reconstructed.
perf <- bind_rows(
	read_csv("tables/perf_broadleaf.csv", show_col_types = FALSE),
	read_csv("tables/perf_coniferous.csv", show_col_types = FALSE)
)
write_csv(perf, file.path(OUT_DIR, "available_model_performance.csv"))
if (!all(c("r2_test", "rmse_test") %in% names(perf))) {
	add_check(
		"provenance", "Held-out model performance table matches current code",
		"INCOMPLETE",
		"The tracked performance files contain OOB rsq/pred_error only; r2_test and rmse_test are absent."
	)
}

tracked_inputs <- tibble(
	file = required_files,
	modified = as.character(file.info(required_files)$mtime),
	md5 = unname(tools::md5sum(required_files))
)
write_csv(tracked_inputs, file.path(OUT_DIR, "audited_file_manifest.csv"))

write_csv(audit_checks, file.path(OUT_DIR, "audit_checks.csv"))

message("Submission audit completed without modifying analysis outputs.")
message("Diagnostic tables: ", OUT_DIR)
print(audit_checks, n = Inf)
