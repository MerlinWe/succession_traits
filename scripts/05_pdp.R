################################################################################
## succession_traits: 05 — Partial dependence analysis
## Tests for spatio-temporal interactions between environmental filtering and
## successional filtering (RQ2). For each environmental variable, data are
## stratified into high/low quantile groups and new RF models are fitted to
## each stratum. Partial dependence curves for stand age are then computed
## and compared between strata.
##
## Rationale for refitting models on strata rather than using global models:
## Random forests do not explicitly model interaction terms. Stratified models
## restrict predictions to observed feature combinations (avoiding extrapolation),
## capture heterogeneous effects where environmental predictors have opposing
## effects at different gradient extremes, and allow direct assessment of model
## sensitivity across environmental conditions.
##
## Bootstrap uncertainty (n = N_BOOT):
## The bootstrap resamples which observations enter each stratum, giving a
## distribution of slope and intercept differences across resamples. Inference
## is based on whether the 95% bootstrap CI of the signed difference excludes
## zero — a difference is considered robust if the CI does not include zero
## across resamples. This avoids the problems of formal significance testing
## on bootstrapped quantities (inflated power, p-values that reflect resampling
## variance rather than ecological effect size).
##
## Input:  data_processed/fia_traits_clean.rds
##         tables/perf_broadleaf.csv / perf_coniferous.csv  (hyperparameters)
##         tables/shap_per_var.rds                          (for SHAP weighting)
##
## Output: tables/pdp_raw.rds          (all bootstrap PDP curves)
##         tables/pdp_stats.rds        (slope/intercept diffs per iteration)
##         tables/pdp_summary.rds      (medians + CIs + robust flags)
##         plots/pdp/fig2b_ellipse.png (RQ2 main figure candidate)
##
## Author: M. Weiss @ Maynard Lab UCL / ETH Zürich
################################################################################

rm(list = ls())
set.seed(42)

# ── Libraries ─────────────────────────────────────────────────────────────────
library(ranger)
library(doParallel)
library(foreach)
library(ggh4x)
library(tidyverse)

# ── Configuration ─────────────────────────────────────────────────────────────
N_BOOT      <- 100L
PROBS_MAIN  <- c(0.25, 0.75)
PROBS_SENS  <- c(0.10, 0.90)
RUN_SENS    <- FALSE
PARALLEL    <- TRUE
N_CORES     <- 32L

PATH_DATA   <- "data_processed/fia_traits_clean.rds"
PATH_TABLES <- "tables"
PATH_PLOTS  <- "figures/supplementary/pdp"

source("scripts/functions.R")
source("scripts/plot_theme.R")

# ── Vocabulary ────────────────────────────────────────────────────────────────

TRAITS <- c("bark_thickness", "conduit_diam", "height", "leaf_density",
						 "leaf_k", "root_depth", "seed_dry_mass", "shade_tolerance",
						 "specific_leaf_area")
COVARIATES <- c("standage", "temp_pc", "soil_pc", "rain_pc", "elevation", "soil_ph")
LEAF_TYPES <- c("broadleaf", "coniferous")
ENV_VARS   <- c("temp_pc", "soil_pc", "rain_pc", "elevation", "soil_ph")

# ══════════════════════════════════════════════════════════════════════════════
# 1. Load data and hyperparameters
# ══════════════════════════════════════════════════════════════════════════════

data <- read_rds(PATH_DATA)

data <- data %>%
	mutate(
		leaf_type = case_when(
			biome_boreal_forests_or_taiga    == 1 |
				biome_temperate_conifer_forests  == 1 ~ "coniferous",
			biome_temperate_broadleaf_forests == 1 |
				biome_mediterranean_woodlands    == 1 ~ "broadleaf",
			TRUE ~ NA_character_
		)
	) %>%
	filter(!is.na(leaf_type))

message(sprintf("Data: %d plots (%d broadleaf / %d coniferous)",
								nrow(data),
								sum(data$leaf_type == "broadleaf"),
								sum(data$leaf_type == "coniferous")))

hyper_params <- map(LEAF_TYPES, function(lt) {
	read_csv(file.path(PATH_TABLES, sprintf("perf_%s.csv", lt)),
					 show_col_types = FALSE) %>%
		dplyr::select(trait, num_trees, mtry, min_node_size)
}) %>% set_names(LEAF_TYPES)

shap_per_var <- read_rds(file.path(PATH_TABLES, "shap_per_var.rds"))


# ══════════════════════════════════════════════════════════════════════════════
# 2. Core PDP bootstrap function
# ══════════════════════════════════════════════════════════════════════════════

run_pdp_once <- function(boot_data, covariates, traits, hyper_grid, probs) {
	
	qtxt <- if (identical(probs, c(0.10, 0.90))) "10%" else "25%"
	
	q_labels <- list(
		temp_pc   = c(paste0("Cold (lower ",       qtxt, ")"),
									paste0("Warm (upper ",        qtxt, ")")),
		soil_pc   = c(paste0("Sandy (lower ",       qtxt, ")"),
									paste0("Water-retentive (upper ", qtxt, ")")),
		rain_pc   = c(paste0("Low precip (lower ",  qtxt, ")"),
									paste0("High precip (upper ", qtxt, ")")),
		elevation = c(paste0("Low elevation (lower ",  qtxt, ")"),
									paste0("High elevation (upper ", qtxt, ")")),
		soil_ph   = c(paste0("Low pH (lower ",  qtxt, ")"),
									paste0("High pH (upper ", qtxt, ")"))
	)
	
	map_dfr(names(q_labels), function(env_var) {
		
		labels <- q_labels[[env_var]]
		strata <- stratify_within(boot_data, env_var, probs)
		
		if (nrow(strata$lower) < 50 || nrow(strata$upper) < 50) {
			warning(sprintf("Skipping %s: stratum too small (lower=%d, upper=%d)",
											env_var, nrow(strata$lower), nrow(strata$upper)))
			return(NULL)
		}
		
		lower_models <- map(traits, ~ {
			fit_rf_model(.x, strata$lower, covariates, hyper_grid,
									 num_threads = 1L)$trait_mod
		}) %>% set_names(traits)
		
		upper_models <- map(traits, ~ {
			fit_rf_model(.x, strata$upper, covariates, hyper_grid,
									 num_threads = 1L)$trait_mod
		}) %>% set_names(traits)
		
		map_dfr(traits, function(tr) {
			bind_rows(
				compute_pdp(lower_models[[tr]], strata$lower, "standage") %>%
					mutate(group = "low",  group_label = labels[1]),
				compute_pdp(upper_models[[tr]], strata$upper, "standage") %>%
					mutate(group = "high", group_label = labels[2])
			) %>%
				mutate(trait = tr, variable = env_var)
		})
	})
}


# ══════════════════════════════════════════════════════════════════════════════
# 3. Run bootstrap
# ══════════════════════════════════════════════════════════════════════════════

run_bootstrap <- function(data, leaf_type, hyper_grid, probs, n_boot, n_cores) {
	
	df_lt <- data %>% filter(leaf_type == !!leaf_type)
	message(sprintf("\n  %s: %d plots, %d bootstrap iterations",
									leaf_type, nrow(df_lt), n_boot))
	
	if (n_cores > 1) {
		cl <- makeCluster(n_cores)
		registerDoParallel(cl)
		message(sprintf("  Workers registered: %d", foreach::getDoParWorkers()))
		on.exit(stopCluster(cl), add = TRUE)
	} else {
		registerDoSEQ()
	}
	
	foreach(
		iter      = seq_len(n_boot),
		.combine  = bind_rows,
		.packages = c("ranger", "dplyr", "purrr", "tidyr"),
		.export   = c("run_pdp_once", "stratify_within", "fit_rf_model",
									"compute_pdp", "predict_fn", "COVARIATES", "TRAITS")
	) %dopar% {
		
		set.seed(iter)
		boot_data <- df_lt[sample(nrow(df_lt), replace = TRUE), ]
		
		tryCatch(
			run_pdp_once(boot_data, COVARIATES, TRAITS, hyper_grid, probs) %>%
				mutate(iteration = iter, leaf_type = leaf_type),
			error = function(e) {
				warning(sprintf("Boot iter %d failed: %s", iter, conditionMessage(e)))
				NULL
			}
		)
	}
}

message(sprintf("\nRunning PDP bootstrap (n = %d, probs = [%.2f, %.2f])...",
								N_BOOT, PROBS_MAIN[1], PROBS_MAIN[2]))

pdp_raw <- map_dfr(LEAF_TYPES, function(lt) {
	run_bootstrap(
		data       = data,
		leaf_type  = lt,
		hyper_grid = hyper_params[[lt]],
		probs      = PROBS_MAIN,
		n_boot     = N_BOOT,
		n_cores    = if (PARALLEL) N_CORES else 1L
	)
})

write_rds(pdp_raw, file.path(PATH_TABLES, "pdp_raw.rds"))
message(sprintf("PDP bootstrap complete: %d rows saved to tables/pdp_raw.rds",
								nrow(pdp_raw)))

if (RUN_SENS) {
	message(sprintf("\nRunning sensitivity analysis (probs = [%.2f, %.2f])...",
									PROBS_SENS[1], PROBS_SENS[2]))
	pdp_raw_sens <- map_dfr(LEAF_TYPES, function(lt) {
		run_bootstrap(data, lt, hyper_params[[lt]], PROBS_SENS, N_BOOT,
									if (PARALLEL) N_CORES else 1L)
	})
	write_rds(pdp_raw_sens, file.path(PATH_TABLES, "pdp_raw_sens.rds"))
	message("Sensitivity analysis saved to tables/pdp_raw_sens.rds")
}


# ══════════════════════════════════════════════════════════════════════════════
# 4. Extract slopes and intercepts per bootstrap iteration
# ══════════════════════════════════════════════════════════════════════════════
# Slope  = rate of trait change through succession in that environmental context
# Intercept = initial trait expression at the onset of succession
# Signed difference (high - low) preserves directionality:
#   positive slope_diff  → high-env group has steeper successional trajectory
#   negative slope_diff  → low-env group has steeper trajectory
#   positive intercept_diff → high-env group starts higher

pdp_stats <- pdp_raw %>%
	group_by(iteration, leaf_type, trait, variable, group) %>%
	summarise(
		slope     = coef(lm(yhat ~ standage))[2],
		intercept = coef(lm(yhat ~ standage))[1],
		.groups   = "drop"
	) %>%
	pivot_wider(
		names_from  = group,
		values_from = c(slope, intercept)
	) %>%
	mutate(
		slope_diff     = slope_high     - slope_low,
		intercept_diff = intercept_high - intercept_low,
		trait_label    = recode(trait,    !!!TRAIT_LABELS),
		variable_label = recode(variable, !!!ENV_LABELS)
	)

write_rds(pdp_stats, file.path(PATH_TABLES, "pdp_stats.rds"))


# ══════════════════════════════════════════════════════════════════════════════
# 5. SHAP-weighted mean summary
# ══════════════════════════════════════════════════════════════════════════════
# Weight each environmental variable's slope/intercept difference by its
# relative SHAP importance (from script 04), giving an importance-weighted
# overall measure of spatio-temporal interaction per trait × leaf type.

shap_weights <- shap_per_var %>%
	filter(variable %in% ENV_VARS) %>%
	group_by(trait, leaf_type) %>%
	mutate(weight = sum_abs_shap / sum(sum_abs_shap)) %>%
	ungroup() %>%
	dplyr::select(trait, leaf_type, variable, weight)

pdp_weighted <- pdp_stats %>%
	left_join(shap_weights, by = c("trait", "leaf_type", "variable")) %>%
	filter(!is.na(weight)) %>%
	group_by(iteration, leaf_type, trait, trait_label) %>%
	summarise(
		slope_diff     = sum(slope_diff     * weight, na.rm = TRUE),
		intercept_diff = sum(intercept_diff * weight, na.rm = TRUE),
		variable       = "SHAP-weighted mean",
		variable_label = "SHAP-Weighted Mean",
		.groups        = "drop"
	)

pdp_stats_full <- bind_rows(pdp_stats, pdp_weighted) %>%
	mutate(
		variable_label = factor(variable_label,
														levels = c("SHAP-Weighted Mean",
																			 "Temperature PC",
																			 "Soil water retention PC",
																			 "Precipitation PC",
																			 "Elevation",
																			 "Soil pH"))
	)


# ══════════════════════════════════════════════════════════════════════════════
# 6. Summary: medians, 95% CIs, and CI-excludes-zero criterion
# ══════════════════════════════════════════════════════════════════════════════
# Inference is based on whether the 95% bootstrap CI of the signed difference
# excludes zero. A difference is "robust" if the entire CI lies on one side of
# zero — meaning the direction of the effect is consistent across all resamples.
#
# This is preferred over formal significance testing (Wilcoxon, t-test) on
# bootstrapped quantities because:
#   - Formal tests on bootstrap distributions conflate resampling variance with
#     ecological effect size, producing inflated power
#   - The CI-excludes-zero criterion directly reflects what the bootstrap is
#     designed to show: stability of an effect across data resamples
#   - It maps cleanly onto result sentences: "X of 9 traits showed robust slope
#     differences under temperature stratification (95% CI excludes zero)"
#
# We additionally report the median absolute slope difference per variable
# (averaged across traits) as a summary of overall interaction strength.

pdp_summary <- pdp_stats %>%
	group_by(leaf_type, trait, trait_label, variable, variable_label) %>%
	summarise(
		# Medians and 95% bootstrap CIs for signed differences
		slope_median     = median(slope_diff,       na.rm = TRUE),
		slope_lwr        = quantile(slope_diff,     0.025, na.rm = TRUE),
		slope_upr        = quantile(slope_diff,     0.975, na.rm = TRUE),
		intercept_median = median(intercept_diff,   na.rm = TRUE),
		intercept_lwr    = quantile(intercept_diff, 0.025, na.rm = TRUE),
		intercept_upr    = quantile(intercept_diff, 0.975, na.rm = TRUE),
		
		# CI-excludes-zero: robust if CI lies entirely above or below zero
		# This is the primary inference criterion
		slope_robust     = (slope_lwr     > 0 | slope_upr     < 0),
		intercept_robust = (intercept_lwr > 0 | intercept_upr < 0),
		
		.groups = "drop"
	)

write_rds(pdp_summary, file.path(PATH_TABLES, "pdp_summary.rds"))


# ══════════════════════════════════════════════════════════════════════════════
# 7. Figures
# ══════════════════════════════════════════════════════════════════════════════

# Raw PDP curves (yhat vs standage) for high vs low groups across all bootstrap
# iterations. Points show iteration-level scatter; lines show linear fits.
# One panel per trait × variable, faceted by leaf type.

p_lines <- pdp_raw %>%
	mutate(
		leaf_type      = str_to_title(leaf_type),
		trait_label    = recode(trait,    !!!TRAIT_LABELS),
		variable_label = recode(variable, !!!ENV_LABELS)
	) %>%
	ggplot(aes(x = standage, y = yhat)) +
	geom_point(aes(fill = group), alpha = 0.04, size = 0.4,
						 shape = 21, stroke = 0.15, colour = "black") +
	geom_smooth(aes(colour = group), method = "lm",
							se = FALSE, linewidth = 0.6, formula = y ~ x) +
	ggh4x::facet_nested(
		rows = vars(leaf_type, trait_label),
		cols = vars(variable_label),
		scales = "free_y"
	) +
	scale_fill_manual(
		values = c("low" = "lightcyan", "high" = "lightcoral"),
		labels = c("low" = "Lower quantile", "high" = "Upper quantile"),
		name   = NULL
	) +
	scale_colour_manual(
		values = c("low" = "navy", "high" = "darkred"),
		labels = c("low" = "Lower quantile", "high" = "Upper quantile"),
		name   = NULL
	) +
	labs(x = "Stand age (years)", y = "Predicted trait expression") +
	theme_bw(base_size = 8) +
	theme(
		legend.position  = "top",
		strip.text       = element_text(size = 7, face = "bold"),
		strip.background = element_rect(fill = "white", colour = "black",
																		linewidth = 0.4),
		panel.grid.minor = element_blank()
	)

ggsave(
	file.path("figures/supplementary/pdp", "supp_pdp_lines.png"),
	plot   = p_lines,
	width  = 260, height = 360, units = "mm", dpi = 300
)
message("Supplementary PDP lines saved")


# ══════════════════════════════════════════════════════════════════════════════
# 8. Console summary
# ══════════════════════════════════════════════════════════════════════════════

message("\n── RQ2: Robust slope differences (95% CI excludes zero) ────────────")
pdp_summary %>%
	filter(slope_robust) %>%
	dplyr::select(leaf_type, trait_label, variable_label,
								slope_median, slope_lwr, slope_upr) %>%
	mutate(across(where(is.numeric), ~ round(., 5))) %>%
	arrange(leaf_type, variable_label, desc(abs(slope_median))) %>%
	print(n = Inf)

message("\n── Slope differences by variable (median |slope| across traits) ─────")
pdp_summary %>%
	group_by(leaf_type, variable_label) %>%
	summarise(
		median_abs_slope = median(abs(slope_median), na.rm = TRUE),
		n_robust_traits  = sum(slope_robust,         na.rm = TRUE),
		.groups = "drop"
	) %>%
	arrange(leaf_type, desc(median_abs_slope)) %>%
	mutate(across(where(is.numeric), ~ round(., 5))) %>%
	print(n = Inf)

message("\n── Intercept differences by variable (initial filtering strength) ───")
pdp_summary %>%
	group_by(leaf_type, variable_label) %>%
	summarise(
		median_abs_intercept = median(abs(intercept_median), na.rm = TRUE),
		n_robust_traits      = sum(intercept_robust,         na.rm = TRUE),
		.groups = "drop"
	) %>%
	arrange(leaf_type, desc(median_abs_intercept)) %>%
	mutate(across(where(is.numeric), ~ round(., 5))) %>%
	print(n = Inf)

message("\n── 05_pdp.R complete ───────────────────────────────────────────────")
message(sprintf("Tables: %s/", PATH_TABLES))
message(sprintf("Plots:  %s/", PATH_PLOTS))