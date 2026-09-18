################################################################################
## succession_traits: post-run audit for grouped VEcv outputs
##
## Validates the completed production run of scripts/06_vecv.R, reconstructs
## combined objects from independent checkpoints, and writes manuscript-facing
## summaries. It does not refit models or modify production outputs.
##
## Run from the project root:
##   Rscript scripts/audit_grouped_vecv.R
################################################################################

rm(list = ls())

if (!requireNamespace("tidyverse", quietly = TRUE)) {
	stop("Missing package: tidyverse", call. = FALSE)
}
suppressPackageStartupMessages(library(tidyverse))

RUN_DIR <- "tables/vecv_grouped_pid/normal_pid_f10_r30_seed42"
CHECKPOINT_DIR <- file.path(RUN_DIR, "checkpoints")
OUT_DIR <- "tables/diagnostics/submission_audit_grouped_vecv"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

TRAITS <- c(
	"bark_thickness", "conduit_diam", "height", "leaf_density", "leaf_k",
	"root_depth", "seed_dry_mass", "shade_tolerance", "specific_leaf_area"
)
LEAF_TYPES <- c("broadleaf", "coniferous")
N_REPEATS <- 30L
N_FOLDS <- 10L
MIN_BIN_N <- 30L
BASE_SEED <- 42L
COMMON_AGE_MIN <- 15
COMMON_AGE_MAX <- 125
EARLY_AGES <- c(15, 25)
LATE_AGES <- c(105, 115, 125)

required_files <- c(
	file.path(RUN_DIR, "run_metadata.rds"),
	file.path(RUN_DIR, "vecv_raw_unfiltered.rds"),
	file.path(RUN_DIR, "vecv_raw.rds"),
	file.path(RUN_DIR, "vecv_summary.rds"),
	file.path(RUN_DIR, "vecv_divergence_raw.rds"),
	file.path(RUN_DIR, "vecv_divergence.rds")
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
	stop("Required grouped-VEcv file(s) missing:\n  ",
		paste(missing_files, collapse = "\n  "), call. = FALSE)
}

checks <- tibble(
	section = character(), check = character(), status = character(),
	detail = character()
)
add_check <- function(section, check, status, detail) {
	checks <<- bind_rows(
		checks,
		tibble(section = section, check = check, status = status, detail = detail)
	)
}
same_data <- function(x, y, keys) {
	x <- arrange(x, across(all_of(keys)))
	y <- arrange(y, across(all_of(keys)))
	isTRUE(all.equal(x, y, check.attributes = TRUE))
}
max_abs_diff <- function(x, y) {
	z <- abs(x - y)
	if (!any(is.finite(z))) return(NA_real_)
	max(z, na.rm = TRUE)
}
q025 <- function(x) unname(quantile(x, 0.025, na.rm = TRUE, type = 8))
q975 <- function(x) unname(quantile(x, 0.975, na.rm = TRUE, type = 8))

# ── Metadata and checkpoint completeness ─────────────────────────────────────

metadata <- read_rds(file.path(RUN_DIR, "run_metadata.rds"))
metadata_ok <- identical(metadata$analysis_version, "grouped-pid-v1") &&
	identical(metadata$config_id, "normal_pid_f10_r30_seed42") &&
	identical(metadata$resampling_unit, "PID") &&
	identical(metadata$n_folds, N_FOLDS) &&
	identical(metadata$n_repeats, N_REPEATS) &&
	identical(metadata$minimum_bin_n, MIN_BIN_N) &&
	identical(metadata$base_seed, BASE_SEED) &&
	identical(metadata$smoke_test, FALSE)
add_check(
	"metadata", "Production configuration matches the grouped run",
	if_else(metadata_ok, "PASS", "FAIL"),
	sprintf("config=%s; unit=%s; folds=%d; repeats=%d; minimum n=%d",
		metadata$config_id, metadata$resampling_unit, metadata$n_folds,
		metadata$n_repeats, metadata$minimum_bin_n)
)

expected_jobs <- expand_grid(
	trait = TRAITS,
	leaf_type = LEAF_TYPES,
	repeat_id = seq_len(N_REPEATS)
) %>%
	mutate(
		job_id = sprintf("%s__%s__rep%03d", leaf_type, trait, repeat_id),
		expected_seed = as.integer(
			BASE_SEED + 1000000L * match(trait, TRAITS) +
				10000L * match(leaf_type, LEAF_TYPES) + 100L * repeat_id
		),
		path = file.path(CHECKPOINT_DIR, paste0(job_id, ".rds"))
	)

checkpoint_files <- sort(list.files(
	CHECKPOINT_DIR, pattern = "\\.rds$", full.names = TRUE
))
payloads <- lapply(expected_jobs$path, function(path) {
	if (!file.exists(path)) return(NULL)
	tryCatch(read_rds(path), error = function(e) NULL)
})
valid_payload <- map2_lgl(payloads, seq_len(nrow(expected_jobs)), function(x, i) {
	!is.null(x) && is.list(x) && identical(x$status, "complete") &&
		identical(x$analysis_version, metadata$analysis_version) &&
		identical(x$config_id, metadata$config_id) &&
		identical(x$input_signature, metadata$input_signature) &&
		identical(x$settings_signature, metadata$settings_signature) &&
		identical(x$job_id, expected_jobs$job_id[i]) &&
		identical(x$seed, expected_jobs$expected_seed[i]) &&
		is.data.frame(x$result) && nrow(x$result) == 150L
})
checkpoint_ok <- length(checkpoint_files) == nrow(expected_jobs) &&
	all(valid_payload)
add_check(
	"checkpoints", "All expected independent jobs are present and valid",
	if_else(checkpoint_ok, "PASS", "FAIL"),
	sprintf("%d/%d valid checkpoints (18 trait-forest cells x 30 repeats)",
		sum(valid_payload), nrow(expected_jobs))
)
if (!all(valid_payload)) {
	write_csv(
		expected_jobs %>% mutate(valid = valid_payload) %>% filter(!valid),
		file.path(OUT_DIR, "invalid_or_missing_checkpoints.csv")
	)
	stop("Checkpoint audit failed; see diagnostic output.", call. = FALSE)
}

error_files <- if (dir.exists(file.path(RUN_DIR, "errors"))) {
	list.files(file.path(RUN_DIR, "errors"), full.names = TRUE)
} else {
	character()
}
add_check(
	"checkpoints", "No worker error payloads remain",
	if_else(length(error_files) == 0L, "PASS", "WARNING"),
	sprintf("%d error payload(s)", length(error_files))
)

# ── Reconstruction and saved summary integrity ───────────────────────────────

vecv_raw_unfiltered <- read_rds(file.path(RUN_DIR, "vecv_raw_unfiltered.rds"))
vecv_raw <- read_rds(file.path(RUN_DIR, "vecv_raw.rds"))
vecv_summary <- read_rds(file.path(RUN_DIR, "vecv_summary.rds"))
vecv_divergence_saved <- read_rds(file.path(RUN_DIR, "vecv_divergence.rds"))

reconstructed <- map_dfr(payloads, "result")
raw_keys <- c("job_id", "variable", "env_group", "standage_mid")
reconstruction_ok <- same_data(reconstructed, vecv_raw_unfiltered, raw_keys)
add_check(
	"raw", "Unfiltered raw table reconstructs exactly from checkpoints",
	if_else(reconstruction_ok, "PASS", "FAIL"),
	sprintf("%s rows reconstructed; expected 81,000",
		format(nrow(reconstructed), big.mark = ","))
)

filtered_recalc <- filter(vecv_raw_unfiltered, n >= MIN_BIN_N)
filter_ok <- same_data(filtered_recalc, vecv_raw, raw_keys)
add_check(
	"raw", "Minimum-bin-size filtering reproduces exactly",
	if_else(filter_ok, "PASS", "FAIL"),
	sprintf("%s/%s rows retained; n range %d-%d; group-count range %d-%d",
		format(nrow(vecv_raw), big.mark = ","),
		format(nrow(vecv_raw_unfiltered), big.mark = ","),
		min(vecv_raw$n), max(vecv_raw$n),
		min(vecv_raw$n_groups), max(vecv_raw$n_groups))
)

finite_ok <- all(is.finite(vecv_raw$VEcv)) && all(is.finite(vecv_raw$E1)) &&
	all(vecv_raw$n_groups <= vecv_raw$n)
add_check(
	"raw", "Retained skill estimates and sample counts are valid",
	if_else(finite_ok, "PASS", "FAIL"),
	sprintf("VEcv range %.3f to %.3f; E1 range %.3f to %.3f",
		min(vecv_raw$VEcv), max(vecv_raw$VEcv),
		min(vecv_raw$E1), max(vecv_raw$E1))
)

summary_recalc <- vecv_raw %>%
	group_by(
		trait, trait_label, leaf_type, variable, variable_label,
		env_group, standage_bin, standage_mid
	) %>%
	summarise(
		n_med = median(n),
		n_groups_med = median(n_groups),
		VEcv_med = median(VEcv),
		VEcv_lwr = quantile(VEcv, 0.025),
		VEcv_upr = quantile(VEcv, 0.975),
		E1_med = median(E1),
		E1_lwr = quantile(E1, 0.025),
		E1_upr = quantile(E1, 0.975),
		n_repeats = n_distinct(repeat_id),
		.groups = "drop"
	)
summary_ok <- same_data(
	summary_recalc, vecv_summary,
	c("trait", "leaf_type", "variable", "env_group", "standage_mid")
)
add_check(
	"summary", "Saved trait-stratum summaries reproduce from raw rows",
	if_else(summary_ok, "PASS", "FAIL"),
	sprintf("%d cells; every retained cell has %d repeats",
		nrow(summary_recalc), min(summary_recalc$n_repeats))
)

paired_raw <- vecv_raw %>%
	select(
		trait, trait_label, leaf_type, variable, variable_label,
		standage_bin, standage_mid, repeat_id, env_group, VEcv, E1
	) %>%
	pivot_wider(names_from = env_group, values_from = c(VEcv, E1)) %>%
	mutate(
		delta_VEcv = VEcv_high - VEcv_low,
		abs_delta_VEcv = abs(delta_VEcv),
		delta_E1 = E1_high - E1_low,
		abs_delta_E1 = abs(delta_E1)
	)

divergence_complete <- paired_raw %>%
	filter(is.finite(delta_VEcv)) %>%
	group_by(
		trait, trait_label, leaf_type, variable, variable_label,
		standage_bin, standage_mid
	) %>%
	summarise(
		delta_med = median(delta_VEcv),
		delta_lwr = quantile(delta_VEcv, 0.025),
		delta_upr = quantile(delta_VEcv, 0.975),
		abs_delta_med = median(abs_delta_VEcv),
		abs_delta_lwr = quantile(abs_delta_VEcv, 0.025),
		abs_delta_upr = quantile(abs_delta_VEcv, 0.975),
		n_paired_repeats = n_distinct(repeat_id),
		.groups = "drop"
	) %>%
	mutate(
		direction_stable = delta_lwr > 0 | delta_upr < 0,
		sig_divergence = direction_stable
	)

divergence_compare <- inner_join(
	divergence_complete,
	vecv_divergence_saved,
	by = c(
		"trait", "trait_label", "leaf_type", "variable", "variable_label",
		"standage_bin", "standage_mid"
	),
	suffix = c("_calc", "_saved")
)
divergence_max_diff <- max(
	max_abs_diff(divergence_compare$delta_med_calc,
		divergence_compare$delta_med_saved),
	max_abs_diff(divergence_compare$delta_lwr_calc,
		divergence_compare$delta_lwr_saved),
	max_abs_diff(divergence_compare$delta_upr_calc,
		divergence_compare$delta_upr_saved),
	max_abs_diff(divergence_compare$abs_delta_med_calc,
		divergence_compare$abs_delta_med_saved),
	na.rm = TRUE
)
add_check(
	"divergence", "Finite high-low contrasts reproduce from paired repeats",
	if_else(divergence_max_diff < 1e-10, "PASS", "FAIL"),
	sprintf("%d complete contrast cells; maximum discrepancy %.3g",
		nrow(divergence_complete), divergence_max_diff)
)

incomplete_pairs <- paired_raw %>%
	group_by(trait, leaf_type, variable, standage_mid) %>%
	summarise(n_paired_repeats = sum(is.finite(delta_VEcv)), .groups = "drop") %>%
	filter(n_paired_repeats < N_REPEATS)
write_csv(incomplete_pairs,
	file.path(OUT_DIR, "incomplete_high_low_pairs.csv"))
add_check(
	"divergence", "Every divergence row has a complete high-low pair",
	if_else(nrow(incomplete_pairs) == 0L, "PASS", "WARNING"),
	sprintf(
		"%d edge-age cells lack a pair and are omitted from audited contrasts",
		nrow(incomplete_pairs)
	)
)
write_rds(divergence_complete,
	file.path(OUT_DIR, "vecv_divergence_complete_pairs.rds"))
write_csv(divergence_complete,
	file.path(OUT_DIR, "vecv_divergence_complete_pairs.csv"))

# ── Manuscript-facing summaries on the complete common age window ────────────

raw_common <- vecv_raw %>%
	filter(between(standage_mid, COMMON_AGE_MIN, COMMON_AGE_MAX))

# Equal-weight traits after averaging their ten environmental strata. Repeat
# intervals quantify partition stability, not population sampling uncertainty.
predictability_by_repeat <- raw_common %>%
	group_by(leaf_type, trait, repeat_id, standage_mid) %>%
	summarise(
		VEcv = mean(VEcv), E1 = mean(E1), n_context_cells = n(),
		.groups = "drop"
	) %>%
	filter(n_context_cells == 10L) %>%
	group_by(leaf_type, repeat_id, standage_mid) %>%
	summarise(
		VEcv = mean(VEcv), E1 = mean(E1),
		n_traits = n_distinct(trait), .groups = "drop"
	) %>%
	filter(n_traits == 9L)

predictability_age <- predictability_by_repeat %>%
	group_by(leaf_type, standage_mid) %>%
	summarise(
		VEcv_med = median(VEcv), VEcv_lwr = q025(VEcv), VEcv_upr = q975(VEcv),
		E1_med = median(E1), E1_lwr = q025(E1), E1_upr = q975(E1),
		.groups = "drop"
	)
write_csv(predictability_age,
	file.path(OUT_DIR, "predictability_by_age.csv"))

predictability_stage_repeat <- predictability_by_repeat %>%
	mutate(stage = case_when(
		standage_mid %in% EARLY_AGES ~ "early",
		standage_mid %in% LATE_AGES ~ "late",
		TRUE ~ NA_character_
	)) %>%
	filter(!is.na(stage)) %>%
	group_by(leaf_type, repeat_id, stage) %>%
	summarise(VEcv = mean(VEcv), E1 = mean(E1), .groups = "drop") %>%
	pivot_wider(names_from = stage, values_from = c(VEcv, E1)) %>%
	mutate(
		VEcv_change = VEcv_late - VEcv_early,
		E1_change = E1_late - E1_early
	)

predictability_stage <- predictability_stage_repeat %>%
	group_by(leaf_type) %>%
	summarise(
		VEcv_early = median(VEcv_early), VEcv_late = median(VEcv_late),
		VEcv_change_med = median(VEcv_change),
		VEcv_change_lwr = q025(VEcv_change),
		VEcv_change_upr = q975(VEcv_change),
		E1_early = median(E1_early), E1_late = median(E1_late),
		E1_change_med = median(E1_change),
		E1_change_lwr = q025(E1_change),
		E1_change_upr = q975(E1_change),
		.groups = "drop"
	)
write_csv(predictability_stage,
	file.path(OUT_DIR, "predictability_early_late.csv"))

predictability_trends <- predictability_by_repeat %>%
	group_by(leaf_type, repeat_id) %>%
	summarise(
		VEcv_slope = unname(coef(lm(VEcv ~ standage_mid))[2]),
		E1_slope = unname(coef(lm(E1 ~ standage_mid))[2]),
		.groups = "drop"
	) %>%
	group_by(leaf_type) %>%
	summarise(
		VEcv_slope_med = median(VEcv_slope),
		VEcv_slope_lwr = q025(VEcv_slope),
		VEcv_slope_upr = q975(VEcv_slope),
		VEcv_fraction_positive = mean(VEcv_slope > 0),
		E1_slope_med = median(E1_slope),
		E1_slope_lwr = q025(E1_slope),
		E1_slope_upr = q975(E1_slope),
		E1_fraction_positive = mean(E1_slope > 0),
		.groups = "drop"
	)
write_csv(predictability_trends,
	file.path(OUT_DIR, "predictability_trend_stability.csv"))

paired_common <- paired_raw %>%
	filter(
		between(standage_mid, COMMON_AGE_MIN, COMMON_AGE_MAX),
		is.finite(abs_delta_VEcv), is.finite(abs_delta_E1)
	)

dependence_by_repeat <- paired_common %>%
	group_by(leaf_type, variable, variable_label, repeat_id, standage_mid) %>%
	summarise(
		abs_delta_VEcv = mean(abs_delta_VEcv),
		abs_delta_E1 = mean(abs_delta_E1),
		n_traits = n_distinct(trait), .groups = "drop"
	) %>%
	filter(n_traits == 9L)

dependence_age <- dependence_by_repeat %>%
	group_by(leaf_type, variable, variable_label, standage_mid) %>%
	summarise(
		abs_delta_VEcv_med = median(abs_delta_VEcv),
		abs_delta_VEcv_lwr = q025(abs_delta_VEcv),
		abs_delta_VEcv_upr = q975(abs_delta_VEcv),
		abs_delta_E1_med = median(abs_delta_E1),
		abs_delta_E1_lwr = q025(abs_delta_E1),
		abs_delta_E1_upr = q975(abs_delta_E1),
		.groups = "drop"
	)
write_csv(dependence_age,
	file.path(OUT_DIR, "environmental_dependence_by_age.csv"))

dependence_stage_repeat <- dependence_by_repeat %>%
	mutate(stage = case_when(
		standage_mid %in% EARLY_AGES ~ "early",
		standage_mid %in% LATE_AGES ~ "late",
		TRUE ~ NA_character_
	)) %>%
	filter(!is.na(stage)) %>%
	group_by(leaf_type, variable, variable_label, repeat_id, stage) %>%
	summarise(
		abs_delta_VEcv = mean(abs_delta_VEcv),
		abs_delta_E1 = mean(abs_delta_E1), .groups = "drop"
	) %>%
	pivot_wider(
		names_from = stage,
		values_from = c(abs_delta_VEcv, abs_delta_E1)
	) %>%
	mutate(
		VEcv_change = abs_delta_VEcv_late - abs_delta_VEcv_early,
		E1_change = abs_delta_E1_late - abs_delta_E1_early
	)

dependence_stage <- dependence_stage_repeat %>%
	group_by(leaf_type, variable, variable_label) %>%
	summarise(
		VEcv_early = median(abs_delta_VEcv_early),
		VEcv_late = median(abs_delta_VEcv_late),
		VEcv_change_med = median(VEcv_change),
		VEcv_change_lwr = q025(VEcv_change),
		VEcv_change_upr = q975(VEcv_change),
		E1_early = median(abs_delta_E1_early),
		E1_late = median(abs_delta_E1_late),
		E1_change_med = median(E1_change),
		E1_change_lwr = q025(E1_change),
		E1_change_upr = q975(E1_change),
		.groups = "drop"
	)
write_csv(dependence_stage,
	file.path(OUT_DIR, "environmental_dependence_early_late.csv"))

dependence_trends <- dependence_by_repeat %>%
	group_by(leaf_type, variable, variable_label, repeat_id) %>%
	summarise(
		VEcv_slope = unname(coef(lm(abs_delta_VEcv ~ standage_mid))[2]),
		E1_slope = unname(coef(lm(abs_delta_E1 ~ standage_mid))[2]),
		.groups = "drop"
	) %>%
	group_by(leaf_type, variable, variable_label) %>%
	summarise(
		VEcv_slope_med = median(VEcv_slope),
		VEcv_slope_lwr = q025(VEcv_slope),
		VEcv_slope_upr = q975(VEcv_slope),
		VEcv_fraction_increasing = mean(VEcv_slope > 0),
		E1_slope_med = median(E1_slope),
		E1_slope_lwr = q025(E1_slope),
		E1_slope_upr = q975(E1_slope),
		E1_fraction_increasing = mean(E1_slope > 0),
		.groups = "drop"
	)
write_csv(dependence_trends,
	file.path(OUT_DIR, "environmental_dependence_trend_stability.csv"))

# Persistence cannot be established merely because an absolute gap is positive.
# Use the signed contrasts to quantify how often the repeated-partition interval
# excludes zero in the three late-successional bins.
late_signed_support <- divergence_complete %>%
	filter(standage_mid %in% LATE_AGES) %>%
	group_by(leaf_type, variable, variable_label) %>%
	summarise(
		n_trait_age_cells = n(),
		n_direction_stable = sum(direction_stable),
		proportion_direction_stable = mean(direction_stable),
		median_abs_delta_VEcv = median(abs_delta_med),
		.groups = "drop"
	)
write_csv(late_signed_support,
	file.path(OUT_DIR, "late_signed_contrast_support.csv"))

late_signed_by_forest <- late_signed_support %>%
	group_by(leaf_type) %>%
	summarise(
		n_trait_age_cells = sum(n_trait_age_cells),
		n_direction_stable = sum(n_direction_stable),
		proportion_direction_stable =
			n_direction_stable / n_trait_age_cells,
		.groups = "drop"
	)
write_csv(late_signed_by_forest,
	file.path(OUT_DIR, "late_signed_contrast_support_by_forest.csv"))

add_check(
	"narrative", "Mean grouped predictive skill is higher in later succession",
	if_else(
		all(predictability_stage$VEcv_change_lwr > 0) &&
			all(predictability_stage$E1_change_lwr > 0),
		"PASS", "FAIL"
	),
	paste0(
		paste(
			predictability_stage$leaf_type,
			sprintf("VEcv %.3f to %.3f",
				predictability_stage$VEcv_early,
				predictability_stage$VEcv_late),
			collapse = "; "
		),
		"; E1 agrees"
	)
)

add_check(
	"narrative", "Signed environmental predictability contrasts persist late",
	if_else(
		all(late_signed_support$proportion_direction_stable >= 0.75),
		"PASS", "WARNING"
	),
	paste0(
		paste(
			late_signed_by_forest$leaf_type,
			sprintf("%d/%d (%.1f%%) stable trait-gradient-age cells",
				late_signed_by_forest$n_direction_stable,
				late_signed_by_forest$n_trait_age_cells,
				100 * late_signed_by_forest$proportion_direction_stable),
			collapse = "; "
		),
		sprintf(
			"; late mean |delta VEcv| range %.3f-%.3f",
			min(dependence_stage$VEcv_late), max(dependence_stage$VEcv_late)
		)
	)
)

temp_conifer <- dependence_stage %>%
	filter(leaf_type == "coniferous", variable == "temp_pc")
add_check(
	"narrative", "Legacy claim that temperature divergence uniquely increases",
	"FAIL",
	sprintf(
		"Grouped coniferous temperature |delta VEcv| changes from %.3f to %.3f (late minus early %.3f); directions are heterogeneous",
		temp_conifer$VEcv_early, temp_conifer$VEcv_late,
		temp_conifer$VEcv_change_med
	)
)

# ── Legacy comparison and claim ledger ───────────────────────────────────────

legacy_path <- "tables/vecv_summary.rds"
if (file.exists(legacy_path)) {
	legacy_summary <- read_rds(legacy_path)
	legacy_vs_grouped <- inner_join(
		legacy_summary %>%
			select(
				trait, leaf_type, variable, env_group, standage_mid,
				VEcv_legacy = VEcv_med
			),
		vecv_summary %>%
			select(
				trait, leaf_type, variable, env_group, standage_mid,
				VEcv_grouped = VEcv_med
			),
		by = c("trait", "leaf_type", "variable", "env_group", "standage_mid")
	) %>%
	filter(is.finite(VEcv_legacy), is.finite(VEcv_grouped)) %>%
	group_by(leaf_type) %>%
	summarise(
		n_cells = n(),
		legacy_mean = mean(VEcv_legacy),
		grouped_mean = mean(VEcv_grouped),
		mean_shift = mean(VEcv_grouped - VEcv_legacy),
		cellwise_correlation = cor(VEcv_legacy, VEcv_grouped),
		.groups = "drop"
	)
	write_csv(legacy_vs_grouped,
		file.path(OUT_DIR, "legacy_vs_grouped_predictability.csv"))
}

claim_ledger <- tribble(
	~claim, ~audit_status, ~grouped_evidence, ~submission_action,
	"Predictability is higher in later than early succession",
	"SUPPORTED",
	sprintf(
		"Early-to-late VEcv change is %.3f in broadleaf and %.3f in coniferous forests; repeated-partition intervals exclude zero and E1 agrees.",
		predictability_stage$VEcv_change_med[predictability_stage$leaf_type == "broadleaf"],
		predictability_stage$VEcv_change_med[predictability_stage$leaf_type == "coniferous"]
	),
	"Retain this qualified claim, replace legacy point estimates with grouped estimates, and define intervals as partition-stability intervals; do not imply a monotonic rise.",
	"Predictability remains contingent on environmental context in later succession",
	"SUPPORTED",
	paste0(
		"Signed repeated-partition intervals exclude zero in ",
		sprintf(
			"%.1f%% of broadleaf and %.1f%% of coniferous late trait-gradient-age cells; late mean absolute gaps span %.3f-%.3f.",
			100 * late_signed_by_forest$proportion_direction_stable[
				late_signed_by_forest$leaf_type == "broadleaf"
			],
			100 * late_signed_by_forest$proportion_direction_stable[
				late_signed_by_forest$leaf_type == "coniferous"
			],
			min(dependence_stage$VEcv_late), max(dependence_stage$VEcv_late)
		)
	),
	"Retain as the principal Figure 4 message; use absolute gaps in the main figure and signed contrasts in the supplement.",
	"Environmental predictability gaps generally narrow, except for temperature",
	"NOT SUPPORTED",
	"Grouped resampling shows both strengthening and weakening; broadleaf temperature strengthens, whereas coniferous temperature weakens modestly.",
	"Delete the single-exception narrative from the Abstract, Results, Discussion, and implications paragraph.",
	"Coniferous forests are less predictable than broadleaf forests",
	"NOT SUPPORTED BY GROUPED CV",
	sprintf(
		"Grouped early-to-late averages are %.3f to %.3f for broadleaf and %.3f to %.3f for coniferous forests.",
		predictability_stage$VEcv_early[predictability_stage$leaf_type == "broadleaf"],
		predictability_stage$VEcv_late[predictability_stage$leaf_type == "broadleaf"],
		predictability_stage$VEcv_early[predictability_stage$leaf_type == "coniferous"],
		predictability_stage$VEcv_late[predictability_stage$leaf_type == "coniferous"]
	),
	"Remove or reverse the legacy comparison; avoid causal interpretation of this descriptive difference."
)
write_csv(claim_ledger, file.path(OUT_DIR, "claim_ledger.csv"))

write_csv(checks, file.path(OUT_DIR, "audit_checks.csv"))

message("Grouped VEcv audit complete.")
message("  Production outputs were not modified.")
message("  Diagnostics: ", OUT_DIR)
print(checks, n = Inf)
