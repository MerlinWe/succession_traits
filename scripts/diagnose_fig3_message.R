################################################################################
## Diagnostics for the conceptual redesign of Figure 3
##
## Run from the Trait Succession project root:
##   Rscript scripts/diagnose_fig3_message.R
##
## Outputs are stored under tables/diagnostics/fig3/.
################################################################################

rm(list = ls())

required_packages <- c("dplyr", "tidyr", "readr", "stringr")
missing_packages <- required_packages[
	!vapply(required_packages, requireNamespace, quietly = TRUE,
				FUN.VALUE = logical(1))
]
if (length(missing_packages) > 0L) {
	stop("Missing required package(s): ", paste(missing_packages, collapse = ", "),
			 call. = FALSE)
}

library(dplyr)
library(tidyr)
library(readr)
library(stringr)

input_files <- c(
	pdp_raw = "tables/pdp_raw.rds",
	pdp_stats = "tables/pdp_stats.rds",
	pdp_summary = "tables/pdp_summary.rds"
)
missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files) > 0L) {
	stop("Missing input file(s): ", paste(missing_files, collapse = ", "),
			 call. = FALSE)
}

output_dir <- "tables/diagnostics/fig3"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

pdp_raw <- read_rds(input_files[["pdp_raw"]])
pdp_stats <- read_rds(input_files[["pdp_stats"]])
pdp_summary <- read_rds(input_files[["pdp_summary"]])

LATE_AGE <- 100L

# Paired change in the magnitude of environmental differentiation.
gap_iteration <- pdp_raw %>%
	group_by(iteration, leaf_type, trait, variable, group) %>%
	summarise(
		yhat_early = yhat[which.min(abs(standage - 0))],
		yhat_late = yhat[which.min(abs(standage - LATE_AGE))],
		.groups = "drop"
	) %>%
	pivot_wider(names_from = group, values_from = c(yhat_early, yhat_late)) %>%
	mutate(
		gap_early = abs(yhat_early_high - yhat_early_low),
		gap_late = abs(yhat_late_high - yhat_late_low),
		gap_change = gap_late - gap_early
	)

gap_by_combination <- gap_iteration %>%
	group_by(leaf_type, trait, variable) %>%
	summarise(
		gap_early_median = median(gap_early, na.rm = TRUE),
		gap_late_median = median(gap_late, na.rm = TRUE),
		gap_change_median = median(gap_change, na.rm = TRUE),
		gap_change_lwr = quantile(gap_change, 0.025, na.rm = TRUE),
		gap_change_upr = quantile(gap_change, 0.975, na.rm = TRUE),
		gap_class = case_when(
			gap_change_lwr > 0 ~ "robust_widening",
			gap_change_upr < 0 ~ "robust_narrowing",
			TRUE ~ "uncertain_change"
		),
		.groups = "drop"
	)

gap_by_forest <- gap_by_combination %>%
	count(leaf_type, gap_class, name = "n_combinations") %>%
	complete(
		leaf_type,
		gap_class = c("robust_widening", "robust_narrowing", "uncertain_change"),
		fill = list(n_combinations = 0L)
	) %>%
	group_by(leaf_type) %>%
	mutate(proportion = n_combinations / sum(n_combinations)) %>%
	ungroup()

# Summarise slope-difference magnitude and consistency by environmental axis.
slope_consistency <- pdp_stats %>%
	group_by(leaf_type, trait, variable) %>%
	summarise(
		p_positive = mean(slope_diff > 0, na.rm = TRUE),
		directional_consistency = 2 * abs(p_positive - 0.5),
		slope_high_median = median(slope_high, na.rm = TRUE),
		slope_low_median = median(slope_low, na.rm = TRUE),
		abs_rate_difference = abs(slope_high_median) - abs(slope_low_median),
		.groups = "drop"
	)

slope_by_combination <- pdp_summary %>%
	left_join(
		slope_consistency,
		by = c("leaf_type", "trait", "variable")
	) %>%
	mutate(
		slope_difference_sign = case_when(
			slope_median > 0 ~ "more_positive_upper",
			slope_median < 0 ~ "more_positive_lower",
			TRUE ~ "no_difference"
		),
		absolute_rate_sign = case_when(
			abs_rate_difference > 0 ~ "faster_upper",
			abs_rate_difference < 0 ~ "faster_lower",
			TRUE ~ "equal_rate"
		),
		# TRUE where calling a positive difference "upper steeper" (or a negative
		# difference "lower steeper") gives the wrong absolute-rate conclusion.
		steeper_label_mismatch = case_when(
			slope_median > 0 ~ abs_rate_difference <= 0,
			slope_median < 0 ~ abs_rate_difference >= 0,
			TRUE ~ FALSE
		)
	)

slope_by_environment <- slope_by_combination %>%
	group_by(leaf_type, variable, variable_label) %>%
	summarise(
		n_traits = n(),
		n_robust = sum(slope_robust, na.rm = TRUE),
		median_abs_delta_slope = median(abs(slope_median), na.rm = TRUE),
		median_directional_consistency = median(directional_consistency, na.rm = TRUE),
		n_more_positive_upper = sum(slope_median > 0, na.rm = TRUE),
		n_more_positive_lower = sum(slope_median < 0, na.rm = TRUE),
		n_steeper_label_mismatch = sum(steeper_label_mismatch, na.rm = TRUE),
		.groups = "drop"
	) %>%
	arrange(leaf_type, desc(median_abs_delta_slope))

temperature_detail <- slope_by_combination %>%
	filter(variable == "temp_pc") %>%
	select(
		leaf_type, trait, trait_label,
		slope_high_median, slope_low_median, slope_median,
		abs_rate_difference, slope_lwr, slope_upr, slope_robust,
		directional_consistency, slope_difference_sign, absolute_rate_sign,
		steeper_label_mismatch
	) %>%
	arrange(leaf_type, desc(slope_median))

write_csv(gap_by_combination, file.path(output_dir, "gap_change_by_combination.csv"))
write_csv(gap_by_forest, file.path(output_dir, "gap_change_summary.csv"))
write_csv(slope_by_environment, file.path(output_dir, "slope_summary_by_environment.csv"))
write_csv(temperature_detail, file.path(output_dir, "temperature_slope_detail.csv"))

message("\nPaired gap-change classifications:")
print(gap_by_forest, n = Inf)

message("\nSlope-modification summary by environmental axis:")
print(slope_by_environment, n = Inf)

message("\nPotentially misleading 'steeper' classifications:")
print(
	slope_by_combination %>%
		summarise(
			n_total = n(),
			n_mismatch = sum(steeper_label_mismatch, na.rm = TRUE),
			proportion_mismatch = mean(steeper_label_mismatch, na.rm = TRUE)
		)
)

message("\nDiagnostics saved under ", output_dir, "/")
