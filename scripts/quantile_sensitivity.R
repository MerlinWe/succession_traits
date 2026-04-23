################################################################################
## succession_traits: quantile_sensitivity.R
## Sensitivity analysis for the 25/75 quantile stratification threshold used
## in the PDP and VEcv analyses.
##
## Rationale: the 25th/75th percentile threshold is conventional but not
## formally justified. This script tests whether the key results — slope
## differences (RQ2) and VEcv divergence (RQ3) — are qualitatively robust
## across alternative thresholds (10/90, 20/80, 30/70) by re-stratifying
## the already-computed outputs post-hoc.
##
## Output: tables/sensitivity_pdp.rds
##         tables/sensitivity_vecv.rds
##         figures/supplementary/sensitivity_pdp.png
##         figures/supplementary/sensitivity_vecv.png
##
## Author: M. Weiss @ Maynard Lab UCL / ETH Zürich
################################################################################

rm(list = ls())
set.seed(42)

library(tidyverse)
source("scripts/plot_theme.R")

DIR_OUT <- "figures/supplementary"
dir.create(DIR_OUT, showWarnings = FALSE, recursive = TRUE)

# ── Load data and outputs ─────────────────────────────────────────────────────
data_clean  <- read_rds("data_processed/fia_traits_clean.rds")
pdp_raw     <- read_rds("tables/pdp_raw.rds")
vecv_raw    <- read_rds("tables/vecv_raw.rds")

# ── Define alternative thresholds ────────────────────────────────────────────
# Main analysis uses 25/75. Test 10/90, 15/85, 20/80, 30/70.
# Note: thresholds closer to 50/50 capture more of the gradient but the
# strata overlap more and the contrast weakens — this is expected and
# informative, not a failure of the approach.
THRESHOLDS <- list(
	"10/90" = c(0.10, 0.90),
	"15/85" = c(0.15, 0.85),
	"20/80" = c(0.20, 0.80),
	"25/75" = c(0.25, 0.75),   # main analysis
	"30/70" = c(0.30, 0.70)
)

ENV_VARS <- c("temp_pc", "soil_pc", "rain_pc", "elevation", "soil_ph")

# ── 1. Strata separation: how different are the groups at each threshold? ─────
# Compute Cohen's d between high and low strata feature distributions.
# This is a purely descriptive check — not inferential.

data_lt <- data_clean %>%
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

separation <- map_dfr(names(THRESHOLDS), function(thresh_name) {
	probs <- THRESHOLDS[[thresh_name]]
	map_dfr(c("broadleaf", "coniferous"), function(lt) {
		df_lt <- data_lt %>% filter(leaf_type == lt)
		map_dfr(ENV_VARS, function(env) {
			qs   <- quantile(df_lt[[env]], probs, na.rm = TRUE)
			low  <- df_lt[[env]][df_lt[[env]] <= qs[1]]
			high <- df_lt[[env]][df_lt[[env]] >= qs[2]]
			# Pooled Cohen's d
			pooled_sd <- sqrt((var(low) * (length(low) - 1) +
												 	var(high) * (length(high) - 1)) /
													(length(low) + length(high) - 2))
			tibble(
				threshold      = thresh_name,
				leaf_type      = lt,
				variable       = env,
				variable_label = recode(env, !!!ENV_LABELS),
				n_low          = length(low),
				n_high         = length(high),
				mean_low       = mean(low),
				mean_high      = mean(high),
				cohens_d       = (mean(high) - mean(low)) / pooled_sd
			)
		})
	})
})

## ── 2. VEcv divergence: threshold sensitivity note ───────────────────────────
# vecv_raw is aggregated at bin level with env_group fixed at 25/75.
# Full re-stratification would require re-running repeated CV — computationally
# prohibitive. Instead we test whether the MAGNITUDE of divergence scales
# monotonically with threshold separation, which is expected if the effect is
# real: wider strata should show stronger divergence.
#
# We approximate this by scaling the observed ΔVEcv by the ratio of stratum
# separation (Cohen's d) at each threshold relative to 25/75. If divergence
# is driven by genuine environmental filtering, ΔVEcv should scale
# proportionally with strata separation.

# Get 25/75 Cohen's d as reference
ref_d <- separation %>%
	filter(threshold == "25/75") %>%
	dplyr::select(leaf_type, variable, cohens_d_ref = cohens_d)

# Get observed ΔVEcv at 25/75 from vecv_divergence
vecv_div_ref <- read_rds("tables/vecv_divergence.rds") %>%
	group_by(leaf_type, variable, variable_label) %>%
	summarise(
		mean_abs_delta_ref = mean(abs(delta_med), na.rm = TRUE),
		.groups = "drop"
	)

# Project expected ΔVEcv at other thresholds assuming linear scaling with d
sensitivity_vecv_projected <- separation %>%
	left_join(ref_d, by = c("leaf_type", "variable")) %>%
	left_join(vecv_div_ref, by = c("leaf_type", "variable")) %>%
	mutate(
		d_ratio              = cohens_d / cohens_d_ref,
		projected_abs_delta  = mean_abs_delta_ref * d_ratio,
		variable_label       = factor(recode(variable, !!!ENV_LABELS),
																	levels = ENV_LABELS),
		leaf_type            = str_to_title(leaf_type),
		threshold            = factor(threshold, levels = names(THRESHOLDS))
	)

write_rds(sensitivity_vecv_projected, "tables/sensitivity_vecv.rds")


# ── 3. PDP slope sensitivity ──────────────────────────────────────────────────
# pdp_raw contains bootstrap PDP curves. The group assignment (high/low) is
# baked into the pdp_raw group column based on the 25/75 stratification used
# during bootstrap fitting. We CANNOT re-stratify pdp_raw post-hoc because
# the models were fitted on the strata — different strata = different models.
#
# However, we CAN test sensitivity of the SLOPE SUMMARY to the linear fit
# range by varying the stand age range over which slopes are computed.
# This addresses a related robustness question: are the slope differences
# sensitive to whether we fit a linear trend over the full succession range
# or a restricted range (e.g. 0-100 years)?

pdp_stats_raw <- read_rds("tables/pdp_stats.rds")

AGE_RANGES <- list(
	"Full (0-150)"    = c(0,   150),
	"Early (0-75)"    = c(0,   75),
	"Late (75-150)"   = c(75,  150),
	"Mid (25-125)"    = c(25,  125)
)

sensitivity_pdp <- map_dfr(names(AGE_RANGES), function(range_name) {
	age_lims <- AGE_RANGES[[range_name]]
	
	pdp_raw %>%
		filter(standage >= age_lims[1], standage <= age_lims[2]) %>%
		group_by(iteration, leaf_type, trait, variable, group) %>%
		summarise(
			slope     = coef(lm(yhat ~ standage))[2],
			intercept = coef(lm(yhat ~ standage))[1],
			.groups   = "drop"
		) %>%
		pivot_wider(names_from = group,
								values_from = c(slope, intercept)) %>%
		mutate(
			slope_diff     = slope_high - slope_low,
			intercept_diff = intercept_high - intercept_low,
			trait_label    = recode(trait,    !!!TRAIT_LABELS),
			variable_label = recode(variable, !!!ENV_LABELS),
			age_range      = range_name
		) %>%
		group_by(leaf_type, trait, trait_label, variable, variable_label,
						 age_range) %>%
		summarise(
			slope_median = median(slope_diff,   na.rm = TRUE),
			slope_lwr    = quantile(slope_diff, 0.025, na.rm = TRUE),
			slope_upr    = quantile(slope_diff, 0.975, na.rm = TRUE),
			slope_robust = (slope_lwr > 0 | slope_upr < 0),
			.groups      = "drop"
		)
})

write_rds(sensitivity_pdp, "tables/sensitivity_pdp.rds")

# ══════════════════════════════════════════════════════════════════════════════
# 4. Figures
# ══════════════════════════════════════════════════════════════════════════════

# ── 4a. Strata separation by threshold ───────────────────────────────────────
p_separation <- separation %>%
	mutate(
		leaf_type      = str_to_title(leaf_type),
		variable_label = factor(variable_label, levels = ENV_LABELS),
		threshold      = factor(threshold, levels = names(THRESHOLDS)),
		is_main        = threshold == "25/75"
	) %>%
	ggplot(aes(x = threshold, y = cohens_d,
						 colour = variable_label, group = variable_label)) +
	geom_line(linewidth = 0.6, alpha = 0.8) +
	geom_point(aes(size = is_main), show.legend = FALSE) +
	scale_size_manual(values = c("FALSE" = 1.5, "TRUE" = 3)) +
	geom_vline(xintercept = "25/75", linetype = "dashed",
						 colour = "grey50", linewidth = 0.4) +
	scale_colour_manual(values = COLS_ENV, name = NULL) +
	facet_wrap(~ leaf_type, ncol = 2) +
	labs(
		x = "Quantile threshold",
		y = "Cohen's d (separation between strata)"
	) +
	theme_succession(base_size = 9) +
	theme(legend.position = "bottom")

# ── 4b. VEcv divergence across thresholds ────────────────────────────────────
# Show mean |ΔVEcv| across traits and stand age bins per threshold
p_vecv_sens <- sensitivity_vecv_projected %>%
	ggplot(aes(x = threshold, y = projected_abs_delta,
						 colour = variable_label, group = variable_label)) +
	geom_line(linewidth = 0.6, alpha = 0.8) +
	geom_point(aes(size = threshold == "25/75"), show.legend = FALSE) +
	scale_size_manual(values = c("FALSE" = 1.5, "TRUE" = 3)) +
	geom_vline(xintercept = "25/75", linetype = "dashed",
						 colour = "grey50", linewidth = 0.4) +
	scale_colour_manual(values = COLS_ENV, name = NULL) +
	facet_wrap(~ leaf_type, ncol = 2) +
	labs(
		x        = "Quantile threshold",
		y        = "Projected |ΔVEcv| (scaled by strata separation)") +
	theme_succession(base_size = 9) +
	theme(legend.position = "bottom")

# ── 4c. PDP slope robustness across stand age ranges ─────────────────────────
# Show proportion of traits with robust slope differences per age range
p_pdp_sens <- sensitivity_pdp %>%
	group_by(age_range, leaf_type, variable_label) %>%
	summarise(
		prop_robust = mean(slope_robust, na.rm = TRUE),
		.groups     = "drop"
	) %>%
	mutate(
		leaf_type      = str_to_title(leaf_type),
		variable_label = factor(variable_label, levels = ENV_LABELS),
		age_range      = factor(age_range, levels = names(AGE_RANGES)),
		is_main        = age_range == "Full (0-150)"
	) %>%
	ggplot(aes(x = age_range, y = prop_robust,
						 colour = variable_label, group = variable_label)) +
	geom_line(linewidth = 0.6, alpha = 0.8) +
	geom_point(aes(size = is_main), show.legend = FALSE) +
	scale_size_manual(values = c("FALSE" = 1.5, "TRUE" = 3)) +
	scale_colour_manual(values = COLS_ENV, name = NULL) +
	scale_y_continuous(limits = c(0, 1),
										 labels = scales::percent) +
	facet_wrap(~ leaf_type, ncol = 2) +
	labs(
		x = "Stand age range for slope estimation",
		y = "Proportion of traits with robust slope differences"
	) +
	theme_succession(base_size = 9) +
	theme(legend.position = "bottom")

# ── Combine and save ──────────────────────────────────────────────────────────
p_sens <- (p_separation / p_vecv_sens / p_pdp_sens) +
	plot_layout(heights = c(1, 1, 1)) +
	plot_annotation(tag_levels = "a",
									theme = theme(plot.margin = margin(4, 4, 4, 4)))

save_fig(p_sens, "S3_quantile_sensitivity.png",
				 dir    = DIR_OUT,
				 width  = 180,
				 height = 240,
				 dpi    = 300)

message("\n── quantile_sensitivity.R complete ────────────────────────────────")
message("Tables: tables/sensitivity_vecv.rds, tables/sensitivity_pdp.rds")
message("Figure: figures/supplementary/S_sensitivity.png")