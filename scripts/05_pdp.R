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
## The bootstrap (n = N_BOOT iterations) resamples which observations enter
## each stratum, providing uncertainty estimates for slope and intercept
## differences that are the primary test statistics for RQ2.
##
## Input:  data_processed/fia_traits_clean.rds
##         tables/perf_broadleaf.csv / perf_coniferous.csv  (hyperparameters)
##         tables/shap_per_var.rds                          (for SHAP weighting)
##
## Output: tables/pdp_raw.rds          (all bootstrap PDP curves)
##         tables/pdp_stats.rds        (slope/intercept diffs per iteration)
##         tables/pdp_summary.rds      (medians + CIs + tests)
##         plots/pdp/fig2b_ellipse.png (RQ2 main figure candidate)
##         plots/pdp/supp_pdp_lines.png
##
## Author: M. Weiss @ Maynard Lab UCL / ETH Zürich
################################################################################

rm(list = ls())
set.seed(42)

# ── Libraries ─────────────────────────────────────────────────────────────────
library(ranger)
library(doParallel)
library(foreach)
library(effsize)
library(ggh4x)
library(tidyverse)

# ── Configuration ─────────────────────────────────────────────────────────────
N_BOOT      <- 100L          # bootstrap iterations
PROBS_MAIN  <- c(0.25, 0.75) # main analysis quantile thresholds
PROBS_SENS  <- c(0.10, 0.90) # sensitivity analysis thresholds
RUN_SENS    <- FALSE         # set TRUE to also run 10/90 sensitivity
PARALLEL    <- TRUE
N_CORES     <- parallel::detectCores(logical = FALSE) - 1L

PATH_DATA   <- "data_processed/fia_traits_clean.rds"
PATH_TABLES <- "tables"
PATH_PLOTS  <- "plots/pdp"

dir.create(PATH_PLOTS,  recursive = TRUE, showWarnings = FALSE)
dir.create(PATH_TABLES, showWarnings = FALSE)

source("scripts/functions.R")

# ── Vocabulary ────────────────────────────────────────────────────────────────
TRAIT_LABELS <- c(
	"bark_thickness"     = "Bark Thickness",
	"conduit_diam"       = "Conduit Diameter",
	"height"             = "Tree Height",
	"leaf_density"       = "Leaf Density",
	"leaf_k"             = "Leaf Potassium",
	"root_depth"         = "Root Depth",
	"seed_dry_mass"      = "Seed Dry Mass",
	"shade_tolerance"    = "Shade Tolerance",
	"specific_leaf_area" = "Specific Leaf Area"
)
TRAITS     <- names(TRAIT_LABELS)
COVARIATES <- c("standage", "temp_pc", "soil_pc", "rain_pc", "elevation", "soil_ph")
LEAF_TYPES <- c("broadleaf", "coniferous")
ENV_VARS   <- c("temp_pc", "soil_pc", "rain_pc", "elevation", "soil_ph")

ENV_LABELS <- c(
	"temp_pc"   = "Temperature PC",
	"soil_pc"   = "Soil water retention PC",
	"rain_pc"   = "Precipitation PC",
	"elevation" = "Elevation",
	"soil_ph"   = "Soil pH"
)

VAR_LABELS <- c(
	"standage"  = "Stand age",
	"temp_pc"   = "Temperature PC",
	"soil_pc"   = "Soil water retention PC",
	"rain_pc"   = "Precipitation PC",
	"elevation" = "Elevation",
	"soil_ph"   = "Soil pH"
)


# ══════════════════════════════════════════════════════════════════════════════
# 1. Load data and hyperparameters
# ══════════════════════════════════════════════════════════════════════════════

data <- read_rds(PATH_DATA)

# Derive leaf type from biome dummies (consistent with 03_rf_fit.R)
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

# Load tuned hyperparameters (one table per leaf type, used for all strata)
hyper_params <- map(LEAF_TYPES, function(lt) {
	read_csv(file.path(PATH_TABLES, sprintf("perf_%s.csv", lt)),
					 show_col_types = FALSE) %>%
		dplyr::select(trait, num_trees, mtry, min_node_size)
}) %>% set_names(LEAF_TYPES)

# Load SHAP per-variable importance for SHAP-weighted mean panel
shap_per_var <- read_rds(file.path(PATH_TABLES, "shap_per_var.rds"))


# ══════════════════════════════════════════════════════════════════════════════
# 2. Core PDP bootstrap function
# ══════════════════════════════════════════════════════════════════════════════
# For one bootstrap iteration and one leaf type:
#   - resample data with replacement (standard bootstrap)
#   - for each environmental variable, stratify into high/low quantile groups
#   - fit RF models on each stratum
#   - compute partial dependence of stand age for each trait in each stratum
#   - return long table of PDP curves with iteration label

run_pdp_once <- function(boot_data, covariates, traits, hyper_grid, probs) {
	
	qtxt <- if (identical(probs, c(0.10, 0.90))) "10%" else "25%"
	
	# Build quantile labels programmatically to avoid repetition
	q_labels <- list(
		temp_pc   = c(paste0("Cold (lower ",  qtxt, ")"),
									paste0("Warm (upper ",  qtxt, ")")),
		soil_pc   = c(paste0("Sandy (lower ", qtxt, ")"),
									paste0("Water-retentive (upper ", qtxt, ")")),
		rain_pc   = c(paste0("Low precip (lower ",  qtxt, ")"),
									paste0("High precip (upper ", qtxt, ")")),
		elevation = c(paste0("Low elevation (lower ",  qtxt, ")"),
									paste0("High elevation (upper ", qtxt, ")")),
		soil_ph   = c(paste0("Low pH (lower ",  qtxt, ")"),
									paste0("High pH (upper ", qtxt, ")"))
	)
	
	map_dfr(names(q_labels), function(env_var) {
		
		labels  <- q_labels[[env_var]]
		strata  <- stratify_within(boot_data, env_var, probs)
		
		# Minimum cell size guard — skip if too few plots in either stratum
		if (nrow(strata$lower) < 50 || nrow(strata$upper) < 50) {
			warning(sprintf("Skipping %s: stratum too small (lower=%d, upper=%d)",
											env_var, nrow(strata$lower), nrow(strata$upper)))
			return(NULL)
		}
		
		# Fit models on lower and upper strata
		# num_threads = 1: parallelism is handled at the bootstrap level
		lower_models <- map(traits, ~ {
			m <- fit_rf_model(.x, strata$lower, covariates, hyper_grid, num_threads = 1L)
			m$trait_mod
		}) %>% set_names(traits)
		
		upper_models <- map(traits, ~ {
			m <- fit_rf_model(.x, strata$upper, covariates, hyper_grid, num_threads = 1L)
			m$trait_mod
		}) %>% set_names(traits)
		
		# Compute PDP for standage in each stratum for each trait
		map_dfr(traits, function(tr) {
			
			pdp_lower <- compute_pdp(lower_models[[tr]], strata$lower, "standage")
			pdp_upper <- compute_pdp(upper_models[[tr]], strata$upper, "standage")
			
			bind_rows(
				pdp_lower %>% mutate(group = "low",  group_label = labels[1]),
				pdp_upper %>% mutate(group = "high", group_label = labels[2])
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
		on.exit(stopCluster(cl), add = TRUE)
	} else {
		registerDoSEQ()
	}
	
	results <- foreach(
		iter      = seq_len(n_boot),
		.combine  = bind_rows,
		.packages = c("ranger", "dplyr", "purrr", "tidyr"),
		.export   = c("run_pdp_once", "stratify_within", "fit_rf_model",
									"compute_pdp", "predict_fn", "COVARIATES", "TRAITS")
	) %dopar% {
		
		set.seed(iter)
		
		# Standard bootstrap: resample with replacement at 100%
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
	
	results
}

# ── Run for both leaf types ───────────────────────────────────────────────────

message(sprintf("\nRunning PDP bootstrap (n = %d, probs = [%.2f, %.2f])...",
								N_BOOT, PROBS_MAIN[1], PROBS_MAIN[2]))

pdp_raw <- map_dfr(LEAF_TYPES, function(lt) {
	run_bootstrap(
		data      = data,
		leaf_type = lt,
		hyper_grid = hyper_params[[lt]],
		probs     = PROBS_MAIN,
		n_boot    = N_BOOT,
		n_cores   = if (PARALLEL) N_CORES else 1L
	)
})

write_rds(pdp_raw, file.path(PATH_TABLES, "pdp_raw.rds"))
message(sprintf("PDP bootstrap complete: %d rows saved to tables/pdp_raw.rds",
								nrow(pdp_raw)))

# ── Optional sensitivity analysis (10/90 quantiles) ──────────────────────────
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
# For each PDP curve (one per trait × variable × group × iteration), fit a
# linear model of yhat ~ standage to extract slope and intercept.
# Slope = rate of trait change through succession in that environmental context.
# Intercept = initial trait expression at the onset of succession.
# Signed difference (high - low) is retained to preserve directionality:
#   positive slope_diff = high-environment group has steeper successional trajectory
#   negative slope_diff = low-environment group has steeper trajectory

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
		# Signed differences: high minus low (direction is meaningful)
		slope_diff     = slope_high     - slope_low,
		intercept_diff = intercept_high - intercept_low,
		# Labels
		trait_label    = recode(trait,    !!!TRAIT_LABELS),
		variable_label = recode(variable, !!!ENV_LABELS)
	)

write_rds(pdp_stats, file.path(PATH_TABLES, "pdp_stats.rds"))


# ══════════════════════════════════════════════════════════════════════════════
# 5. SHAP-weighted mean summary
# ══════════════════════════════════════════════════════════════════════════════
# Compute a single summary per trait × leaf type × iteration that weights each
# environmental variable's slope/intercept difference by its relative SHAP
# importance. This gives an importance-weighted overall measure of
# spatio-temporal interaction that is comparable across traits.
#
# SHAP weights come from shap_per_var.rds (script 04). We normalise within
# trait × leaf type so weights sum to 1 across environmental variables.

shap_weights <- shap_per_var %>%
	filter(variable %in% ENV_VARS) %>%
	group_by(trait, leaf_type) %>%
	mutate(weight = sum_abs_shap / sum(sum_abs_shap)) %>%
	ungroup() %>%
	dplyr::select(trait, leaf_type, variable, weight)

# Join weights to PDP stats and compute weighted mean per iteration
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

# Combine per-variable and weighted mean for plotting
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
# 6. Summary statistics and significance tests
# ══════════════════════════════════════════════════════════════════════════════
# Paired Wilcoxon test on bootstrapped slope_high vs slope_low distributions,
# with BH correction across traits within each variable × leaf type.
# Cohen's dz quantifies effect size on the paired differences.

pdp_summary <- pdp_stats %>%
	group_by(leaf_type, trait, trait_label, variable, variable_label) %>%
	summarise(
		# Medians and 95% bootstrap CIs
		slope_median     = median(slope_diff,     na.rm = TRUE),
		slope_lwr        = quantile(slope_diff,   0.025, na.rm = TRUE),
		slope_upr        = quantile(slope_diff,   0.975, na.rm = TRUE),
		intercept_median = median(intercept_diff, na.rm = TRUE),
		intercept_lwr    = quantile(intercept_diff, 0.025, na.rm = TRUE),
		intercept_upr    = quantile(intercept_diff, 0.975, na.rm = TRUE),
		
		# Paired Wilcoxon test: are slope_high and slope_low distributions different?
		slope_p          = wilcox.test(slope_high, slope_low,
																	 paired = TRUE, exact = FALSE)$p.value,
		intercept_p      = wilcox.test(intercept_high, intercept_low,
																	 paired = TRUE, exact = FALSE)$p.value,
		
		# Cohen's dz effect size on paired differences
		slope_dz         = effsize::cohen.d(slope_high, slope_low,
																				paired = TRUE,
																				hedges.correction = TRUE)$estimate,
		intercept_dz     = effsize::cohen.d(intercept_high, intercept_low,
																				paired = TRUE,
																				hedges.correction = TRUE)$estimate,
		.groups = "drop"
	) %>%
	
	# BH correction within each variable × leaf type
	group_by(leaf_type, variable) %>%
	mutate(
		slope_p_adj     = p.adjust(slope_p,     method = "BH"),
		intercept_p_adj = p.adjust(intercept_p, method = "BH")
	) %>%
	ungroup() %>%
	mutate(
		slope_sig     = slope_p_adj     < 0.05,
		intercept_sig = intercept_p_adj < 0.05
	)

write_rds(pdp_summary, file.path(PATH_TABLES, "pdp_summary.rds"))

message("\n── RQ2: Significant slope differences (BH-adjusted p < 0.05) ──────")
pdp_summary %>%
	filter(slope_sig) %>%
	dplyr::select(leaf_type, trait_label, variable_label,
								slope_median, slope_dz, slope_p_adj) %>%
	mutate(across(where(is.numeric), ~ round(., 4))) %>%
	arrange(leaf_type, desc(abs(slope_dz))) %>%
	print(n = Inf)


# ══════════════════════════════════════════════════════════════════════════════
# 7. Figures
# ══════════════════════════════════════════════════════════════════════════════

# Summary medians for point overlay
pdp_summary_full <- pdp_stats_full %>%
	group_by(leaf_type, trait, trait_label, variable_label) %>%
	summarise(
		intercept_median = median(intercept_diff, na.rm = TRUE),
		slope_median     = median(slope_diff,     na.rm = TRUE),
		.groups = "drop"
	)

# ── 7a. Ellipse plot: Δintercept vs Δslope (Figure 2b candidate) ──────────────
# Each point = one bootstrap iteration (one trait, one variable, one leaf type).
# Ellipse = bootstrapped uncertainty cloud per trait.
# Median point overlaid.
# Split by leaf type (two panels) to show replication.
# Dashed lines at zero: quadrant tells the ecological story.
#   Top-right:  high env has both higher initial traits AND steeper trajectory
#   Bottom-right: high env has higher initial traits but flatter trajectory
#   Top-left:   high env has lower initial traits but steeper trajectory

# Shape palette: 9 traits need 9 distinguishable filled shapes
shape_vals <- c(21, 22, 23, 24, 25, 21, 22, 23, 24)
names(shape_vals) <- TRAIT_LABELS

p_ellipse <- pdp_stats_full %>%
	mutate(leaf_type = str_to_title(leaf_type)) %>%
	ggplot(aes(x = intercept_diff, y = slope_diff,
						 fill = trait_label, shape = trait_label)) +
	geom_vline(xintercept = 0, linetype = "dashed",
						 linewidth = 0.3, colour = "grey55") +
	geom_hline(yintercept = 0, linetype = "dashed",
						 linewidth = 0.3, colour = "grey55") +
	stat_ellipse(aes(group = interaction(leaf_type, trait_label)),
							 geom = "polygon", alpha = 0.12, colour = NA) +
	geom_point(
		data = pdp_summary_full %>% mutate(leaf_type = str_to_title(leaf_type)),
		aes(x = intercept_median, y = slope_median,
				fill = trait_label, shape = trait_label),
		size = 2.8, colour = "black", stroke = 0.4, inherit.aes = FALSE
	) +
	scale_fill_viridis_d(name = NULL) +
	scale_shape_manual(values = shape_vals, name = NULL) +
	guides(fill  = guide_legend(nrow = 3, override.aes = list(size = 2.5)),
				 shape = guide_legend(nrow = 3)) +
	facet_grid(leaf_type ~ variable_label, scales = "fixed") +
	labs(
		x = expression("Initial environmental filtering (" * Delta * " intercept, high \u2212 low)"),
		y = expression("Spatio-temporal interaction (" * Delta * " slope, high \u2212 low)")
	) +
	theme_bw(base_size = 9) +
	theme(
		strip.text       = element_text(face = "bold", size = 8),
		strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.5),
		legend.position  = "bottom",
		legend.key       = element_rect(fill = "white"),
		panel.grid.minor = element_blank()
	)

ggsave(
	file.path(PATH_PLOTS, "fig2b_ellipse.png"),
	plot = p_ellipse,
	width = 260, height = 180, units = "mm", dpi = 400
)
message("Figure 2b ellipse saved")

# ── 7b. Supplementary: PDP fit lines ─────────────────────────────────────────
# Shows the actual PDP curves (yhat vs standage) for high vs low groups.
# One panel per trait × variable. Faceted by leaf type.
# Points = bootstrap iterations; lines = linear fits.

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
	file.path(PATH_PLOTS, "supp_pdp_lines.png"),
	plot   = p_lines,
	width  = 260, height = 360, units = "mm", dpi = 300
)
message("Supplementary PDP lines saved")


# ══════════════════════════════════════════════════════════════════════════════
# 8. Console summary
# ══════════════════════════════════════════════════════════════════════════════

message("\n── Slope differences by variable (median across traits) ────────────")
pdp_summary %>%
	group_by(leaf_type, variable_label) %>%
	summarise(
		median_slope_diff = median(slope_median, na.rm = TRUE),
		n_sig_traits      = sum(slope_sig, na.rm = TRUE),
		.groups = "drop"
	) %>%
	arrange(leaf_type, desc(abs(median_slope_diff))) %>%
	print(n = Inf)

message("\n── 05_pdp.R complete ───────────────────────────────────────────────")
message(sprintf("Tables: %s/", PATH_TABLES))
message(sprintf("Plots:  %s/", PATH_PLOTS))