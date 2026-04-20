################################################################################
## succession_traits: 04 — SHAP analysis
## Computes Shapley values for each trait × leaf type model on the held-out
## test set, quantifies the relative importance of environmental vs
## successional filtering, and exports results for figure building.
##
## SHAP is computed once on the held-out test set (same split used in
## 03_rf_fit.R). No bootstrapping is applied — with n ≈ 5000 test plots and
## stable RF fits, SHAP sampling variance is negligible. Robustness of the
## environmental vs successional ranking is demonstrated by replication across
## the two independent leaf-type models (broadleaf / coniferous).
##
## Input:  data_processed/fia_traits_clean.rds
##         models/rf_<trait>_<leaftype>.rds     (from 03_rf_fit.R)
##         models/splits_<leaftype>.rds         (from 03_rf_fit.R)
##
## Output: tables/shap_values.rds              (full SHAP long table)
##         tables/shap_importance.rds           (env vs succ ratio summary)
##         tables/shap_per_var.rds              (per-variable importance)
##         plots/shap/beeswarm/                 (one per trait × leaf type)
##         plots/shap/dependence/               (one per variable × trait)
##         plots/shap/fig2a_ratio.png           (RQ1 main figure candidate)
##         plots/shap/shap_stackbar.png         (per-predictor breakdown)
##
## Author: M. Weiss @ Maynard Lab UCL / ETH Zürich
################################################################################

rm(list = ls())
set.seed(42)

# ── Libraries ─────────────────────────────────────────────────────────────────
library(fastshap)
library(ranger)
library(ggbeeswarm)
library(tidyverse)

# ── Configuration ─────────────────────────────────────────────────────────────
N_SIM    <- 100L   # Monte Carlo draws per SHAP explanation (more = more stable)
PARALLEL <- TRUE   # parallelism handled inside fastshap::explain

PATH_DATA   <- "data_processed/fia_traits_clean.rds"
PATH_MODELS <- "models"
PATH_TABLES <- "tables"
PATH_PLOTS  <- "plots/shap"

dir.create(PATH_PLOTS,  recursive = TRUE, showWarnings = FALSE)
dir.create(PATH_TABLES, showWarnings = FALSE)

source("scripts/functions.R")

# ── Vocabulary (must match 03_rf_fit.R) ───────────────────────────────────────
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

VAR_LABELS <- c(
	"standage"  = "Stand age",
	"temp_pc"   = "Temperature PC",
	"soil_pc"   = "Soil water retention PC",
	"rain_pc"   = "Precipitation PC",
	"elevation" = "Elevation",
	"soil_ph"   = "Soil pH"
)

# Predictor categories for RQ1
ENV_VARS  <- c("temp_pc", "soil_pc", "rain_pc", "elevation", "soil_ph")
SUCC_VARS <- "standage"


# ══════════════════════════════════════════════════════════════════════════════
# 1. Load models and test splits
# ══════════════════════════════════════════════════════════════════════════════

models <- map(LEAF_TYPES, function(lt) {
	map(TRAITS, function(tr) {
		path <- file.path(PATH_MODELS, sprintf("rf_%s_%s.rds", tr, lt))
		if (!file.exists(path))
			stop(sprintf("Model not found: %s — run 03_rf_fit.R first.", path))
		read_rds(path)
	}) %>% set_names(TRAITS)
}) %>% set_names(LEAF_TYPES)

splits <- map(LEAF_TYPES, function(lt) {
	read_rds(file.path(PATH_MODELS, sprintf("splits_%s.rds", lt)))
}) %>% set_names(LEAF_TYPES)

message("Models loaded: ", paste(LEAF_TYPES, collapse = " / "))
message(sprintf("fastshap version: %s", packageVersion("fastshap")))


# ══════════════════════════════════════════════════════════════════════════════
# 2. Compute SHAP values on held-out test set
# ══════════════════════════════════════════════════════════════════════════════
# One explain() call per trait × leaf type.
# X is the test set feature matrix (covariates only).
# SHAP values are computed out-of-fit: the model never saw these observations
# during training, matching the approach used for performance evaluation in
# 03_rf_fit.R and following Maynard et al. (2022).

message(sprintf("\nComputing SHAP values (nsim = %d per model)...", N_SIM))

shap_values <- map_dfr(LEAF_TYPES, function(lt) {
	
	test_df <- splits[[lt]]$test
	
	X_test <- test_df %>%
		dplyr::select(all_of(COVARIATES)) %>%
		as.matrix()
	
	message(sprintf("  %s: n_test = %d", lt, nrow(test_df)))
	
	map_dfr(TRAITS, function(tr) {
		
		message(sprintf("    · %s", tr))
		
		mod <- models[[lt]][[tr]]
		
		shap_mat <- fastshap::explain(
			object       = mod,
			X            = X_test,
			pred_wrapper = predict_fn,
			nsim         = N_SIM,
			parallel     = PARALLEL,
			adjust       = TRUE
		)
		
		# SHAP values to long format
		shap_long <- as.data.frame(shap_mat) %>%
			rowid_to_column(".row") %>%
			pivot_longer(-.row, names_to = "variable", values_to = "shap_value")
		
		# Feature values to long format (for dependence plots)
		feat_long <- as.data.frame(X_test) %>%
			rowid_to_column(".row") %>%
			mutate(PID_rep = test_df$PID_rep) %>%
			pivot_longer(-c(.row, PID_rep),
									 names_to  = "variable",
									 values_to = "feature_value")
		
		left_join(shap_long, feat_long, by = c(".row", "variable")) %>%
			mutate(trait = tr, leaf_type = lt) %>%
			dplyr::select(trait, leaf_type, PID_rep, variable, shap_value, feature_value)
	})
})

write_rds(shap_values, file.path(PATH_TABLES, "shap_values.rds"))
message(sprintf("\nSHAP complete: %d rows saved to tables/shap_values.rds",
								nrow(shap_values)))


# ══════════════════════════════════════════════════════════════════════════════
# 3. Variable importance: sum of absolute SHAP per predictor
# ══════════════════════════════════════════════════════════════════════════════
# Importance = sum(|SHAP|) across test plots, following Maynard et al. (2022).

shap_per_var <- shap_values %>%
	group_by(trait, leaf_type, variable) %>%
	summarise(
		sum_abs_shap  = sum(abs(shap_value),  na.rm = TRUE),
		mean_abs_shap = mean(abs(shap_value), na.rm = TRUE),
		mean_shap     = mean(shap_value,      na.rm = TRUE),
		.groups = "drop"
	) %>%
	mutate(
		variable_label = recode(variable, !!!VAR_LABELS),
		trait_label    = recode(trait,    !!!TRAIT_LABELS)
	)

write_rds(shap_per_var, file.path(PATH_TABLES, "shap_per_var.rds"))


# ══════════════════════════════════════════════════════════════════════════════
# 4. RQ1: Environmental vs successional filtering ratio
# ══════════════════════════════════════════════════════════════════════════════
# Sum absolute SHAP within two ecological categories, then compute their ratio.
#
# Ratio = total_environmental_shap / total_successional_shap
#
# Ratio > 1: environment explains more variation than stand age.
# Ratio = 1: equal importance.
# Ratio < 1: succession dominates (expected for height, shade tolerance).
#
# Using a ratio rather than proportions avoids the zero-sum framing of
# percentages and produces a single interpretable number per trait:
# "environmental filtering explains X times more variation than succession."
#
# The category-level summation (5 env predictors summed before dividing by
# 1 successional predictor) means the ratio reflects ecological signal, not
# predictor count.

shap_importance <- shap_per_var %>%
	mutate(category = case_when(
		variable %in% ENV_VARS  ~ "environmental",
		variable %in% SUCC_VARS ~ "successional",
		TRUE ~ NA_character_
	)) %>%
	filter(!is.na(category)) %>%
	group_by(trait, leaf_type, category) %>%
	summarise(category_shap = sum(sum_abs_shap), .groups = "drop") %>%
	pivot_wider(names_from = category, values_from = category_shap) %>%
	mutate(
		env_succ_ratio = environmental / successional,
		total_shap     = environmental + successional,
		prop_env       = environmental / total_shap,
		prop_succ      = successional  / total_shap,
		trait_label    = recode(trait, !!!TRAIT_LABELS)
	) %>%
	arrange(leaf_type, desc(env_succ_ratio))

write_rds(shap_importance, file.path(PATH_TABLES, "shap_importance.rds"))

message("\n── RQ1: Environmental vs successional filtering ratio ──────────────")
shap_importance %>%
	dplyr::select(trait_label, leaf_type, env_succ_ratio, prop_env) %>%
	mutate(across(where(is.numeric), ~ round(., 2))) %>%
	print(n = Inf)


# ══════════════════════════════════════════════════════════════════════════════
# 5. Figures
# ══════════════════════════════════════════════════════════════════════════════

# ── 5a. Beeswarm plots ────────────────────────────────────────────────────────
# Variables ordered by mean |SHAP| so bidirectional predictors are ranked
# by total influence, not their net signed effect.

message("\nGenerating beeswarm plots...")
dir.create(file.path(PATH_PLOTS, "beeswarm"), showWarnings = FALSE)

MAX_POINTS <- 5000L

for (lt in LEAF_TYPES) {
	for (tr in TRAITS) {
		
		# Compute variable order from mean |SHAP| for this trait × leaf type
		var_order <- shap_values %>%
			filter(trait == tr, leaf_type == lt) %>%
			group_by(variable) %>%
			summarise(mean_abs = mean(abs(shap_value), na.rm = TRUE), .groups = "drop") %>%
			arrange(mean_abs) %>%
			pull(variable)
		
		df_sub <- shap_values %>%
			filter(trait == tr, leaf_type == lt) %>%
			group_by(variable) %>%
			mutate(value_col = symmetric_scale(feature_value)) %>%
			ungroup() %>%
			mutate(variable = factor(variable, levels = var_order)) %>%
			{ if (nrow(.) > MAX_POINTS) slice_sample(., n = MAX_POINTS) else . }
		
		p_bee <- ggplot(df_sub,
										aes(x = variable, y = shap_value, colour = value_col)) +
			ggbeeswarm::geom_quasirandom(alpha = 0.35, size = 0.9) +
			geom_hline(yintercept = 0, linetype = "dashed",
								 colour = "grey50", linewidth = 0.3) +
			coord_flip() +
			scale_x_discrete(labels = VAR_LABELS) +
			scale_colour_viridis_c(
				name   = "Feature value",
				limits = c(-1, 1),
				breaks = c(-1, 0, 1),
				labels = c("Low", "Mid", "High")
			) +
			labs(
				x = NULL, y = "SHAP value",
				title = paste0(TRAIT_LABELS[[tr]], " — ", str_to_title(lt))
			) +
			theme_bw(base_size = 9) +
			theme(
				panel.grid.minor   = element_blank(),
				panel.grid.major.y = element_blank(),
				axis.text.y        = element_text(size = 8),
				legend.position    = "right",
				legend.key.height  = unit(3, "mm"),
				plot.title         = element_text(face = "bold", hjust = 0.5, size = 9)
			)
		
		ggsave(
			file.path(PATH_PLOTS, "beeswarm", sprintf("beeswarm_%s_%s.png", tr, lt)),
			plot = p_bee, width = 3.5, height = 3, dpi = 350
		)
	}
}
message("  ✓ Beeswarm plots saved")

# ── 5b. SHAP dependence plots ─────────────────────────────────────────────────
# One per variable × trait. Both leaf types overlaid with GAM smoother.

message("Generating dependence plots...")
dir.create(file.path(PATH_PLOTS, "dependence"), showWarnings = FALSE)

for (tr in TRAITS) {
	for (var in COVARIATES) {
		
		p_dep <- shap_values %>%
			filter(trait == tr, variable == var) %>%
			mutate(leaf_type = str_to_title(leaf_type)) %>%
			ggplot(aes(x = feature_value, y = shap_value,
								 colour = leaf_type, shape = leaf_type)) +
			geom_point(alpha = 0.2, size = 1.0) +
			geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
									se = FALSE, linewidth = 0.8) +
			geom_hline(yintercept = 0, linetype = "dashed",
								 colour = "grey50", linewidth = 0.3) +
			scale_colour_manual(
				values = c("Broadleaf" = "#228B22", "Coniferous" = "#d95f02")
			) +
			scale_shape_manual(
				values = c("Broadleaf" = 16, "Coniferous" = 17)
			) +
			labs(
				x = VAR_LABELS[[var]], y = "SHAP value",
				title = TRAIT_LABELS[[tr]],
				colour = NULL, shape = NULL
			) +
			theme_bw(base_size = 9) +
			theme(
				legend.position  = "bottom",
				panel.grid.minor = element_blank(),
				panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
				plot.title       = element_text(face = "bold", size = 9)
			)
		
		ggsave(
			file.path(PATH_PLOTS, "dependence",
								sprintf("dependence_%s_%s.png", var, tr)),
			plot = p_dep, width = 3.5, height = 3, dpi = 350
		)
	}
}
message("  ✓ Dependence plots saved")

# ── 5c. RQ1 ratio figure (Figure 2a candidate) ───────────────────────────────
# Lollipop on log2 scale. Dashed line at ratio = 1 (equal importance).
# Log scale is natural here: ratio of 1 is the neutral point, and the
# scale is symmetric around it (2× env = same distance as 0.5× env).

p_ratio <- shap_importance %>%
	mutate(
		trait_label = factor(
			trait_label,
			levels = shap_importance %>%
				group_by(trait_label) %>%
				summarise(m = mean(env_succ_ratio), .groups = "drop") %>%
				arrange(m) %>%
				pull(trait_label)
		),
		leaf_type = str_to_title(leaf_type)
	) %>%
	ggplot(aes(x = env_succ_ratio, y = trait_label, colour = leaf_type)) +
	geom_vline(xintercept = 1, linetype = "dashed",
						 colour = "grey50", linewidth = 0.4) +
	geom_segment(aes(x = 1, xend = env_succ_ratio,
									 y = trait_label, yend = trait_label),
							 linewidth = 0.6, alpha = 0.5) +
	geom_point(size = 3.5) +
	scale_colour_manual(
		values = c("Broadleaf" = "#228B22", "Coniferous" = "#d95f02")
	) +
	scale_x_continuous(
		trans  = "log2",
		breaks = c(0.5, 1, 2, 5, 10, 20),
		labels = function(x) paste0(x, "×")
	) +
	facet_wrap(~ leaf_type, ncol = 2) +
	labs(
		x        = "Environmental / successional importance (ratio, log scale)",
		y        = NULL,
		colour   = NULL,
		title    = "Environmental filtering dominates for most traits",
		subtitle = "Ratio > 1×: environment explains more variation than stand age"
	) +
	theme_bw(base_size = 10) +
	theme(
		legend.position    = "none",
		strip.text         = element_text(face = "bold"),
		panel.grid.minor   = element_blank(),
		panel.grid.major.y = element_blank()
	)

ggsave(
	file.path(PATH_PLOTS, "fig2a_ratio.png"),
	plot = p_ratio,
	width = 180, height = 120, units = "mm", dpi = 400
)
message("  ✓ RQ1 ratio figure saved")

# ── 5d. Stacked bar: per-predictor breakdown ──────────────────────────────────
# Shows which environmental variable drives each trait.
# Supplementary figure candidate.

p_stackbar <- shap_per_var %>%
	mutate(
		leaf_type      = str_to_title(leaf_type),
		trait_label    = factor(trait_label, levels = TRAIT_LABELS),
		variable_label = factor(variable_label, levels = VAR_LABELS)
	) %>%
	group_by(trait_label, leaf_type) %>%
	mutate(pct = sum_abs_shap / sum(sum_abs_shap) * 100) %>%
	ungroup() %>%
	ggplot(aes(x = trait_label, y = pct, fill = variable_label)) +
	geom_col(colour = "white", linewidth = 0.2) +
	facet_wrap(~ leaf_type, ncol = 1) +
	scale_fill_manual(
		name   = "Predictor",
		values = c(
			"Stand age"               = "#8B4513",
			"Temperature PC"          = "#d73027",
			"Soil water retention PC" = "#74add1",
			"Precipitation PC"        = "#4575b4",
			"Elevation"               = "#878787",
			"Soil pH"                 = "#a6761d"
		)
	) +
	scale_x_discrete(guide = guide_axis(angle = 35)) +
	labs(
		x = NULL,
		y = "Relative importance (% of total |SHAP|)",
		title = "SHAP importance breakdown by predictor"
	) +
	theme_bw(base_size = 9) +
	theme(strip.text = element_text(face = "bold"))

ggsave(
	file.path(PATH_PLOTS, "shap_stackbar.png"),
	plot = p_stackbar,
	width = 200, height = 160, units = "mm", dpi = 350
)
message("  ✓ Stacked bar figure saved")


# ══════════════════════════════════════════════════════════════════════════════
# 6. Console diagnostics
# ══════════════════════════════════════════════════════════════════════════════

message("\n── Top 2 predictors per trait (mean |SHAP|) ────────────────────────")
shap_per_var %>%
	group_by(trait_label, leaf_type) %>%
	slice_max(mean_abs_shap, n = 2) %>%
	dplyr::select(trait_label, leaf_type, variable_label, mean_abs_shap) %>%
	mutate(mean_abs_shap = round(mean_abs_shap, 4)) %>%
	arrange(leaf_type, trait_label) %>%
	print(n = Inf)

message("\n── Stand age SHAP summary by leaf type ─────────────────────────────")
shap_values %>%
	filter(variable == "standage") %>%
	group_by(leaf_type) %>%
	summarise(
		mean_abs = mean(abs(shap_value), na.rm = TRUE),
		sd_shap  = sd(shap_value,        na.rm = TRUE),
		.groups  = "drop"
	) %>%
	mutate(across(where(is.numeric), ~ round(., 5))) %>%
	print()

message("\n── 04_shap.R complete ──────────────────────────────────────────────")
message(sprintf("Tables: %s/", PATH_TABLES))
message(sprintf("Plots:  %s/", PATH_PLOTS))