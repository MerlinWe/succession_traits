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
##
## Author: M. Weiss @ Maynard Lab UCL / ETH Zürich
################################################################################

rm(list = ls())
set.seed(42)

# ── Libraries ─────────────────────────────────────────────────────────────────
library(fastshap)
library(ranger)
library(tidyverse)

# ── Configuration ─────────────────────────────────────────────────────────────
N_SIM    <- 100L   # Monte Carlo draws per SHAP explanation (more = more stable)
PARALLEL <- TRUE   # parallelism handled inside fastshap::explain

PATH_DATA   <- "data_processed/fia_traits_clean.rds"
PATH_MODELS <- "models"
PATH_TABLES <- "tables"
PATH_PLOTS  <- "figures/supplementary"

source("scripts/functions.R")

# ── Vocabulary (must match 03_rf_fit.R) ───────────────────────────────────────

TRAITS <- c("bark_thickness", "conduit_diam", "height", "leaf_density",
						"leaf_k", "root_depth", "seed_dry_mass", "shade_tolerance",
						"specific_leaf_area")

COVARIATES <- c("standage", "temp_pc", "soil_pc", "rain_pc", "elevation", "soil_ph")
LEAF_TYPES <- c("broadleaf", "coniferous")
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
# 5. Console diagnostics
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