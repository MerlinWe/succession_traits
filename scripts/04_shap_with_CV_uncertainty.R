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
##         tables/shap_importance_ci.rds        (env/succ ratio 95% CI via repeated CV; optional)
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
# 4b. Uncertainty on the env/succ ratio via repeated cross-validation (optional)
# ══════════════════════════════════════════════════════════════════════════════
# The ratio above is a single point estimate from one held-out test set. To show
# how robust it is to sampling and model-fit variability, we repeat k-fold CV a
# number of times, refit the RF on each training partition, recompute SHAP on the
# held-out fold, and form the env/succ ratio from that fold's summed |SHAP|.
# Folds are assigned at the PLOT level (PID), so all repeated inventories of a
# plot stay together and never leak between train and test. The 2.5/97.5
# percentiles across the repeat × fold ratios give a 95% interval per trait ×
# leaf type. The point estimates reported in the text are unchanged; this only
# adds error bars to Fig. 2a.
#
# Cost: refits / re-explains 18 models × N_REPEATS_CI × N_FOLDS_CI times. nsim is
# lower than the main run (aggregation across folds compensates) and the explained
# set can be capped (EXPLAIN_MAX). On a 32-core node the defaults run in a few
# hours; lower N_REPEATS_CI / EXPLAIN_MAX to trade precision for speed, or set
# RUN_SHAP_CI <- FALSE to skip entirely.

RUN_SHAP_CI  <- TRUE
N_FOLDS_CI   <- 10L
N_REPEATS_CI <- 5L
N_SIM_CI     <- 50L
EXPLAIN_MAX  <- 1500L   # max plots explained per fold (NULL = use all)

if (RUN_SHAP_CI) {
	
	library(doParallel)
	library(foreach)
	
	# Full modelling dataset with leaf type derived as in 03/05/06
	data_ci <- read_rds(PATH_DATA) %>%
		mutate(leaf_type = case_when(
			biome_boreal_forests_or_taiga    == 1 |
				biome_temperate_conifer_forests  == 1 ~ "coniferous",
			biome_temperate_broadleaf_forests == 1 |
				biome_mediterranean_woodlands    == 1 ~ "broadleaf",
			TRUE ~ NA_character_)) %>%
		filter(!is.na(leaf_type))
	
	# Tuned hyperparameters (same as the global models)
	hyper_ci <- map(LEAF_TYPES, function(lt)
		read_csv(file.path(PATH_TABLES, sprintf("perf_%s.csv", lt)),
						 show_col_types = FALSE) %>%
			dplyr::select(trait, num_trees, mtry, min_node_size)) %>%
		set_names(LEAF_TYPES)
	
	jobs <- expand.grid(trait = TRAITS, leaf_type = LEAF_TYPES,
											rep = seq_len(N_REPEATS_CI),
											stringsAsFactors = FALSE)
	
	n_cores <- min(32L, parallel::detectCores())
	cl <- makeCluster(n_cores)
	registerDoParallel(cl)
	message(sprintf("\n4b. SHAP CV uncertainty: %d jobs on %d cores...",
									nrow(jobs), n_cores))
	
	ratio_folds <- foreach(
		i         = seq_len(nrow(jobs)),
		.combine  = bind_rows,
		.packages = c("ranger", "fastshap", "dplyr", "purrr", "tibble", "tidyr"),
		.export   = c("fit_rf_model", "predict_fn", "COVARIATES", "ENV_VARS",
									"SUCC_VARS", "TRAITS", "data_ci", "hyper_ci",
									"N_FOLDS_CI", "N_SIM_CI", "EXPLAIN_MAX")
	) %dopar% {
		
		tr <- jobs$trait[i]; lt <- jobs$leaf_type[i]; rp <- jobs$rep[i]
		
		df_lt <- dplyr::filter(data_ci, leaf_type == lt)
		hg    <- hyper_ci[[lt]]
		
		# Plot-level fold assignment (reproducible per trait × repeat)
		set.seed(1000L * rp + match(tr, TRAITS))
		pids <- unique(df_lt$PID)
		fold_of_pid <- setNames(
			sample(rep(seq_len(N_FOLDS_CI), length.out = length(pids))),
			pids)
		df_lt$fold_id <- fold_of_pid[as.character(df_lt$PID)]
		
		purrr::map_dfr(seq_len(N_FOLDS_CI), function(k) {
			
			train_df <- df_lt %>%
				dplyr::filter(fold_id != k) %>%
				dplyr::select(dplyr::all_of(c(tr, COVARIATES))) %>%
				tidyr::drop_na()
			test_df  <- df_lt %>%
				dplyr::filter(fold_id == k) %>%
				tidyr::drop_na(dplyr::all_of(c(tr, COVARIATES)))
			
			if (nrow(train_df) < 200L || nrow(test_df) < 50L) return(NULL)
			
			mod <- fit_rf_model(tr, train_df, COVARIATES, hg,
													num_threads = 1L)$trait_mod
			
			X_test <- test_df %>% dplyr::select(dplyr::all_of(COVARIATES))
			if (!is.null(EXPLAIN_MAX) && nrow(X_test) > EXPLAIN_MAX)
				X_test <- X_test[sample(nrow(X_test), EXPLAIN_MAX), , drop = FALSE]
			
			sm <- fastshap::explain(mod, X = as.matrix(X_test),
															pred_wrapper = predict_fn, nsim = N_SIM_CI,
															parallel = FALSE, adjust = TRUE)
			abs_sum  <- colSums(abs(as.data.frame(sm)), na.rm = TRUE)
			env_sum  <- sum(abs_sum[ENV_VARS],  na.rm = TRUE)
			succ_sum <- sum(abs_sum[SUCC_VARS], na.rm = TRUE)
			
			tibble::tibble(
				trait = tr, leaf_type = lt, rep = rp, fold = k,
				env_succ_ratio = env_sum / succ_sum,
				prop_env       = env_sum / (env_sum + succ_sum))
		})
	}
	
	stopCluster(cl)
	
	shap_importance_ci <- ratio_folds %>%
		filter(is.finite(env_succ_ratio)) %>%
		group_by(trait, leaf_type) %>%
		summarise(
			n_cv         = n(),
			ratio_cv_med = median(env_succ_ratio),
			ratio_lwr    = quantile(env_succ_ratio, 0.025),
			ratio_upr    = quantile(env_succ_ratio, 0.975),
			prop_env_med = median(prop_env),
			prop_env_lwr = quantile(prop_env, 0.025),
			prop_env_upr = quantile(prop_env, 0.975),
			.groups = "drop"
		) %>%
		mutate(trait_label = recode(trait, !!!TRAIT_LABELS))
	
	write_rds(shap_importance_ci, file.path(PATH_TABLES, "shap_importance_ci.rds"))
	write_rds(ratio_folds,        file.path(PATH_TABLES, "shap_importance_cv_raw.rds"))
	
	message("\n── RQ1 env/succ ratio with CV 95% intervals ────────────────────────")
	shap_importance_ci %>%
		dplyr::select(trait_label, leaf_type,
									ratio_cv_med, ratio_lwr, ratio_upr, n_cv) %>%
		mutate(across(where(is.numeric), ~ round(., 2))) %>%
		arrange(leaf_type, desc(ratio_cv_med)) %>%
		print(n = Inf)
	
} else {
	message("\n4b. Skipping SHAP CV uncertainty (RUN_SHAP_CI = FALSE)")
}


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