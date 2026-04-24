################################################################################
## results_numbers.R
## Pull and prints all quantitative values needed for the results
################################################################################

rm(list = ls())
library(tidyverse)
source("scripts/plot_theme.R")

perf        <- read_csv("tables/perf_combined.csv", show_col_types = FALSE)
shap_imp    <- read_rds("tables/shap_importance.rds")
shap_pervar <- read_rds("tables/shap_per_var.rds")
pdp_summary <- read_rds("tables/pdp_summary.rds")
vecv_sum    <- read_rds("tables/vecv_summary.rds")
vecv_div    <- read_rds("tables/vecv_divergence.rds")

ENV_VARS <- c("temp_pc", "soil_pc", "rain_pc", "elevation", "soil_ph")

cat("\n════════════════════════════════════════════════════════════════\n")
cat("MODEL PERFORMANCE\n")
cat("════════════════════════════════════════════════════════════════\n")

perf %>%
	summarise(
		mean_r2   = mean(r2_test),
		min_r2    = min(r2_test),
		max_r2    = max(r2_test),
		mean_bl   = mean(r2_test[leaf_type == "broadleaf"]),
		mean_co   = mean(r2_test[leaf_type == "coniferous"])
	) %>%
	mutate(across(everything(), ~ round(., 3))) %>%
	print()

cat("\nWorst performing models:\n")
perf %>% arrange(r2_test) %>%
	dplyr::select(trait, leaf_type, r2_test) %>%
	slice(1:4) %>% print()

cat("\n════════════════════════════════════════════════════════════════\n")
cat("RQ1: ENV/SUCC RATIOS\n")
cat("════════════════════════════════════════════════════════════════\n")

cat("\nFull ratio table:\n")
shap_imp %>%
	dplyr::select(trait_label, leaf_type, env_succ_ratio, prop_env) %>%
	mutate(across(where(is.numeric), ~ round(., 2))) %>%
	arrange(leaf_type, env_succ_ratio) %>%
	print(n = Inf)

cat("\nMean ratio by forest type:\n")
shap_imp %>%
	group_by(leaf_type) %>%
	summarise(
		mean_ratio   = mean(env_succ_ratio),
		median_ratio = median(env_succ_ratio),
		min_ratio    = min(env_succ_ratio),
		max_ratio    = max(env_succ_ratio),
		.groups = "drop"
	) %>%
	mutate(across(where(is.numeric), ~ round(., 2))) %>%
	print()

cat("\nTraits where env_succ_ratio < 5 (succession-relevant):\n")
shap_imp %>%
	filter(env_succ_ratio < 5) %>%
	dplyr::select(trait_label, leaf_type, env_succ_ratio) %>%
	mutate(across(where(is.numeric), ~ round(., 2))) %>%
	print()

cat("\nTop environmental predictor per trait × leaf type:\n")
shap_pervar %>%
	filter(variable %in% ENV_VARS) %>%
	group_by(trait_label, leaf_type) %>%
	slice_max(mean_abs_shap, n = 1) %>%
	dplyr::select(trait_label, leaf_type, variable_label, mean_abs_shap) %>%
	mutate(mean_abs_shap = round(mean_abs_shap, 4)) %>%
	arrange(leaf_type, desc(mean_abs_shap)) %>%
	print(n = Inf)

cat("\nTemp PC dominance in broadleaf (% of traits where temp is top predictor):\n")
shap_pervar %>%
	filter(variable %in% ENV_VARS, leaf_type == "broadleaf") %>%
	group_by(trait_label) %>%
	slice_max(mean_abs_shap, n = 1) %>%
	summarise(prop_temp = mean(variable == "temp_pc")) %>%
	summarise(prop_temp_dominant = mean(prop_temp)) %>%
	print()

cat("\nElevation dominance in coniferous:\n")
shap_pervar %>%
	filter(variable %in% ENV_VARS, leaf_type == "coniferous") %>%
	group_by(trait_label) %>%
	slice_max(mean_abs_shap, n = 1) %>%
	count(variable_label) %>%
	arrange(desc(n)) %>%
	print()

cat("\n════════════════════════════════════════════════════════════════\n")
cat("RQ2: PDP SLOPE DIFFERENCES\n")
cat("════════════════════════════════════════════════════════════════\n")

cat("\nNumber of robust slope differences per env variable × leaf type:\n")
pdp_summary %>%
	group_by(leaf_type, variable_label) %>%
	summarise(
		n_robust         = sum(slope_robust, na.rm = TRUE),
		median_abs_slope = round(median(abs(slope_median), na.rm = TRUE), 5),
		.groups = "drop"
	) %>%
	arrange(leaf_type, desc(n_robust)) %>%
	print(n = Inf)

cat("\nLargest robust slope differences (top 10 by |median slope|):\n")
pdp_summary %>%
	filter(slope_robust) %>%
	arrange(desc(abs(slope_median))) %>%
	dplyr::select(leaf_type, trait_label, variable_label,
								slope_median, slope_lwr, slope_upr) %>%
	mutate(across(where(is.numeric), ~ round(., 5))) %>%
	slice(1:10) %>%
	print()

cat("\nDirection summary — resource-acquisitive vs stress-tolerance traits:\n")
cat("(positive slope_diff = high env steeper; negative = low env steeper)\n")
pdp_summary %>%
	filter(variable_label == "Temperature PC", slope_robust) %>%
	dplyr::select(leaf_type, trait_label, slope_median) %>%
	mutate(
		direction = if_else(slope_median > 0, "high_env_steeper", "low_env_steeper"),
		slope_median = round(slope_median, 5)
	) %>%
	arrange(leaf_type, direction) %>%
	print(n = Inf)

cat("\nIntercept differences (initial environmental filtering):\n")
pdp_summary %>%
	filter(slope_robust) %>%
	group_by(leaf_type, variable_label) %>%
	summarise(
		median_abs_intercept = round(median(abs(intercept_median), na.rm = TRUE), 4),
		n_robust_intercept   = sum(intercept_robust, na.rm = TRUE),
		.groups = "drop"
	) %>%
	arrange(leaf_type, desc(median_abs_intercept)) %>%
	print(n = Inf)

cat("\n════════════════════════════════════════════════════════════════\n")
cat("RQ3: VECV PREDICTABILITY AND DIVERGENCE\n")
cat("════════════════════════════════════════════════════════════════\n")

cat("\nMean VEcv by trait and forest type:\n")
vecv_sum %>%
	group_by(trait_label, leaf_type) %>%
	summarise(mean_VEcv = round(mean(VEcv_med, na.rm = TRUE), 3),
						.groups = "drop") %>%
	pivot_wider(names_from = leaf_type, values_from = mean_VEcv) %>%
	print(n = Inf)

cat("\nOverall mean VEcv by forest type:\n")
vecv_sum %>%
	group_by(leaf_type) %>%
	summarise(mean_VEcv = round(mean(VEcv_med, na.rm = TRUE), 3),
						.groups = "drop") %>%
	print()

cat("\nVEcv at early vs late succession (mean across traits):\n")
vecv_sum %>%
	mutate(succession_stage = case_when(
		standage_mid <= 30  ~ "early (0-30)",
		standage_mid >= 100 ~ "late (100+)",
		TRUE ~ "mid"
	)) %>%
	filter(succession_stage != "mid") %>%
	group_by(leaf_type, succession_stage) %>%
	summarise(mean_VEcv = round(mean(VEcv_med, na.rm = TRUE), 3),
						.groups = "drop") %>%
	pivot_wider(names_from = succession_stage, values_from = mean_VEcv) %>%
	print()

cat("\nTraits with strongest predictability divergence (mean |ΔVEcv|):\n")
vecv_div %>%
	group_by(trait_label, leaf_type) %>%
	summarise(
		mean_abs_delta = round(mean(abs(delta_med), na.rm = TRUE), 3),
		prop_sig       = round(mean(sig_divergence, na.rm = TRUE), 3),
		.groups = "drop"
	) %>%
	arrange(desc(mean_abs_delta)) %>%
	print(n = Inf)

cat("\nEnvironmental variables ranked by divergence strength:\n")
vecv_div %>%
	group_by(leaf_type, variable_label) %>%
	summarise(
		mean_abs_delta = round(mean(abs(delta_med), na.rm = TRUE), 3),
		prop_sig       = round(mean(sig_divergence, na.rm = TRUE), 3),
		.groups = "drop"
	) %>%
	arrange(leaf_type, desc(mean_abs_delta)) %>%
	print(n = Inf)

cat("\nBroadleaf vs coniferous divergence comparison:\n")
vecv_div %>%
	group_by(leaf_type) %>%
	summarise(
		mean_abs_delta = round(mean(abs(delta_med), na.rm = TRUE), 3),
		prop_sig       = round(mean(sig_divergence, na.rm = TRUE), 3),
		.groups = "drop"
	) %>%
	print()

cat("\nDoes divergence increase or decrease through succession?\n")
cat("(slope of |delta_med| ~ standage_mid, averaged across traits)\n")
vecv_div %>%
	filter(!is.na(delta_med)) %>%
	group_by(leaf_type, variable_label) %>%
	summarise(
		divergence_trend = coef(lm(abs(delta_med) ~ standage_mid))[2],
		.groups = "drop"
	) %>%
	mutate(
		direction = if_else(divergence_trend > 0, "increasing", "decreasing"),
		divergence_trend = round(divergence_trend, 6)
	) %>%
	group_by(leaf_type, direction) %>%
	summarise(n = n(), .groups = "drop") %>%
	print()

cat("\n════════════════════════════════════════════════════════════════\n")
cat("results_numbers.R complete\n")
cat("════════════════════════════════════════════════════════════════\n")