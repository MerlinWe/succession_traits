############################################################################################################################
########################################  succession_traits: main script (analysis) ########################################  
############################################################################################################################

rm(list = ls())   # make sure environment is clean 
set.seed(42)      # set seed for reproducibility

# ----- Session set-up -----

# Load necessary libraries
library(caret)
library(ranger) 
library(pdp)
library(mgcv)
library(car)
library(broom)
library(rsample)
library(doParallel)
library(foreach)
library(forcats)
library(effsize)
library(glue)
library(gtable)
library(cowplot)
library(gridExtra)
library(ggh4x)
library(ggpubr)
library(patchwork)
library(scales)
library(ggforce)
library(tidyverse)
library(beepr)

export   = TRUE     # Export plots?
tuning   = FALSE     # Model tuning or predefined parameters?
parallel = TRUE      # Run in parallel? 

node_name <- Sys.info()["nodename"] # Check which device is running

# Get functions
source("functions.R")

# Set parallel cluster 
if (parallel) { 
	num_cores <-  ifelse(node_name == "sea", 4, 8)
	cl <- makeCluster(num_cores)
	registerDoParallel(cl, cores = num_cores)
	getDoParWorkers()
}

# ---------- DATA ----------

data <- read_csv("data/fia_traits_clean.csv") %>%
	# Define leaf type
	mutate(
		leaf_type = case_when(
			biome_boreal_forests_or_taiga == 1 |
				biome_temperate_conifer_forests == 1 ~ "coniferous",
			biome_temperate_broadleaf_forests == 1 |
				biome_mediterranean_woodlands == 1 ~ "broadleaf",
			TRUE ~ NA_character_
		)
	) %>%
	# Keep only forests
	filter(!is.na(leaf_type))

# Define full trait name mapping
trait_labels <- c(
	"bark_thickness"    = "Bark Thickness",
	"conduit_diam"      = "Conduit Diameter",
	"height"            = "Tree Height",
	"leaf_density"      = "Leaf Density",
	"leaf_k"            = "Leaf Potassium",
	"root_depth"        = "Root Depth",
	"seed_dry_mass"     = "Seed Dry Mass",
	"shade_tolerance"   = "Shade Tolerance",
	"specific_leaf_area"= "Specific Leaf Area"
)

traits <- names(trait_labels)

all_covariates <- c(
	"standage", "temp_pc", "soil_pc", "rain_pc",
	"elevation", "soil_ph"
)

covariates_cont <- all_covariates
group_var        <- "leaf_type"
groups_to_analyze <- c("broadleaf", "coniferous")

## ---------- Model tuning, training, and validation ----------

#  Define grid and tuning flag 
#  Define grid and tuning flag 
if (tuning) {
	base_grid <- expand.grid(
		num.trees     = c(500, 1000, 1500),
		mtry          = 2:4,
		min.node.size = c(1, 10, 20)
	)
} else {
	# === Load tuned hyperparameters for each leaf type ===
	tuned_params <- list(
		broadleaf   = read_csv("tables/perf_broadleaf.csv", show_col_types = FALSE),
		coniferous  = read_csv("tables/perf_coniferous.csv", show_col_types = FALSE)
	)
}

# Initialize outputs 
fits_by_leaf <- list()
perf_by_leaf <- list()

## Loop over leaf types 
for (lt in groups_to_analyze) {
	message("\n=== Leaf type: ", lt, " ===")
	
	df_lt <- data %>%
		filter(leaf_type == lt) %>%
		mutate(.id = row_number())
	
	if (nrow(df_lt) < 1000)
		warning("Leaf type '", lt, "' has only ", nrow(df_lt), " rows; results may be unstable.")
	
	splits <- split_train_test(df_lt, ntest = min(5000, floor(0.2 * nrow(df_lt))), seed = 42)
	train_data <- splits$train
	test_data  <- splits$test
	
	# === TUNING =====================================
	
	if (tuning) {
		message("  → Running tuning phase (OOB error)…")
		
		tuned_list <- list()
		
		for (tr in traits) {
			message("    - Tuning ", tr)
			tune_res <- tune_rf_model(
				trait      = tr,
				data       = train_data,
				covariates = covariates_cont,
				hyper_grid = base_grid
			)
			
			tuned_list[[tr]] <- tune_res$best_hyperparameters
			
			# Save diagnostic plots
			ggsave(
				filename = file.path("plots", paste0("tuning_", lt, "_", tr, ".png")),
				plot = tune_res$error_plot, width = 4, height = 3, dpi = 300
			)
		}
		
		hyper_grid <- bind_rows(tuned_list)
		
	} else {
		
		# Skip tuning — use fixed base grid
		hyper_grid <- tuned_params[[lt]] %>%
			dplyr::select(trait, num_trees, mtry, min_node_size)
		
	}
	
	# ================================================================
	
	# ---- Fit random forest models per trait ----
	message("  → Fitting RF models…")
	
	best_models <- foreach(tr = traits, .packages = c("ranger", "dplyr")) %dopar% {
		fit_rf_model(
			trait            = tr,
			df_train         = train_data,
			covariates       = covariates_cont,
			hyper_parameters = hyper_grid,
			num_threads      = 4L
		)
	}
	names(best_models) <- traits
	
	# ---- Collect performance ----
	perf_tbl <- lapply(best_models, `[[`, "performance") %>%
		bind_rows() %>%
		mutate(leaf_type = lt)
	
	models_only <- map(best_models, "trait_mod")
	
	fits_by_leaf[[lt]] <- list(
		models     = models_only,
		test_data  = test_data,
		train_data = train_data,
		hyper_grid = hyper_grid
	)
	
	perf_by_leaf[[lt]] <- perf_tbl
	
	# Save performance summary 
	write_csv(perf_tbl, paste0("tables/perf_", lt, ".csv"))
	message("  ✓ Completed: ", lt)
}

## ---------- SHAP ----------

shap_by_leaf <- imap_dfr(fits_by_leaf, function(fit_obj, lt) {
	
	models <- fit_obj$models
	test_data <- fit_obj$test_data
	X_test <- test_data %>%
		select(.id, all_of(covariates_cont))
	
	imap_dfr(models, function(model, trait_name) {
		
		fs <- fastshap::explain(
			object       = model,
			X            = X_test %>% select(-.id) %>% as.matrix(),
			pred_wrapper = predict_fn,
			nsim         = 100,  # how many monte carlo draws?
			parallel     = TRUE,
			adjust       = TRUE
		)
		
		shap_long <- as.data.frame(fs) %>%
			rowid_to_column(".row") %>%
			pivot_longer(-.row, names_to = "variable", values_to = "shap_value")
		
		X_long <- X_test %>%
			rowid_to_column(".row") %>%
			pivot_longer(-c(.row, .id), names_to = "variable", values_to = "feature_value")
		
		left_join(shap_long, X_long, by = c(".row", "variable")) %>%
			mutate(trait = trait_name, leaf_type = lt) %>%
			select(trait, leaf_type, .id, variable, shap_value, feature_value)
	})
})

# export
beep(4)
if (export) write_csv(shap_by_leaf, "tables/shap_values_leaftype.csv")

shap_by_leaf <- read_csv("tables/shap_values_leaftype.csv")

#  NOTE! SHAP magnitudes are not directly comparable across leaf types (different baseline predictions). 
#  What is robust are the relative rankings of variables and the directional patterns of feature effects.
#	 Report SHAP with sample sizes per biome so reviewers don’t ding us for small-N effects.

# sanity check 
shap_by_leaf %>%
	group_by(leaf_type, variable) %>%
	summarise(mean_abs = mean(abs(shap_value), na.rm = TRUE), .groups="drop") %>%
	arrange(desc(mean_abs)) %>%
	group_split(leaf_type)

# Check beeswarms
for (tr in traits) {
	message(" → ", tr)
		p <- plot_beeswarm_by_leaf(shap_by_leaf, trait_label = tr)
		print(p)
		Sys.sleep(2) 
}

# Check standage scatter
for (tr in traits) {
	message(" → ", tr)
	p <- plot_shap_dependence(shap_by_leaf, var = "standage", trait_label = tr)
	print(p)
	Sys.sleep(2) 
}

# Check SD
for (tr in traits) {
	print(
		data %>%
			group_by(leaf_type) %>%
			summarise(sd = sd(!!sym(tr), na.rm = TRUE)) %>%
			mutate(trait = tr)
	)
}

# Check variance & mean SHAP for standage
shap_by_leaf %>%
	filter(variable == "standage") %>%
	group_by(leaf_type) %>%
	summarise(
		n          = n(),
		mean_shap  = mean(shap_value, na.rm = TRUE),
		sd_shap    = sd(shap_value, na.rm = TRUE),
		var_shap   = var(shap_value, na.rm = TRUE),
		mean_abs   = mean(abs(shap_value), na.rm = TRUE),
		.groups = "drop"
	)

# Broadleaves show ≈ 30–40 % greater SHAP variance and mean absolute SHAP for stand age.
# Stand age is a systematically stronger and more variable driver of trait prediction in broadleaves than in conifers.

# Check distributions
ggplot(shap_by_leaf %>% filter(variable == "standage"),
			 aes(x = shap_value, fill = leaf_type)) +
	geom_density(alpha = 0.5) +
	theme_bw(base_size = 11) +
	labs(title = "Distribution of stand-age SHAP values",
			 x = "SHAP value", y = "Density")

# Check variance
leveneTest(shap_value ~ leaf_type,
					 data = shap_by_leaf %>% 
					 	filter(variable == "standage"))

# Export plots to build compound figure (2)
if (export) {
	
	# ---- Pretty labels ----
	var_labels <- c(
		"standage"   = "Stand age",
		"temp_pc"    = "Temperature PC",
		"rain_pc"    = "Precipitation PC",
		"soil_pc"    = "Soil water retention PC",
		"elevation"  = "Elevation",
		"soil_ph"    = "Soil pH"
	)
	
	trait_labels <- c(
		"bark_thickness"     = "Bark thickness",
		"conduit_diam"       = "Conduit diameter",
		"height"             = "Height",
		"leaf_density"       = "Leaf density",
		"leaf_k"             = "Leaf potassium (K)",
		"root_depth"         = "Root depth",
		"seed_dry_mass"      = "Seed dry mass",
		"shade_tolerance"    = "Shade tolerance",
		"specific_leaf_area" = "Specific leaf area (SLA)"
	)
	
	# ---- SHAP Beeswarms ----
	max_points <- 5000
	
	for (tr in traits) {
		
		message("==== Processing trait: ", tr, " ====")
		tr_label <- trait_labels[[tr]]
		
		for (lt in c("Broadleaf", "Coniferous")) {
			
			df_sub <- shap_by_leaf %>%
				mutate(leaf_type = recode(leaf_type,
																	"broadleaf"   = "Broadleaf",
																	"coniferous"  = "Coniferous"
				)) %>%
				filter(trait == tr, leaf_type == lt) %>%
				group_by(variable) %>%
				mutate(value_col = symmetric_scale(feature_value)) %>%
				ungroup() %>%
				{ if (nrow(.) > max_points) slice_sample(., n = max_points) else . }
			
			p_bee <- ggplot(df_sub,
											aes(x = reorder(variable, shap_value, FUN = median),
													y = shap_value,
													color = value_col)) +
				ggbeeswarm::geom_quasirandom(alpha = 0.35, size = 0.9) +
				coord_flip() +
				scale_x_discrete(labels = var_labels) +
				scale_color_viridis_c(
					name = "Feature value",
					limits = c(-1, 1),
					breaks = c(-1, 0, 1),
					labels = c("Low", "Mid", "High")
				) +
				labs(
					x = NULL,
					y = "SHAP value",
					title = paste0(tr_label, " — ", lt)
				) +
				theme_bw(base_size = 9) +
				theme(
					panel.grid.minor = element_blank(),
					panel.grid.major.y = element_blank(),
					axis.text.y = element_text(size = 8, hjust = 1),
					axis.text.x = element_text(size = 8),
					strip.text = element_text(face = "bold", size = 8),
					legend.position = "none",
					plot.title = element_text(face = "bold", hjust = 0.5, size = 9)
				)
			
			ggsave(
				filename = file.path("plots/beeswarm",
														 paste0("beeswarm_", tr, "_", tolower(lt), ".png")),
				plot = p_bee,
				width = 3.5, height = 3, dpi = 350
			)
		}
		
		# ---- SHAP dependence (scatter) plots ----
		vars_tr <- shap_by_leaf %>%
			filter(trait == tr) %>%
			pull(variable) %>%
			unique()
		
		for (var in vars_tr) {
			
			var_label <- ifelse(var %in% names(var_labels),
													var_labels[[var]],
													paste0(var, " (feature value)"))
			
			p_scatter <- shap_by_leaf %>%
				mutate(leaf_type = recode(leaf_type,
																	"broadleaf"   = "Broadleaf",
																	"coniferous"  = "Coniferous"
				)) %>%
				filter(trait == tr, variable == var) %>%
				ggplot(aes(x = feature_value, y = shap_value,
									 color = leaf_type, shape = leaf_type)) +
				geom_point(alpha = 0.25, size = 1.4) +
				scale_color_manual(
					name = "Leaf type",
					values = c("Broadleaf" = "#228B22", "Coniferous" = "#d95f02")
				) +
				scale_shape_manual(
					name = "Leaf type",
					values = c("Broadleaf" = 16, "Coniferous" = 17)
				) +
				labs(
					x = var_label,
					y = "SHAP value"
				) +
				theme_bw(base_size = 9) +
				theme(
					legend.position = c(0.82, 0.16),
					legend.title = element_blank(),
					legend.background = element_rect(fill = "white", color = "grey70"),
					panel.grid.minor = element_blank(),
					panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
					axis.text = element_text(size = 8),
					axis.title = element_text(size = 8),
					plot.margin = margin(5, 5, 5, 5)
				)
			
			ggsave(
				filename = file.path("plots/scatter",
														 paste0("shap_scatter_", var, "_", tr, ".png")),
				plot = p_scatter,
				width = 3.5, height = 3, dpi = 350
			)
		}
	}
}


















# =========== PDPs by BIOME ==========

biome_col <- "biome_group"  
env_vars  <- c("temp_pc","soil_pc","rain_pc","elevation","soil_ph")

quantile_labels <- list(
	temp_pc = c("Cold Temperatures (Lower 25%)", "Warm Temperatures (Upper 25%)"),
	soil_pc = c("Sandy Soils (Lower 25%)", "Water-Retentive Soils (Upper 25%)"),
	rain_pc = c("Low Precipitation (Lower 25%)", "High Precipitation (Upper 25%)"),
	elevation = c("Low Elevation (Lower 25%)", "High Elevation (Upper 25%)"),
	soil_ph = c("Low Soil pH (Lower 25%)", "High Soil pH (Upper 25%)")
)

# switch to c(0.10, 0.90) for sensitivity (10/90):
probs_main <- c(0.25, 0.75)
probs_sens <- c(0.10, 0.90)

# Use cluster if set; else fallback to sequential
`%dopar_or_do%` <- if (foreach::getDoParWorkers() > 1) `%dopar%` else `%do%`

# RUN (25/75 main; 10/90 sensitivity optional)

num_bootstrap <- 100
num_threads   <- ifelse(exists("node_name") && node_name == "threadeast", 10, 1)

bootstrap_pdp_biome_main <- run_pdp_bootstrap_by_biome(
	data, traits, covariates, hyper_grid,
	num_bootstrap = num_bootstrap, num_threads = num_threads, probs = probs_main
)

# 90%/10% sensitivity:
bootstrap_pdp_biome_sens <- run_pdp_bootstrap_by_biome(
  data, traits, covariates, hyper_grid,
  num_bootstrap = num_bootstrap, num_threads = num_threads, probs = probs_sens
)













# -------------------------
# STATS: slopes & intercepts (signed high–low) per biome
# -------------------------
pdp_stats_biome <- bootstrap_pdp_biome_main %>%
	dplyr::mutate(group = recode_group(group)) %>%
	dplyr::filter(!is.na(group)) %>%
	dplyr::group_by(biome, trait, variable, group, iteration) %>%
	dplyr::summarise(
		slope = coef(lm(yhat ~ standage, data = dplyr::cur_data()))[2],
		intercept = coef(lm(yhat ~ standage, data = dplyr::cur_data()))[1],
		.groups = "drop"
	) %>%
	tidyr::pivot_wider(names_from = group, values_from = c(slope, intercept)) %>%
	dplyr::mutate(
		# keep the sign to show direction (high - low)
		slope_diff     = slope_high - slope_low,
		intercept_diff = intercept_high - intercept_low
	)

# If you have trait_labels (named vector) for prettier facets:
if (exists("trait_labels", inherits = FALSE)) {
	pdp_stats_biome <- pdp_stats_biome %>%
		dplyr::mutate(trait = factor(trait, levels = names(trait_labels), labels = trait_labels))
}

# -------------------------
# SIMPLE TESTS: paired t-tests and BH-adjusted p-values
# -------------------------
pdp_tests_biome <- pdp_stats_biome %>%
	dplyr::group_by(biome, trait, variable) %>%
	dplyr::summarise(
		slope_t_p     = t.test(slope_high, slope_low, paired = TRUE)$p.value,
		intercept_t_p = t.test(intercept_high, intercept_low, paired = TRUE)$p.value,
		slope_diff_mean     = mean(slope_diff, na.rm = TRUE),
		intercept_diff_mean = mean(intercept_diff, na.rm = TRUE),
		slope_diff_median   = median(slope_diff, na.rm = TRUE),
		intercept_diff_median = median(intercept_diff, na.rm = TRUE),
		slope_diff_lwr = quantile(slope_diff, 0.025, na.rm = TRUE),
		slope_diff_upr = quantile(slope_diff, 0.975, na.rm = TRUE),
		int_diff_lwr   = quantile(intercept_diff, 0.025, na.rm = TRUE),
		int_diff_upr   = quantile(intercept_diff, 0.975, na.rm = TRUE),
		.groups = "drop"
	) %>%
	dplyr::group_by(trait, variable) %>%
	dplyr::mutate(
		slope_p_adj     = p.adjust(slope_t_p, method = "BH"),
		intercept_p_adj = p.adjust(intercept_t_p, method = "BH")
	) %>%
	dplyr::ungroup()

# -------------------------
# FIGURE 3: ellipse plot of Δ-intercept vs Δ-slope, per biome
# -------------------------
pdp_stats_summary_biome <- pdp_stats_biome %>%
	dplyr::group_by(biome, trait, variable) %>%
	dplyr::summarise(
		intercept_median = median(intercept_diff, na.rm = TRUE),
		intercept_lwr    = quantile(intercept_diff, 0.025, na.rm = TRUE),
		intercept_upr    = quantile(intercept_diff, 0.975, na.rm = TRUE),
		slope_median     = median(slope_diff, na.rm = TRUE),
		slope_lwr        = quantile(slope_diff, 0.025, na.rm = TRUE),
		slope_upr        = quantile(slope_diff, 0.975, na.rm = TRUE),
		.groups = "drop"
	)

pdp_plot_biome <- ggplot(pdp_stats_biome,
												 aes(x = intercept_diff, y = slope_diff, fill = trait, shape = trait)) +
	geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey55") +
	geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey55") +
	ggplot2::stat_ellipse(aes(group = interaction(biome, trait)),
												geom = "polygon", alpha = 0.15, color = NA) +
	geom_point(data = pdp_stats_summary_biome,
						 aes(x = intercept_median, y = slope_median, fill = trait, shape = trait),
						 size = 2.8, color = "black", stroke = 0.4, inherit.aes = FALSE) +
	scale_fill_viridis_d(name = NULL) +
	scale_shape_manual(values = c(
		"Bark Thickness" = 21, "Conduit Diameter" = 22, "Tree Height" = 23,
		"Leaf Nitrogen" = 24, "Seed Dry Mass" = 25, "Shade Tolerance" = 21,
		"Specific Leaf Area" = 22, "Wood Density" = 23), name = NULL) +
	guides(fill = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
	facet_grid(biome ~ variable, scales = "fixed") +
	labs(
		x = expression("Initial Environmental Filtering (" * Delta * " intercept; high - low)"),
		y = expression("Spatio-temporal Interaction (" * Delta * " slope; high - low)")
	) +
	theme_bw() +
	theme(
		text = element_text(size = 12),
		strip.text = element_text(face = "bold"),
		strip.background = element_rect(fill = "white", colour = "black", linewidth = .75),
		legend.position = "top",
		legend.box = "horizontal",
		legend.key = element_rect(fill = "white", colour = "black")
	)

# -------------------------
# SUPPLEMENT: PDP fit lines (per biome × trait × variable)
# -------------------------
pdp_fit_biome <- bootstrap_pdp_biome_main %>%
	dplyr::mutate(group = recode_group(group)) %>%
	dplyr::filter(!is.na(group)) %>%
	ggplot(aes(x = standage, y = yhat)) +
	geom_point(aes(fill = group), alpha = 0.05, size = 0.5, shape = 21, stroke = 0.2, color = "black") +
	geom_smooth(aes(color = group), method = "lm", se = FALSE, linewidth = 0.7) +
	ggh4x::facet_nested(
		rows = vars(biome, trait),
		cols = vars(variable),
		scales = "free",
		labeller = if (exists("trait_labels", inherits = FALSE)) labeller(trait = trait_labels) else label_value
	) +
	scale_fill_manual(values = c("low" = "lightcyan", "high" = "lightcoral"),
										labels = c("high" = "Upper quantile", "low" = "Lower quantile"),
										name = NULL) +
	scale_color_manual(values = c("low" = "navy", "high" = "darkred"),
										 labels = c("high" = "Upper quantile", "low" = "Lower quantile"),
										 name = NULL) +
	labs(x = "Stand Age (years)", y = "Predicted Trait Expression") +
	theme_bw() +
	theme(
		text = element_text(size = 12),
		legend.position = "top",
		strip.text = element_text(size = 9, face = "bold"),
		legend.title = element_blank()
	)

# -------------------------
# EXPORTS (optional)
# -------------------------
# path_out must exist
# readr::write_csv(bootstrap_pdp_biome_main, file = file.path(path_out, "output/data/bootstrap_pdp_biome_main.csv"))
# readr::write_csv(pdp_tests_biome,          file = file.path(path_out, "output/data/pdp_tests_biome.csv"))

# ggsave(filename = file.path(path_out, "output/plots/fig3_biome.png"),
#        plot = pdp_plot_biome, bg = "transparent", width = 260, height = 180, units = "mm", dpi = 800)
# ggsave(filename = file.path(path_out, "output/plots/supplementary/pdp_fit_biome.png"),
#        plot = pdp_fit_biome, bg = "transparent", width = 260, height = 275, units = "mm", dpi = 800)




























## ---------- Scenario Models for spatio-temporal interaction per biome ----------

# Set Quantile values; covariates, num_threads; and quantile labels
quantiles_25 <- c(0.25, 0.75)

quantile_labels <- list(
	temp_pc = c("Cold Temperatures (Lower 25%)", "Warm Temperatures (Upper 25%)"),
	soil_pc = c("Sandy Soils (Lower 25%)", "Water-Retentive Soils (Upper 25%)"),
	rain_pc = c("Low Precipitation (Lower 25%)", "High Precipitation (Upper 25%)"),
	elevation = c("Low Elevation (Lower 25%)", "High Elevation (Upper 25%)"),
	soil_ph = c("Low Soil pH (Lower 25%)", "High Soil pH (Upper 25%)"))

# Set up parallel backend
num_bootstrap <- 100  # Number of bootstrap iterations
num_threads <- ifelse(node_name == "threadeast", 10, 1)

# Ensure all functions and variables are available in each worker
clusterExport(cl, c("fit_rf_model", "fit_models_on_strata", "calculate_pdp_for_scenario", 
										"calculate_partial_dependence", "stratify", "traits", "cont_covariates", 
										"hyper_grid", "quantile_labels"))

# Bootstrap fitting stratified scenario models and pdp calculation in parallel
bootstrap_pdp <- foreach(iter = 1:num_bootstrap, .combine = bind_rows, .packages = c("ranger", "foreach", "doParallel", "rsample", "tidyverse")) %dopar% {
	
	# Draw bootstrap sample (60-80% of the data)
	boot_data <- data %>% sample_frac(runif(1, 0.6, 0.8), replace = TRUE)
	
	# Fit models on stratified bootstrap data
	temp <- fit_models_on_strata(boot_data, "temp_pc", traits, covariates, hyper_grid, num_threads, quantile_labels$temp_pc)
	soil <- fit_models_on_strata(boot_data, "soil_pc", traits, covariates, hyper_grid, num_threads, quantile_labels$soil_pc)
	prcp <- fit_models_on_strata(boot_data, "rain_pc", traits, covariates, hyper_grid, num_threads, quantile_labels$rain_pc)
	elev <- fit_models_on_strata(boot_data, "elevation", traits, covariates, hyper_grid, num_threads, quantile_labels$elevation)
	ph   <- fit_models_on_strata(boot_data, "soil_ph", traits, covariates, hyper_grid, num_threads, quantile_labels$soil_ph)
	
	# Extract R² values for each model
	performance <- bind_rows(temp$performance, soil$performance, prcp$performance,
													 elev$performance, ph$performance) %>%
		dplyr::select(trait, group, rsq, pred_error) %>%
		mutate(iteration = iter)  # Add iteration ID
	
	# Compute Partial Dependence for each trait and scenario
	pdp_temp <- calculate_pdp_for_scenario(temp, traits, "standage", quantile_labels$temp_pc) %>% mutate(variable = "Temperature (PC)")
	pdp_soil <- calculate_pdp_for_scenario(soil, traits, "standage", quantile_labels$soil_pc) %>% mutate(variable = "Soil - Water Retention (PC)")
	pdp_prcp <- calculate_pdp_for_scenario(prcp, traits, "standage", quantile_labels$rain_pc) %>% mutate(variable = "Precipitation (PC)")
	pdp_elev <- calculate_pdp_for_scenario(elev, traits, "standage", quantile_labels$elevation) %>% mutate(variable = "Elevation")
	pdp_ph   <- calculate_pdp_for_scenario(ph, traits, "standage", quantile_labels$soil_ph) %>% mutate(variable = "Soil pH")
	
	# Combine all PDP results
	pdp_data <- bind_rows(pdp_temp, pdp_soil, pdp_prcp, pdp_elev, pdp_ph) %>%
		left_join(performance, by = c("trait", "group")) %>%
		mutate(iteration = iter)  # Add iteration ID
	
	return(pdp_data)
}

stopImplicitCluster()

if (export) {
	write_csv(bootstrap_pdp, file = paste0(path_out, "/output/data/bootstrap_pdp.csv"))
}

# ----- Compare slopes and intercepts -----

bootstrap_pdp <- read_csv("/Users/merlin/Documents/MSc/Thesis/Code/output/data/bootstrap_pdp.csv")

# Calculate slope & intercept, bootstrap confidence intervals, and t-tests
pdp_stats <- bootstrap_pdp %>%
	mutate(group = recode_group(group)) %>%
	
	# Compute bootstrapped slopes and intercepts
	group_by(trait, variable, group, iteration) %>%
	summarise(slope = coef(lm(yhat ~ standage, data = cur_data()))[2],
						intercept = coef(lm(yhat ~ standage, data = cur_data()))[1],
						.groups = "drop") %>%
	
	# Pivot to compare high vs. low environmental groups
	pivot_wider(names_from = group, values_from = c(slope, intercept), names_prefix = "") %>%
	
	# Compute absolute differences in slopes & intercepts
	mutate(slope_diff = abs(slope_high - slope_low),
				 intercept_diff = abs(intercept_high - intercept_low)) %>%
	mutate(trait = factor(trait, levels = names(trait_labels), labels = trait_labels))

## Plot results along a shap weighted mean 
pdp_stats_ellipses <- bootstrap_shap %>%
	
	# 1. Compute total SHAP per trait-variable-iteration
	dplyr::group_by(trait, variable, bootstrap_id) %>%
	dplyr::summarise(total_shap = sum(abs(shap_value)), .groups = "drop") %>%
	
	# 2. Filter to environmental variables and recode names
	filter(variable %in% c("elevation", "rain_pc", "soil_pc", "soil_ph", "temp_pc")) %>%
	mutate(variable = recode(variable,
													 elevation = "Elevation",
													 rain_pc = "Precipitation (PC)",
													 soil_pc = "Soil - Water Retention (PC)",
													 soil_ph = "Soil pH",
													 temp_pc = "Temperature (PC)"),
				 trait = recode(trait,
				 							 bark_thickness = "Bark Thickness",
				 							 conduit_diam = "Conduit Diameter",
				 							 height = "Tree Height",
				 							 leaf_n = "Leaf Nitrogen",
				 							 seed_dry_mass = "Seed Dry Mass",
				 							 shade_tolerance = "Shade Tolerance",
				 							 specific_leaf_area = "Specific Leaf Area",
				 							 wood_density = "Wood Density")) %>%
	
	# 3. Compute per-bootstrap SHAP weights
	group_by(trait, bootstrap_id) %>%
	mutate(weight = total_shap / sum(total_shap)) %>%
	rename(iteration = bootstrap_id) %>%
	select(trait, iteration, variable, weight) %>%
	
	# 4. Join to PDP slope/intercept diffs and compute weighted means
	right_join(pdp_stats, by = c("trait", "variable", "iteration")) %>%
	group_by(trait, iteration) %>%
	summarise(
		slope_diff = sum(slope_diff * weight),
		intercept_diff = sum(intercept_diff * weight),
		variable = "SHAP-Weighted Mean",
		.groups = "drop"
	) %>%
	
	# 5. Combine with unweighted PDP stats and set facet levels
	bind_rows(pdp_stats, .) %>%
	mutate(variable = factor(variable, levels = c(
		"SHAP-Weighted Mean", "Elevation",
		"Temperature (PC)", "Precipitation (PC)",
		"Soil - Water Retention (PC)", "Soil pH"
	)))

# Summarise bootstrap distributions per trait × variable
pdp_stats_summary_extended <- pdp_stats_ellipses %>%
  dplyr::group_by(trait, variable) %>%
  dplyr::summarise(
    intercept_median = median(intercept_diff, na.rm = TRUE),
    intercept_lwr    = quantile(intercept_diff, 0.025, na.rm = TRUE),
    intercept_upr    = quantile(intercept_diff, 0.975, na.rm = TRUE),
    slope_median     = median(slope_diff, na.rm = TRUE),
    slope_lwr        = quantile(slope_diff, 0.025, na.rm = TRUE),
    slope_upr        = quantile(slope_diff, 0.975, na.rm = TRUE),
    r_med            = sqrt(intercept_median^2 + slope_median^2),  # optional magnitude
    .groups = "drop"
  )

# Build the plot (Fig 3)
pdp_plot <- pdp_stats_ellipses %>% ggplot(aes(x = intercept_diff, y = slope_diff, fill = trait, shape = trait)) +
  # reference lines at zero (if you switch to signed deltas this becomes very informative)
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey55") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey55") +
  
	stat_ellipse(aes(group = trait), geom = "polygon", alpha = 0.2, color = NA) +
	geom_point(data = pdp_stats_summary_extended,
						 aes(x = intercept_median, y = slope_median, fill = trait, shape = trait),
						 size = 3, color = "black", stroke = 0.5, inherit.aes = FALSE) +
	scale_fill_viridis_d(name = NULL) +
	scale_shape_manual(values = c(
		"Bark Thickness" = 21, "Conduit Diameter" = 22, "Tree Height" = 23,
		"Leaf Nitrogen" = 24, "Seed Dry Mass" = 25, "Shade Tolerance" = 21,
		"Specific Leaf Area" = 22, "Wood Density" = 23), name = NULL) +
	guides(fill = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) + 
	facet_wrap(~variable, scales = "fixed", ncol = 2) +
	labs(
		x = expression("Initial Environmental Filtering (" * Delta * " intercept)"),
		y = expression("Spatio-temporal Interaction (" * Delta * " slope)")) +
	theme_bw() +
	theme(
		text = element_text(size = 12),
		strip.text = element_text(face = "bold"),
		strip.background = element_rect(fill = "white", colour = "black", linewidth = .75),
		legend.position = "top",
		legend.box = "horizontal",
		legend.key = element_rect(fill = "white", colour = "black"))

# Get summary statistics, compute Wilcoxon test and Cohen’s d for slopes and intercepts
pdp_stats_summary <- pdp_stats %>%
	group_by(trait, variable) %>%
	summarise(
		
		# Wilcoxon paired test
		slope_p_value = wilcox.test(slope_high, slope_low, paired = TRUE, exact = FALSE)$p.value,
		intercept_p_value = wilcox.test(intercept_high, intercept_low, paired = TRUE, exact = FALSE)$p.value,
		
		# Compute Cohen’s d for effect size
		slope_cohen_d = effsize::cohen.d(slope_high, slope_low, paired = TRUE, hedges.correction = TRUE)$estimate,
		intercept_cohen_d = effsize::cohen.d(intercept_high, intercept_low, paired = TRUE, hedges.correction = TRUE)$estimate,
		
		# Significance flags
		significant_slope = slope_p_value < 0.05,
		significant_intercept = intercept_p_value < 0.05,
		
		slope_median = median(slope_diff),
		slope_lower = quantile(slope_diff, 0.025),
		slope_upper = quantile(slope_diff, 0.975),
		intercept_median = median(intercept_diff),
		intercept_lower = quantile(intercept_diff, 0.025),
		intercept_upper = quantile(intercept_diff, 0.975),
		.groups = "drop") 

# Build supplementary figure 
pdp_fit <- bootstrap_pdp %>%
	mutate(group = recode_group(group)) %>%
	
	ggplot(aes(x = standage, y = yhat)) +
	
	# Points: use shape = 21 to allow fill + border
	geom_point(aes(fill = group), 
						 alpha = 0.05, size = 0.5, shape = 21, stroke = 0.2, color = "black") +
	
	# Lines: color by group
	geom_smooth(aes(color = group), method = "lm", se = FALSE, linewidth = 0.7) +
	
	# Facets
	ggh4x::facet_nested(
		rows = vars(trait),
		cols = vars(variable),
		scales = "free",
		labeller = labeller(trait = trait_labels)) +
	
	# Color and fill scales
	scale_fill_manual(values = c("low" = "lightcyan", "high" = "lightcoral"), 
										labels = c("high" = "Upper 25% quantile", "low" = "Lower 25% quantile"),
										name = NULL) +
	scale_color_manual(values = c("low" = "navy", "high" = "darkred"), 
										 labels = c("high" = "Upper 25% quantile", "low" = "Lower 25% quantile"),
										 name = NULL) +
	
	labs(x = "Stand Age (years)", y = "Predicted Trait Expression") +
	
	theme_bw() +
	theme(
		text = element_text(size = 12),
		legend.position = "top",
		strip.text = element_text(size = 9, face = "bold"),
		legend.title = element_blank())

# Export results
if (export) {
	write_csv(pdp_stats_summary, file = paste0(path_out, "/output/data/pdp_stats_summary.csv"))
	
	ggsave(filename = paste0(path_out, "/output/plots/fig3.png"),
				 plot = pdp_plot, 
				 bg = "transparent",
				 width = 260, 
				 height = 180, 
				 units = "mm", 
				 dpi = 800)
	
	ggsave(filename = paste0(path_out, "/output/plots/supplementary/s3.png"),
				 plot = pdp_fit, 
				 bg = "transparent",
				 width = 260, 
				 height = 275, 
				 units = "mm", 
				 dpi = 800)
	
}

### ------ Predictability vs Standage -------

# Apply to all traits
standage_mse <- purrr::map2_dfr(traits, best_models, function(trait, model) {
	compute_mse_by_bin(trait, model, data, env_vars = c("elevation", "rain_pc", "soil_pc", "soil_ph", "temp_pc"))
})

# Summarise MSE by standage bin and environmental group
divergence_by_env <- standage_mse %>%
	filter(env_group %in% c("high", "low")) %>%
	
	# Compute average MSE per variable, trait, standage_bin, and env_group
	group_by(variable, trait, standage_bin, env_group) %>%
	summarise(mean_mse = mean(mse, na.rm = TRUE), .groups = "drop") %>%
	
	# Pivot to wide format to get both high and low MSE in same row
	pivot_wider(names_from = env_group, values_from = mean_mse, names_prefix = "mse_") %>%
	
	# Compute absolute difference and mid standage
	mutate(mse_diff = abs(mse_high - mse_low),
				 standage_mid = as.numeric(str_extract(standage_bin, "(?<=\\[)\\d+")) + 5,
				 trait = recode(trait,
				 							 bark_thickness = "Bark Thickness",
				 							 conduit_diam = "Conduit Diameter",
				 							 height = "Tree Height",
				 							 leaf_n = "Leaf Nitrogen",
				 							 seed_dry_mass = "Seed Dry Mass",
				 							 shade_tolerance = "Shade Tolerance",
				 							 specific_leaf_area = "Specific Leaf Area",
				 							 wood_density = "Wood Density")) 

# Central line (mean of high/low env MSE per trait × standage)
plot_summary <- divergence_by_env %>%
	pivot_longer(cols = starts_with("mse_"),
							 names_to = "env_group",
							 names_prefix = "mse_",
							 values_to = "mse") %>%
	mutate(env_group = recode(env_group, high = "High Env", low = "Low Env")) %>%
	filter(env_group != "diff") %>%
	group_by(trait, standage_bin, env_group) %>%
	summarise(
		mse = mean(mse, na.rm = TRUE),
		standage_mid = as.numeric(str_extract(standage_bin, "(?<=\\[)\\d+")) + 5,
		.groups = "drop"
	)

# Highlight peak divergence (per trait: max absolute diff between high and low lines)
peak_points <- plot_summary %>%
	group_by(trait, standage_mid) %>%
	summarise(diff = abs(diff(mse)), .groups = "drop") %>%
	group_by(trait) %>%
	slice_max(diff, n = 1, with_ties = FALSE) %>%
	left_join(plot_summary, by = c("trait", "standage_mid"))

# Ribbon range
ribbon_summary <- standage_mse %>%
	filter(env_group %in% c("high", "low")) %>%
	mutate(
		standage_mid = as.numeric(str_extract(standage_bin, "(?<=\\[)\\d+")) + 5,
		env_group = recode(env_group, high = "High Env", low = "Low Env"),
		trait = recode(trait,
									 bark_thickness = "Bark Thickness",
									 conduit_diam = "Conduit Diameter",
									 height = "Tree Height",
									 leaf_n = "Leaf Nitrogen",
									 seed_dry_mass = "Seed Dry Mass",
									 shade_tolerance = "Shade Tolerance",
									 specific_leaf_area = "Specific Leaf Area",
									 wood_density = "Wood Density")
	) %>%
	group_by(trait, env_group, standage_mid) %>%
	summarise(
		mse_min = min(mse, na.rm = TRUE),
		mse_max = max(mse, na.rm = TRUE),
		.groups = "drop"
	)

# Final plot
plot_summary$env_group <- factor(plot_summary$env_group,
																 levels = c("Low Env", "High Env"),
																 labels = c("Lower Environmental Quantile", "Upper Environmental Quantile"))
ribbon_summary$env_group <- factor(ribbon_summary$env_group,
																 levels = c("Low Env", "High Env"),
																 labels = c("Lower Environmental Quantile", "Upper Environmental Quantile"))

predictability_plot <- ggplot() +
	geom_ribbon(data = ribbon_summary,
							aes(x = standage_mid, ymin = mse_min, ymax = mse_max, fill = env_group, group = env_group),
							alpha = 0.2) +
	geom_line(data = plot_summary,
						aes(x = standage_mid, y = mse, color = env_group, group = env_group),
						size = 1.2) +
	geom_point(data = peak_points,
						 aes(x = standage_mid, y = mse, color = env_group),
						 shape = 21, fill = "white", size = 2.5, stroke = 1) +
	facet_wrap(~trait, scales = "fixed", ncol = 4, nrow = 2) +
	scale_color_manual(
		values = c("Upper Environmental Quantile" = "#D95F02", 
							 "Lower Environmental Quantile" = "#1B9E77")
	) +
	scale_fill_manual(
		values = c("Upper Environmental Quantile" = "#D95F02", 
							 "Lower Environmental Quantile" = "#1B9E77")
	) +
	labs(
		x = "Stand Age (Years)",
		y = "Prediction Error (MSE)",
		color = "Environmental Group",
		fill = "Environmental Group"
	) +
	theme_bw() +
	theme(
		strip.background = element_rect(fill = "white", colour = "black", linewidth = .75),
		text = element_text(size = 12),
		legend.position = "top",
		strip.text = element_text(size = 9, face = "bold"),
		legend.title = element_blank())

ggsave(filename = paste0(path_out, "/output/plots/fig4.png"),
			 plot = predictability_plot, 
			 bg = "transparent",
			 width = 260, 
			 height = 160, 
			 units = "mm", 
			 dpi = 500)


# Rank divergence
divergence_by_env %>%
	# Now average across traits and bins to rank variables
	group_by(variable) %>%
	summarise(mean_mse_divergence = mean(mse_diff, na.rm = TRUE), .groups = "drop") %>%
	
	ggplot(aes(x = reorder(variable, mean_mse_divergence), 
						 y = mean_mse_divergence)) +
	geom_col(fill = "grey50", colour = "black", alpha = .7) +
	coord_flip() +
	labs(x = "Environmental Variable", 
			 y = "Mean Predictability Divergence (Δ MSE)") +
	theme_classic()

# Find most important environmental vars for divergence 
divergence_by_env %>%
	group_by(trait, variable) %>%
	summarise(avg_div = mean(mse_diff, na.rm = TRUE), .groups = "drop") %>%
	group_by(trait) 


