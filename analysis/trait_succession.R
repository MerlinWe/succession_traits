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
library(broom)
library(rsample)
library(doParallel)
library(foreach)
library(forcats)
library(glue)
library(gtable)
library(cowplot)
library(gridExtra)
library(psych)
library(factoextra)
library(ggh4x)
library(ggpubr)
library(tidyverse)

export = TRUE      # Export plots?
parallel = TRUE    # Parallize? 32 cores if on threadripper - 8 if local
tuning = FALSE     # Model tuning or predefined parameters?
node_name <- Sys.info()["nodename"] # Check which device is running

# Set file paths conditionally
path_in <- "/Users/serpent/Documents/MSc/Thesis/Data/fia_traits"
path_out <- "/Users/serpent/Documents/MSc/Thesis/Code"

# Get functions
if (node_name != "Mac") {
	source(url("https://raw.githubusercontent.com/MerlinWe/succession_traits/main/analysis/functions.R"))
} else {
	source("/Users/serpent/Documents/MSc/Thesis/Code/analysis/functions.R")
}

# Set parallel cluster 
if (parallel) { 
	num_cores <-  ifelse(node_name == "threadeast", 4, 8)
	cl <- makeCluster(num_cores)
	registerDoParallel(cl, cores = num_cores)
	getDoParWorkers()
}

#bootstrap_results <- read_csv("/Users/serpent/Documents/MSc/Thesis/Code/output/data/bootsrap_shap.csv")
#pdp_data <- read_csv("/Users/serpent/Documents/MSc/Thesis/Code/output/data/pdp_data.csv")

# ---------- Read and prepare input data ----------

data <- read_csv(paste0(path_in, "/plotlevel_data_2024-07-15.csv")) %>%
	
	# Filter out standage above upper 10% quantiles and managed plots
	filter(standage < quantile(standage, 0.9)) %>%
	filter(managed == 0)

# Create binary dummy variables based on biome levels with >100 observations; 
# Drop data of remaining biomes.

data <- data %>%
	filter(biome %in% (
		data %>%
			count(biome) %>%
			filter(n > 100) %>%
			pull(biome)))

data <- model.matrix(~ biome - 1, data = data) %>%
	as_tibble() %>%
	rename_with(~ .x %>%
								str_replace_all("biome(?!_)", "biome_") %>%
								str_replace_all(" ", "_") %>%
								str_to_lower()) %>%
	bind_cols(data, .) %>%
	select(-biome)

## ---------- PCA to find environmental axes ----------

## Fit a PCA model to find axes representing underlying climatic and environmental factors. 
## Disregard elevation and soil ph to include in the model and control for their effect
## directly, and since there is no reason to assume underlying factors. 

climate_vars <- data %>% 
	dplyr::select(
		
		annual_mean_temperature,
		max_temperature_of_warmest_month,
		min_temperature_of_coldest_month,
		mean_temperature_of_coldest_quarter,
		mean_temperature_of_warmest_quarter,
		
		annual_precipitation,
		precipitation_of_driest_quarter,
		precipitation_of_wettest_quarter,
		precipitation_of_coldest_quarter,
		precipitation_of_warmest_quarter,
		
		sand_content_015cm,
		sand_content_060cm,
		water_capacity_015cm,
		water_capacity_060cm)

KMO(climate_vars) # check KMO
cortest.bartlett(cor(climate_vars), n = ncol(climate_vars)) # check sphericity

# Run PCA
climate_pca <- prcomp(climate_vars, scale. = TRUE, center = TRUE)
summary(climate_pca) 

# Varimax rotation to link variables with components
varimax_rotation <- varimax(climate_pca$rotation[, 1:3])
varimax_rotation <- varimax_rotation$loadings
print(varimax_rotation)

# Filter scores based on varimax loadings
scores <- round(get_pca_var(climate_pca)$coord[,1:3], digits = 3)

# Assign scores based on varimax rotated loadings
scores <- filter_scores(scores, varimax_rotation, threshold = 0.2) %>%
	as_tibble(rownames = "variable")

## >>> Build PCA plot for appendix <<<
scree_plot <- fviz_eig(climate_pca)
scree_plot <- scree_plot$data %>%
	as_tibble() %>%
	ggplot(aes(x = reorder(dim, -eig), y = eig, group = 1)) +  
	geom_bar(stat = "identity", color = "black", fill = "#1B7E74", alpha = 0.8) +
	geom_point() +
	geom_line(stat = "identity") +
	geom_vline(xintercept = 3.5, linetype = "dashed", colour = "red") +
	labs(y = "Percentage of Explained Variance", x = "Principal Components", title = "PC Eigenvalues") +
	theme_bw() +
	theme(text = element_text(family = "sans", size = 8))

# Contributions based on varimax rotation
varimax_contributions <- varimax_contrib(varimax_rotation)
abbreviations <- c(
	"annual_mean_temperature" = "AnnMeanTemp",
	"max_temperature_of_warmest_month" = "MaxTempWarmest",
	"min_temperature_of_coldest_month" = "MinTempColdest",
	"mean_temperature_of_coldest_quarter" = "MeanTempColdest",
	"mean_temperature_of_warmest_quarter" = "MeanTempWarmest",
	"annual_precipitation" = "AnnPrecip",
	"precipitation_of_driest_quarter" = "PrecipDriest",
	"precipitation_of_wettest_quarter" = "PrecipWettest",
	"precipitation_of_coldest_quarter" = "PrecipColdest",
	"precipitation_of_warmest_quarter" = "PrecipWarmest",
	"sand_content_015cm" = "Sand015",
	"sand_content_060cm" = "Sand060",
	"water_capacity_015cm" = "WaterCap015",
	"water_capacity_060cm" = "WaterCap060")

# Create contribution tibbles
var_contrib_pc1 <- tibble(name = rownames(varimax_contributions), contrib = varimax_contributions[, 1]) %>%
	create_contrib_plot("PC1 (Temperature) - Variable Contributions after Varimax Rotation", "firebrick4", abbreviations)
var_contrib_pc2 <- tibble(name = rownames(varimax_contributions), contrib = varimax_contributions[, 2]) %>%
	create_contrib_plot("PC2 (Soil Water Retention)- Variable Contributions after Varimax Rotation", "tan4", abbreviations)
var_contrib_pc3 <- tibble(name = rownames(varimax_contributions), contrib = varimax_contributions[, 3]) %>%
	create_contrib_plot("PC3 (Precipitation) - Variable Contributions after Varimax Rotation", "dodgerblue3", abbreviations)

# Combine the plots
pca_plot <- plot_grid(scree_plot, var_contrib_pc1, var_contrib_pc2, var_contrib_pc3,
											ncol = 2, nrow = 2, labels = "auto", label_fontfamily = "sans", label_size = 8)

if (export) {
	ggsave(paste0(path_out, "/output/plots/supplementary/s1.png"),
				 plot = pca_plot,
				 bg = "white",
				 width = 200,  
				 height = 130, 
				 units = "mm",
				 dpi = 600)
	}

# Write components
climate_pca <- climate_pca$x %>% 
	as_tibble() %>% 
	select(PC1, PC2, PC3) %>%
	rename(temp_pc = PC1,
				 soil_pc = PC2,
				 rain_pc = PC3) %>%
	mutate(PID_rep = data$PID_rep)

# Build database for further analysis 
data <- data %>% 
	left_join(climate_pca, by = "PID_rep") %>%
	select(starts_with("wmean_"), standage, temp_pc, soil_pc, rain_pc,
				 elevation, soil_ph_015cm, 
				 biome_boreal_forests_or_taiga, biome_flooded_grasslands, biome_mediterranean_woodlands, biome_temperate_broadleaf_forests,
				 biome_temperate_conifer_forests, biome_temperate_grasslands,biome_tundra, biome_xeric_shrublands, LAT, LON) %>%
	rename_with(~ gsub("_015cm", "", .), ends_with("_015cm")) %>%
	rename_with(~ gsub("wmean_", "", .), starts_with("wmean_")) %>%
	filter(complete.cases(.))

rm(climate_pca, climate_vars, pca_plot, scores, scree_plot, var_contrib_pc1,
	 var_contrib_pc2, var_contrib_pc3, abbreviations, varimax_rotation)

## ---------- Model tuning, training, and validation ----------

# Split into training and test sets
split <- initial_split(data, prop = 0.8)

# Extract training and testing datasets
train_data <- training(split)
test_data  <- testing(split)

# Define traits
traits <- c("wood_density", "bark_thickness", "conduit_diam", "leaf_n", 
						"specific_leaf_area", "seed_dry_mass", "shade_tolerance", "height")

# Define covariates 
covariates <- c("standage", "temp_pc", "soil_pc", "rain_pc", "elevation", "soil_ph", 
								"biome_boreal_forests_or_taiga", "biome_flooded_grasslands", 
								"biome_mediterranean_woodlands", "biome_temperate_broadleaf_forests",
								"biome_temperate_conifer_forests", "biome_temperate_grasslands",
								"biome_tundra", "biome_xeric_shrublands")

# Tuning if TRUE, else, fit model with predefined parameters
if (tuning) {
	
	hyper_grid <- expand.grid(
		num.trees = c(500, 1000, 1500),
		mtry = 2:4,
		min.node.size = c(1, 10, 20))
	
	# Get tuning results
	tuning_result <- map(traits, ~ tune_rf_model(.x, train_data, covariates, hyper_grid, num_threads = ifelse(node_name == "threadeast", 4, 1)))
	names(tuning_result) <- traits 
	hyper_grid <- map(tuning_result, "best_hyperparameters") %>% bind_rows()
	
	tuning_error_plot <- plot_grid(plotlist = map(tuning_result, "error_plot"), ncol = 2, nrow = 4)
	
	if (export) {
		ggsave(paste0(path_out, "/output/plots/supplementary/s2.png"),
					 plot = tuning_error_plot,
					 bg = "white",
					 width = 200,
					 height = 130,
					 units = "mm",
					 dpi = 600)
	}
	
	rm(tuning_error_plot, tuning_result)
	
} else { 
	# Use a default hyper_grid with known best values
	hyper_grid <- tibble(
		trait = traits,
		num.trees = rep(500, length(traits)),
		mtry = rep(4, length(traits)),
		min.node.size = rep(1, length(traits)))
}

# Fit models in parallel using foreach
best_models <- foreach(trait = traits, .packages = c('ranger', 'dplyr')) %dopar% {
	fit_rf_model(trait, data, covariates, hyper_grid, num_threads = ifelse(node_name == "threadeast", 4, 1))
}

# Extract performance metrics and models 
names(best_models) <- traits
performance_metrics <- lapply(best_models, `[[`, "performance") %>% bind_rows()
best_models <- map(best_models, "trait_mod")

if (export) {
	write_csv(performance_metrics, file = paste0(path_out, "/output/data/global_performance_metrics.csv"))
}

rm(split)

## ---------- Calculate Shapley values for each trait model ----------

# Define number of bootstrap iterations
num_bootstraps <- 100

bootstrap_results <- foreach(i = 1:num_bootstraps, .combine = 'rbind', .packages = c('dplyr', 'ranger', 'fastshap', 'tidyverse', 'foreach', 'rsample')) %dopar% {
	bootstrap_shap(data, traits, covariates, hyper_grid, num_threads = num_cores, boot_id = i)
}

shap_summary_boot <- bootstrap_results %>%
	mutate(category = case_when(
		variable == "standage" ~ "Successional Filtering",
		TRUE ~ "Environmental Filtering")) %>%
	group_by(trait, category, variable, bootstrap_id) %>%
	summarise(total_shap = sum(abs(shap_value)), .groups = "drop") %>%
	group_by(trait, category, bootstrap_id) %>%
	mutate(weight = total_shap / sum(total_shap)) %>%  # Normalize weights within each bootstrap
	summarise(weighted_shap = sum(total_shap * weight), .groups = "drop") %>%
	pivot_wider(names_from = category, values_from = weighted_shap, values_fill = 0) %>%
	
	# Compute confidence intervals for bootstrapped weighted SHAP values
	pivot_longer(-c(trait, bootstrap_id), names_to = "category", values_to = "shap_value") %>%
	group_by(trait, category) %>%
	summarise(
		median_shap = median(shap_value),  # Use median
		lower_ci = quantile(shap_value, probs = 0.025),
		upper_ci = quantile(shap_value, probs = 0.975),
		.groups = "drop")

# Export bootstrapped results
if (export) {

	write_csv(bootstrap_results, file = paste0(path_out, "/output/data/bootsrap_shap.csv"))
	write_csv(shap_summary_boot, paste0(path_out, "/output/data/shap_summary_boot.csv"))
	
	# Plot results
	shap_plot <- ggplot(shap_summary_boot, aes(x = median_shap, y = trait, color = category, shape = category)) +
		geom_point(size = 3) +
		geom_errorbar(aes(y = trait, xmin = lower_ci, xmax = upper_ci, color = category), width = 0.25) +
		scale_y_discrete(labels = c(
			"wood_density" = "Wood Density",
			"bark_thickness" = "Bark Thickness",
			"conduit_diam" = "Conduit Diameter",
			"leaf_n" = "Leaf Nitrogen",
			"specific_leaf_area" = "Specific Leaf Area",
			"seed_dry_mass" = "Seed Dry Mass",
			"shade_tolerance" = "Shade Tolerance",
			"height" = "Tree Height")) +
		scale_color_manual(values = c("Successional Filtering" = "black", "Environmental Filtering" = "darkred")) +
		scale_shape_manual(values = c("Successional Filtering" = 20, "Environmental Filtering" = 18)) + 
		theme_bw() +
		labs( x = "Total SHAP Importance (Weighted Normalization)",
					y = NULL,
					#caption = "Points represent bootstrapped median values. Error bars indicate 95% confidence intervals from bootstrapping") +
					color = "Predictor",
					shape = "Predictor") +
		theme(text = element_text(size = 11), 
					legend.position = c(.85, .85),
					legend.box.background = element_rect(color = "black"))
	
	ggsave(filename = paste0(path_out, "/output/plots/fig2.png"),
				 plot = shap_plot, 
				 bg = "transparent",
				 width = 220, 
				 height = 120, 
				 units = "mm", 
				 dpi = 800)
}

rm(bootstrap_results, shap_summary_boot) # clean environment

## ---------- Scenario Models for spatio-temporal interaction ----------

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
										"calculate_partial_dependence", "stratify", "traits", "covariates", 
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

# work from here on .. 

og <- bootstrap_pdp

# Recode group labels
bootstrapped_pdp <- bootstrap_pdp %>%
	mutate(group = recode_group(group))

# Compute bootstrapped slopes per trait × variable × iteration
slope_results <- bootstrapped_pdp %>%
	group_by(trait, variable, group, iteration) %>%
	summarise(slope = coef(lm(yhat ~ standage, data = cur_data()))[2], .groups = "drop") %>%
	pivot_wider(names_from = group, values_from = slope, names_prefix = "slope_") %>%
	mutate(slope_diff = abs(slope_high - slope_low))

# Compute bootstrapped intercepts
intercept_results <- bootstrapped_pdp %>%
	group_by(trait, variable, group, iteration) %>%
	summarise(intercept = coef(lm(yhat ~ standage, data = cur_data()))[1], .groups = "drop") %>%
	pivot_wider(names_from = group, values_from = intercept, names_prefix = "intercept_") %>%
	mutate(intercept_diff = abs(intercept_high - intercept_low))

# Merge slope and intercept results
filtering_results <- slope_results %>%
	left_join(intercept_results, by = c("trait", "variable", "iteration")) %>%
	select(trait, variable, iteration, slope_diff, intercept_diff)

# Compute t-test results per trait & environmental variable
t_test_pdp <- filtering_results %>%
	group_by(trait, variable) %>%
	summarise(
		t_p_value_slope = t.test(slope_diff)$p.value,
		t_p_value_intercept = t.test(intercept_diff)$p.value,
		.groups = "drop"
	)

# Compute Wilcoxon test results per trait & environmental variable
wilcox_test_pdp <- filtering_results %>%
	group_by(trait, variable) %>%
	summarise(
		wilcox_p_value_slope = wilcox.test(slope_diff)$p.value,
		wilcox_p_value_intercept = wilcox.test(intercept_diff)$p.value,
		.groups = "drop"
	)

# Merge statistical test results
final_results <- filtering_results %>%
	group_by(trait, variable) %>%
	summarise(
		slope_median = median(slope_diff),
		slope_lower = quantile(slope_diff, 0.025),
		slope_upper = quantile(slope_diff, 0.975),
		intercept_median = median(intercept_diff),
		intercept_lower = quantile(intercept_diff, 0.025),
		intercept_upper = quantile(intercept_diff, 0.975),
		.groups = "drop"
	) %>%
	left_join(t_test_pdp, by = c("trait", "variable")) %>%
	left_join(wilcox_test_pdp, by = c("trait", "variable")) %>%
	mutate(
		significant_slope = t_p_value_slope < 0.05,
		significant_intercept = t_p_value_intercept < 0.05
	)


ggplot(final_results, aes(x = intercept_median, y = slope_median, shape = trait, size = slope_upper - slope_lower)) +
	geom_point(aes(fill = variable, alpha = significant_slope | significant_intercept), color = "black") +
	#geom_errorbarh(aes(xmin = intercept_lower, xmax = intercept_upper), height = 0.05, color = "black") +
	#geom_errorbar(aes(ymin = slope_lower, ymax = slope_upper), width = 0.05, color = "black") +
	scale_shape_manual(values = c("wood_density" = 21, "specific_leaf_area" = 22, "shade_tolerance" = 23,
																"seed_dry_mass" = 24, "leaf_n" = 25, "tree_height" = 21,
																"conduit_diam" = 22, "bark_thickness" = 23)) +
	scale_size_continuous(range = c(.3, .7), name = "CI Range") +
	scale_fill_viridis_d(name = "Environmental Condition") +
	scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3), guide = "none") +
	facet_wrap(~variable, scales = "free", ncol = 2) +  
	labs(
		title = "Bootstrapped Successional and Environmental Filtering",
		x = "Intercept Difference (Initial Environmental Filtering)",
		y = "Slope Difference (Spatio-Temporal Interaction)"
	) +
	theme_bw() +
	theme(
		text = element_text(size = 8),
		legend.position = "right",
		strip.text = element_text(size = 6, face = "bold")
	)




slope_results <- pdp_data %>%
	mutate(group = recode_group(group)) %>%
	group_by(trait, variable, group) %>%
	summarise(
		slope = coef(lm(yhat ~ standage, data = cur_data()))[2], 
		.groups = "drop"
	) %>%
	pivot_wider(names_from = group, values_from = slope, names_prefix = "slope_") %>%
	mutate(slope_diff = abs(slope_high - slope_low))

# Compute intercept differences
intercept_results <- pdp_data %>%
	mutate(group = recode_group(group)) %>%
	group_by(trait, variable, group) %>%
	summarise(intercept = coef(lm(yhat ~ standage, data = cur_data()))[1], .groups = "drop") %>%
	pivot_wider(names_from = group, values_from = intercept, names_prefix = "intercept_") %>%
	mutate(intercept_diff = abs(intercept_high - intercept_low))

# Merge with slope differences
filtering_results <- slope_results %>%
	left_join(intercept_results, by = c("trait", "variable")) %>%
	select(trait, variable, slope_diff, intercept_diff)

ggplot(filtering_results, aes(x = intercept_diff, y = slope_diff, color = trait, group = trait)) +
	geom_point(size = 4, alpha = 0.8) +
	scale_color_viridis_d(name = "Feature") +
	labs(
		x = "Intercept Difference (Initial Environmental Filtering)",
		y = "Slope Difference (Spatio-Temporal Interaction)") +
	facet_wrap(~variable, scales = "fixed") +
	theme_bw() +
	theme(
		text = element_text(size = 14),
		legend.position = "top")









