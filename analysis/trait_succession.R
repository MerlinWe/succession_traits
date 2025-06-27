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
library(effsize)
library(glue)
library(gtable)
library(cowplot)
library(gridExtra)
library(psych)
library(factoextra)
library(ggh4x)
library(ggpubr)
library(tidyverse)

export   = FALSE     # Export plots?
parallel = TRUE      # Parallize? 32 cores if on threadripper - 8 if local
tuning   = FALSE     # Model tuning or predefined parameters?
node_name <- Sys.info()["nodename"] # Check which device is running

# Set file paths conditionally
path_in <- "/Users/merlin/Documents/MSc/Thesis/Data/fia_traits"
path_out <- "/Users/merlin/Documents/MSc/Thesis/Code"

# Get functions
if (node_name != "Mac") {
	source(url("https://raw.githubusercontent.com/MerlinWe/succession_traits/main/analysis/functions.R"))
} else {
	source("/Users/merlin/Documents/MSc/Thesis/Code/analysis/functions.R")
}

# Set parallel cluster 
if (parallel) { 
	num_cores <-  ifelse(node_name == "threadeast", 4, 8)
	cl <- makeCluster(num_cores)
	registerDoParallel(cl, cores = num_cores)
	getDoParWorkers()
}

# Define full trait name mapping
trait_labels <- c(
	"bark_thickness" = "Bark Thickness",
	"conduit_diam" = "Conduit Diameter",
	"height" = "Tree Height",
	"leaf_n" = "Leaf Nitrogen",
	"seed_dry_mass" = "Seed Dry Mass",
	"shade_tolerance" = "Shade Tolerance",
	"specific_leaf_area" = "Specific Leaf Area",
	"wood_density" = "Wood Density")

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

bootstrap_shap <- foreach(i = 1:num_bootstraps, .combine = 'rbind', .packages = c('dplyr', 'ranger', 'fastshap', 'tidyverse', 'foreach', 'rsample')) %dopar% {
	bootstrap_shap(data, traits, covariates, hyper_grid, num_threads = num_cores, boot_id = i)
}

# Summarize SHAP values per trait × bootstrap × variable
shap_summary_boot <- bootstrap_shap %>%
	group_by(trait, variable, bootstrap_id) %>%
	summarise(total_shap = sum(abs(shap_value)), .groups = "drop") %>%
	
	mutate(category = case_when(
		variable == "standage" ~ "Successional Filtering",
		TRUE ~ "Environmental Filtering"
	)) %>%
	
	mutate(trait = recode(trait,
												bark_thickness = "Bark Thickness",
												conduit_diam = "Conduit Diameter",
												height = "Tree Height",
												leaf_n = "Leaf Nitrogen",
												seed_dry_mass = "Seed Dry Mass",
												shade_tolerance = "Shade Tolerance",
												specific_leaf_area = "Specific Leaf Area",
												wood_density = "Wood Density")) %>%
	
	group_by(trait, category, bootstrap_id) %>%
	mutate(weight = total_shap / sum(total_shap)) %>%
	summarise(weighted_shap = sum(total_shap * weight), .groups = "drop")


# Compute paired differences and standardized effect size (dz)
shap_effects <- shap_summary_boot %>%
	
	pivot_wider(names_from = category, values_from = weighted_shap) %>%
	group_by(trait) %>%
	summarise(
		median_successional = median(`Successional Filtering`),
		median_environmental = median(`Environmental Filtering`),
		
		# Bootstrap delta SHAP differences
		shap_diff = list(`Successional Filtering` - `Environmental Filtering`),
		
		difference = median(`Successional Filtering` - `Environmental Filtering`),
		lower_ci = quantile(`Successional Filtering` - `Environmental Filtering`, 0.025),
		upper_ci = quantile(`Successional Filtering` - `Environmental Filtering`, 0.975),
		
		# Paired Cohen’s d_z
		d_z = mean(`Successional Filtering` - `Environmental Filtering`) / 
			sd(`Successional Filtering` - `Environmental Filtering`),
		
		# Bootstrap CI for d_z
		d_z_ci = list({
			diffs <- `Successional Filtering` - `Environmental Filtering`
			replicate(1000, {
				resampled <- sample(diffs, replace = TRUE)
				mean(resampled) / sd(resampled)
			}) %>% quantile(c(0.025, 0.975), na.rm = TRUE)
		}),
		
		.groups = "drop"
	) %>%
	
	# Unpack CI list into columns
	mutate(
		d_z_lower = purrr::map_dbl(d_z_ci, 1),
		d_z_upper = purrr::map_dbl(d_z_ci, 2)
	) %>%
	select(-d_z_ci)

# Plot SHAP importance (Fig 2)
shap_plot <- shap_summary_boot %>%
	
	ggplot(aes(x = weighted_shap, y = trait, fill = category)) +
	
	# Violin plot with proper alignment
	geom_violin(aes(color = category), alpha = 0.3, scale = "width", 
							position = position_identity(), size = 0.3) + 
	
	# Median points aligned exactly at the same y-axis position as violins
	geom_point(data = shap_effects, 
						 aes(x = median_successional, y = trait, color = "Successional Filtering", shape = "Successional Filtering"), 
						 size = 3, inherit.aes = FALSE) +
	geom_point(data = shap_effects, 
						 aes(x = median_environmental, y = trait, color = "Environmental Filtering", shape = "Environmental Filtering"), 
						 size = 3, inherit.aes = FALSE) +
	
	# Custom colours for violins and points
	scale_fill_manual(values = c("Successional Filtering" = "forestgreen", "Environmental Filtering" = "grey50")) +  # Violin fill colors
	scale_color_manual(values = c("Successional Filtering" = "darkgreen", "Environmental Filtering" = "black")) +   # Point & violin outline colors
	scale_shape_manual(values = c("Successional Filtering" = 20, "Environmental Filtering" = 18)) +  # Point shapes
	
	theme_bw() +
	labs(x = "Total SHAP Importance (Weighted Normalization)",
			 y = NULL,
			 fill = "Predictor",
			 color = "Predictor",
			 shape = "Predictor") +
	
	theme(text = element_text(size = 11),
				legend.position = c(.83,.15),
				legend.background = element_rect(colour = "black", linewidth = .2))

shap_plot

# Export bootstrapped results
if (export) {

	write_csv(bootstrap_shap, file = paste0(path_out, "/output/data/bootstrap_shap.csv"))
	write_csv(shap_effects, paste0(path_out, "/output/data/shap_effects.csv"))
	
	ggsave(filename = paste0(path_out, "/output/plots/fig2.png"),
				 plot = shap_plot, 
				 bg = "transparent",
				 width = 170, 
				 height = 120, 
				 units = "mm", 
				 dpi = 800)
}

rm(bootstrap_shap, shap_summary_boot) # clean environment

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
	group_by(trait, variable, bootstrap_id) %>%
	summarise(total_shap = sum(abs(shap_value)), .groups = "drop") %>%
	
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

# Build the plot (Fig 3)
pdp_plot <- pdp_stats_ellipses %>% 
ggplot(aes(x = intercept_diff, y = slope_diff, fill = trait, shape = trait)) +
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
	facet_wrap(~trait, scales = "free_y", ncol = 4, nrow = 2) +
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


