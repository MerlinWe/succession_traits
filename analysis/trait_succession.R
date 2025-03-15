############################################################################################################################
########################################  MSc Diss. Forest Succession Data Analysis ########################################  
############################################################################################################################

rm(list = ls())   # make sure environment is clean 
set.seed(42)      # set seed for reproducibility

# ----- Session set-up -----

# Load necessary libraries
library(caret)
library(ranger)
library(pdp)
library(mgcv)
library(rsample)
library(doParallel)
library(ggbeeswarm)
library(forcats)
library(dggridR)
library(grid)
library(glue)
library(gtable)
library(cowplot)
library(gridExtra)
library(psych)
library(factoextra)
library(tidyverse)

# Export plots?
export = FALSE

# Parallize? 32 cores if on threadripper - 8 if local
parallel = TRUE 

# Spatial downsampling for code development?
downsample = FALSE 
down_res   = 9

# Model tuning or predefined parameters?
tuning = FALSE

# Check which device is running
node_name <- Sys.info()["nodename"]

# Set file paths conditionally
path_in <- ifelse(node_name == "threadeast", "/home/merlin/RDS_drive/merlin/data/fia_traits", 
									"/Users/serpent/Documents/MSc/Thesis/Data/fia_traits")

path_out <- ifelse(node_name == "threadeast", "/home/merlin/traits_output", 
									 "/Users/serpent/Documents/MSc/Thesis/Code")

# Get functions
if (node_name == "threadeast") {
	source(url("https://raw.githubusercontent.com/MerlinWe/succession_traits/main/analysis/scenario_modeling/scenario_functions.R"))
} else {
	source("/Users/serpent/Documents/MSc/Thesis/Code/analysis/scenario_modeling_2.0/scenario_functions.R")
}

# Set parallel cluster 
if (parallel) { 
	num_cores <-  ifelse(node_name == "threadeast", 32, 8)
	cl <- makeCluster(num_cores)
	registerDoParallel(cl, cores = num_cores)
	getDoParWorkers()
}

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
	ggsave(paste0(path_out, "/output/plots/supplementary/s4.png"),
				 plot = pca_plot,
				 bg = "white",
				 width = 200,  
				 height = 130, 
				 units = "mm",
				 dpi = 600)
	write_csv(scores, file = paste0(path_out, "/output/data/pca_scores.csv"))
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

# ----- Spatial downsampling (if TRUE) -----

if (downsample) {
	grid_res <- down_res
	dggs <- dgconstruct(res = grid_res, metric = TRUE, resround = 'nearest')
	data <- data %>%
		mutate(cell = dgGEO_to_SEQNUM(dggs, in_lon_deg = LON, in_lat_deg = LAT)$seqnum) %>%
		dplyr::select(-LAT, -LON) %>%
		group_by(cell) %>%
		dplyr::summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = 'drop')
}

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
		ggsave(paste0(path_out, "/output/plots/supplementary/s5.png"),
					 plot = tuning_error_plot,
					 bg = "white",
					 width = 200,
					 height = 130,
					 units = "mm",
					 dpi = 600)
	}
	
	rm(tuning_error_plot, tuning_result)
	
} else { # Use a default hyper_grid with known best values
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
	write_csv(performance_metrics, file = paste0(path_out, "/output/data/performance_metrics.csv"))
}

rm(dggs, split, varimax_contributions)

## ---------- Calculate Shapley values for each trait model ----------

# Define number of bootstrap iterations
num_bootstraps <- 100

bootstrap_results <- foreach(i = 1:num_bootstraps, .combine = 'rbind', .packages = c('dplyr', 'ranger', 'fastshap', 'tidyverse', 'foreach', 'rsample')) %dopar% {
	bootstrap_shap(data, traits, covariates, hyper_grid, num_threads = num_cores, boot_id = i)
}

bootstrap_results <- read_csv("/Users/serpent/Documents/MSc/Thesis/Code/analysis/scenario_modeling_2.0/output/bootsrap_shap.csv")

shap_summary_boot <- bootstrap_results %>%
	mutate(category = case_when(
		variable == "standage" ~ "Temporal Succession",
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
		median_shap = median(shap_value),  # Use median instead of mean
		lower_ci = quantile(shap_value, probs = 0.025),  # Bootstrapped confidence interval
		upper_ci = quantile(shap_value, probs = 0.975),
		.groups = "drop")

# Plot results
shap_plot <- ggplot(shap_summary_boot, aes(x = median_shap, y = trait, color = category)) +
	geom_point(size = 2) +
	geom_errorbar(aes(y = trait, xmin = lower_ci, xmax = upper_ci, color = category), width = 0.2) +
	scale_y_discrete(labels = c(
		"wood_density" = "Wood Density",
		"bark_thickness" = "Bark Thickness",
		"conduit_diam" = "Conduit Diameter",
		"leaf_n" = "Leaf Nitrogen",
		"specific_leaf_area" = "Specific Leaf Area",
		"seed_dry_mass" = "Seed Dry Mass",
		"shade_tolerance" = "Shade Tolerance",
		"height" = "Tree Height")) +
	scale_color_manual(values = c("Temporal Succession" = "black", "Environmental Filtering" = "darkred")) +
	theme_clean() +
	labs( x = "Total SHAP Importance (Weighted Normalization)",
				y = NULL,
				color = "Predictor",
				caption = "Points represent bootstrapped median values. Error bars indicate 95% confidence intervals from bootstrapping") +
	theme(text = element_text(size = 12), legend.position = "top")

if (export) {
	write_csv(bootstrap_results, file = paste0(path_out, "/bootsrap_shap.csv"))
	write_csv(shap_summary_boot, paste0(path_out, "/shap_summary_boot.csv"))
	
	ggsave(filename = "/Users/serpent/Desktop/shap_plot.png",
				 plot = shap_plot, 
				 bg = "transparent",
				 width = 220, 
				 height = 120, 
				 units = "mm", 
				 dpi = 800)
}








