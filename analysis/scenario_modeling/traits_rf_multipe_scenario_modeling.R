############################################################################################################################
########################################  MSc Diss. Forest Succession Data Analysis ########################################  
############################################################################################################################

rm(list = ls()); gc()   # make sure environment is clean 
set.seed(42)            # set seed for reproducibility

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
export <- TRUE

# Parallize?
parallel <- TRUE 

# Spatial downsampling? If TRUE set resolution!
downsample <- TRUE 
down_res <- 9

# Model tuning or predefined parameters?
tuning = TRUE

# Check which device is running
node_name <- Sys.info()["nodename"]

# Set file paths conditionally
path_in <- ifelse(node_name == "threadeast", 
									"/home/merlin/RDS_drive/merlin/data/fia_traits", 
									"/Volumes/ritd-ag-project-rd01pr-dmayn10/merlin/data/fia_traits")

path_out <- ifelse(node_name == "threadeast", 
									 "/home/merlin/traits_output", 
									 "/Users/serpent/Documents/MSc/Thesis/Code/analysis")

if (parallel) {
	
	# For parallel processing: Register cores - 32 if on threadripper; 8 if local
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
filter_scores <- function(scores, loadings, threshold = threshold) {
	
	strong_loadings <- abs(loadings) >= threshold
	filtered_scores <- matrix(NA, nrow = nrow(scores), ncol = ncol(scores))
	for (i in 1:ncol(loadings)) {
		
		filtered_scores[, i] <- ifelse(strong_loadings[, i], scores[, i], NA)
		
	}
	
	rownames(filtered_scores) <- rownames(scores)
	colnames(filtered_scores) <- colnames(scores)
	return(filtered_scores)
	
}

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

# Function to create contribution plots with abbreviated x-axis labels
create_contrib_plot <- function(var_contrib_data, title, fill_color, abbreviations) {
	var_contrib_data %>%
		as_tibble() %>%
		mutate(name = abbreviations[match(name, names(abbreviations))]) %>%
		ggplot(aes(x = reorder(name, -contrib), y = contrib)) +
		geom_bar(stat = "identity", color = "black", fill = fill_color, alpha = 0.8) +
		labs(y = "Contributions (%)", x = NULL, title = title) +
		theme_bw() +
		theme(text = element_text(family = "sans", size = 8),
					axis.text.x = element_text(angle = 45, hjust = 1))
}

# Create data for varimax rotation contributions
varimax_contrib <- function(loadings) {
	loadings^2 / rowSums(loadings^2) * 100
}

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
	create_contrib_plot("PC1 - Variable Contributions after Varimax Rotation", "firebrick4", abbreviations)
var_contrib_pc2 <- tibble(name = rownames(varimax_contributions), contrib = varimax_contributions[, 2]) %>%
	create_contrib_plot("PC2 - Variable Contributions after Varimax Rotation", "tan4", abbreviations)
var_contrib_pc3 <- tibble(name = rownames(varimax_contributions), contrib = varimax_contributions[, 3]) %>%
	create_contrib_plot("PC3 - Variable Contributions after Varimax Rotation", "dodgerblue3", abbreviations)

# Combine the plots
pca_plot <- plot_grid(scree_plot, var_contrib_pc1, var_contrib_pc2, var_contrib_pc3,
											ncol = 2, nrow = 2, labels = "auto", label_fontfamily = "sans", label_size = 8)

if (export) {
	ggsave(paste0(path_out, "/plots/appendix/pca_plot.png"),
				 plot = pca_plot,
				 bg = "white",
				 width = 200,  
				 height = 130, 
				 units = "mm",
				 dpi = 600)
	
	write_csv(scores, file = paste0(path_out, "/tables/pca_scores.csv"))
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

# Use spatial downsampling to speed up runtime and to deal with spatial autocorrelation? 
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

# Define function to perform cross-validation and find the best hyperparameters
tune_rf_model <- function(trait, data, covariates, hyper_grid, num_threads = 1) {
	
	# Define model formula
	formula <- as.formula(paste(trait, "~", paste(covariates, collapse = " + ")))
	
	# Perform tuning
	prediction_error <- foreach(
		num.trees = hyper_grid$num.trees,
		mtry = hyper_grid$mtry,
		min.node.size = hyper_grid$min.node.size,
		.combine = 'c', 
		.packages = "ranger"
		
	) %dopar% {
		
		# Fit models
		mod <- ranger::ranger(
			formula = formula,
			data = data,
			dependent.variable.name = trait,
			num.trees = num.trees,
			mtry = mtry,
			min.node.size = min.node.size,
			num.threads = num_threads)
		
		# Return prediction error
		return(mod$prediction.error)
	}
	
	# Plot prediction errors
	error_plot <- hyper_grid %>% 
		mutate(prediction_error = prediction_error) %>%
		ggplot(aes(x = mtry, y = as.factor(min.node.size), fill = prediction_error)) + 
		facet_wrap(~ num.trees) + 
		geom_tile() + 
		scale_y_discrete(breaks = c(1, 10, 20)) + 
		scale_fill_viridis_c() + 
		ylab("min.node.size") + 
		ggtitle(trait) + 
		theme(text = element_text(family = "sans", size = 6),
					legend.key.height = unit(3, "mm"),
					legend.key.width = unit(3, "mm"),
					plot.title = element_text(face = "bold", size = 6))
	
	# Get best set 
	best_hyperparameters <- hyper_grid %>% 
		mutate(trait = trait) %>%
		arrange(prediction_error) %>% 
		slice(1)
	
	return(list(error_plot = error_plot, best_hyperparameters = best_hyperparameters))
}

# Define function to fit a random forest with best hyperparameters per trait
fit_rf_model <- function(trait, data, covariates, hyper_parameters, num_threads = 1) {
	
	# Define model formula for each trait
	formula <- as.formula(paste(trait, "~", paste(covariates, collapse = " + ")))
	
	# Get trait-specific hyperparameters
	hyper_grid <- hyper_parameters[hyper_parameters$trait == trait, ]
	
	# Fit model
	trait_mod <- ranger::ranger(
		formula = formula,
		data = data,
		num.trees = hyper_grid$num.trees[1],
		mtry = hyper_grid$mtry[1],
		min.node.size = hyper_grid$min.node.size[1],
		num.threads = num_threads) 
	
	# Get performance metrics 
	performance <- data.frame(
		trait = trait,
		mtry = trait_mod$mtry,
		num_trees = trait_mod$num.trees,
		min_node_size = trait_mod$min.node.size,
		rsq = trait_mod$r.squared,
		pred_error = trait_mod$prediction.error)
	
	# Return model and performance metrics 
	return(list(trait_mod = trait_mod, performance = performance))
}

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
		ggsave(paste0(path_out, "/plots/appendix/tuning_plot.png"),
					 plot = tuning_error_plot,
					 bg = "white",
					 width = 200,
					 height = 130,
					 units = "mm",
					 dpi = 600)
	}
	
	rm(tuning_error_plot, tuning_result)
	
} else {
	
	# Use a default hyper_grid with predefined values
	hyper_grid <- tibble(
		trait = traits,
		num.trees = rep(500, length(traits)),
		mtry = rep(4, length(traits)),
		min.node.size = rep(1, length(traits)))
}

# Now fit trait models using the full dataset
best_models <- map(traits, ~ fit_rf_model(.x, data, covariates, hyper_grid, num_threads = ifelse(node_name == "threadeast", 4, 1)))
performance_metrics <- map(best_models, "performance") %>% bind_rows()
best_models <- map(best_models, "trait_mod")
names(best_models) <- traits

rm(dggs, split, varimax_contributions)

## ---------- Calculate Shapley values for each trait model ----------

# Prediction function for fastshap
predict_fn <- function(object, newdata) {
	predict(object, data = newdata)$prediction
}

# Loop through final trait models to calculate shapley values 
shap_values_list <- list()

# Loop over each model and trait with inner parallelization
for (i in seq_along(best_models)) {
	model <- best_models[[i]]
	trait <- traits[i]
	
	shap_values <- fastshap::explain(model, 
																	 X = test_data %>% select(all_of(covariates)) %>% as.matrix(),
																	 pred_wrapper = predict_fn, 
																	 nsim = 10, 
																	 parallel = TRUE)  
	
	shap_values <- as.data.frame(shap_values)
	shap_values$trait <- trait
	
	# Append to the list
	shap_values_list[[i]] <- shap_values
}

# Combine all results into a single data frame
shap_values <- bind_rows(shap_values_list) %>% as_tibble()
	
## Visualize shapley values; first get a tibble for plotting
shap_long <- shap_values %>%
	pivot_longer(cols = -c(trait), names_to = "feature", values_to = "shap") %>%
	as_tibble() %>%
	mutate(feature = gsub("_", " ", feature),
				 trait_name = case_when(
				 	trait == "wood_density" ~ "Wood Density",
				 	trait == "bark_thickness" ~ "Bark Thickness",
				 	trait == "conduit_diam" ~ "Conduit Diameter",
				 	trait == "leaf_n" ~ "Leaf Nitrogen",
				 	trait == "specific_leaf_area" ~ "Specific Leaf Area",
				 	trait == "seed_dry_mass" ~ "Seed Dry Mass",
				 	trait == "shade_tolerance" ~ "Shade Tolernce",
				 	trait == "height" ~ "Tree Height",
				 	TRUE ~ NA_character_)) 

# Sum absolute Shapley values to determine overall importance
feature_importance <- shap_long %>%
	group_by(trait, feature) %>%
	summarize(importance = sum(abs(shap)), .groups = "drop") %>%
	arrange(trait, importance) %>%
	right_join(shap_long, by = c("trait", "feature")) %>%
	left_join(
		performance_metrics %>%
			mutate(trait_name = case_when(
				trait == "wood_density" ~ "Wood Density",
				trait == "bark_thickness" ~ "Bark Thickness",
				trait == "conduit_diam" ~ "Conduit Diameter",
				trait == "leaf_n" ~ "Leaf Nitrogen",
				trait == "specific_leaf_area" ~ "Specific Leaf Area",
				trait == "seed_dry_mass" ~ "Seed Dry Mass",
				trait == "shade_tolerance" ~ "Shade Tolernce",
				trait == "height" ~ "Tree Height",
				TRUE ~ NA_character_)) %>%
			dplyr::select(trait, rsq) %>%
			mutate(rsq = round(rsq, digits = 3)) %>%
			distinct(),
		by = "trait") 

# Shorten feature labels
feature_labels <- c(
	"standage" = "STAGE",
	"temp pc" = "TEMP (PC)",
	"soil pc" = "SOIL (PC)",
	"rain pc"= "PRCP (PC)",
	"elevation" = "ELEV",
	"soil ph" = "pH",
	"biome temperate conifer forests" = "B-TC",
	"biome temperate broadleaf forests" = "B-TB",
	"biome temperate grasslands" = "B-TG",
	"biome xeric shrublands" = "B-XS",
	"biome mediterranean woodlands" = "B-MW",
	"biome boreal forests or taiga" = "B-BF",
	"biome tundra" = "B-TU",
	"biome flooded grasslands" = "B-FG")

# Function to create shapley plots for a single trait
plot_shapley_for_trait <- function(trait_id, hide_x_axis_labels = FALSE) {
	
	feature_importance %>%
		filter(trait == trait_id) %>%
		{
			shap_range <- shap_max - shap_min
			y_lim_min <- min(.$shap) - 0.1 * shap_range
			y_lim_max <- max(.$shap) + 0.1 * shap_range
			
			# Plot shapley bee swarms
			p1 <- ggplot(., aes(x = reorder(feature, importance), y = shap, color = shap)) +
				geom_quasirandom(alpha = 0.5) +
				scale_color_viridis_c(option = "viridis",guide = "none") +
				scale_x_discrete(labels = feature_labels) +
				scale_y_continuous(n.breaks = 5, limits = c(y_lim_min, y_lim_max)) +
				coord_flip() +
				labs(x = NULL, y = "Shapley Value", title = glue("{distinct(., trait_name)$trait_name}")) +
				theme_bw(base_line_size = .3, base_rect_size = .5) +
				theme(text = element_text(family = "sans", size = 6),
							axis.title.x = element_text(family = "sans", size = 5),
							plot.title = element_text(face = "bold", size = 6),
							strip.text = element_text(face = "bold"),
							plot.margin = margin(t=7.5, b=5, r=1, l=5)) +
				theme(axis.title.x = if(hide_x_axis_labels) element_blank() else element_text())
							
			# Plot absolute shapley importance 
			p2 <- ggplot(distinct(., feature, importance, rsq), aes(x = reorder(feature, importance), y = importance)) +
				geom_bar(stat = "identity", fill = "#1B7E74", color = "black", alpha = 0.7, linewidth = 0.2) +
				scale_y_continuous(n.breaks = 3) +
				labs(x = NULL, y = "Feat. Importance", title = glue("R² = {distinct(., rsq)$rsq}")) +
				coord_flip() +
				theme_bw(base_line_size = .3, base_rect_size = .5) +
				theme(text = element_text(family = "sans", size = 6),
							plot.title = element_text(size = 6),
							axis.text.y = element_blank(),
							axis.ticks.y = element_blank(),
							axis.title.x = element_text(family = "sans", size = 5),
							plot.margin = margin(t = 15, b = 2.5, r = 5, l = 1)) +
				theme(axis.title.x = if(hide_x_axis_labels) element_blank() else element_text())
			
			# Set plot structure 
			plot_grid(p1, p2, 
								ncol = 2, 
								nrow = 1, 
								rel_widths = c(.7, .3),
								align = "h", 
								axis = c("l"), 
								greedy = TRUE)
		}
}

# Add legend and build compound figure 
shap_min <- min(feature_importance$shap)
shap_max <- max(feature_importance$shap)
shap_range <- shap_max - shap_min

legend <- get_legend(
	feature_importance %>%
		filter(trait_name == "Bark Thickness") %>%
		ggplot(aes(x = reorder(feature, importance), y = shap, color = shap)) +
		geom_quasirandom(alpha = 0.5) +
		scale_color_viridis_c(option = "viridis",
													name = "Feature\nValue",
													breaks = c(shap_min + 0.4 * shap_range,
																		 shap_max - 0.15 * shap_range),
													labels = c("low", "high"),
													guide = guide_colorbar(
														label.position = "right",
														barwidth = unit(3, "mm"),
														barheight = unit(100, "mm"))) +
		theme(legend.position = "right", 
					text = element_text(family = "sans", size = 6),
					legend.title = element_text(size = 6),
					legend.text = element_text(size = 6),
					legend.margin = margin(t=0, b=-5, r=0, l=0)))

# Create a list of plots with a flag to hide x-axis labels for the upper row
plots <- map2(traits, rep(c(TRUE, FALSE), each = length(traits) / 2), plot_shapley_for_trait)

shap_plot <- plot_grid(
	plot_grid(
		plotlist = plots, 
		ncol = 4, nrow = 2, align = "v"),
	legend, ncol = 2, rel_widths = c(1, .07))

if (export) {
	ggsave(paste0(path_out, "/plots/shap_plot.png"),
				 plot = shap_plot,
				 bg = "white",
				 width = 200,  
				 height = 130, 
				 units = "mm",
				 dpi = 600)
	
	write.csv(shap_values, file = paste0(path_out, "/tables/shap_values.csv"))
}

# Clean memory before parallel processing
rm(feature_importance, legend, model, shap_long, plots, shap_max, shap_min, shap_range,
	 shap_plot, shap_values, shap_values_list, feature_labels, i, trait); gc()

rm(best_models, test_data, train_data)

## -------- Calculate partial dependence curves ----------

# Here, we stratify the data based on 25% quantiles of our environmental variables and fit new 
# models for every scenario.Thereby we make sure to only consider scenarios observed in the data,
# and implicitly consider interaction between covariates and biomes.

## >>>>> Stratify data based on variable quantiles <<<<<

# Define data stratification
stratify <- function(data, variable, quantiles) {
	
	lower_quantile <- quantile(data[[variable]], quantiles[1])
	upper_quantile <- quantile(data[[variable]], quantiles[2])
	
	lower_data <- data %>% filter(.[[variable]] <= lower_quantile)
	upper_data <- data %>% filter(.[[variable]] >= upper_quantile)
	
	return(list(lower = lower_data, upper = upper_data))
}

# Function to stratify, fit models, and extract performance metrics
fit_models_on_strata <- function(data, variable, traits, covariates, hyper_grid, num_threads, labels) {
	
	stratified_data <- stratify(data, variable, quantiles_25)
	
	lower_models <- map(traits, ~ fit_rf_model(.x, stratified_data$lower, covariates, hyper_grid, num_threads))
	names(lower_models) <- traits
	
	upper_models <- map(traits, ~ fit_rf_model(.x, stratified_data$upper, covariates, hyper_grid, num_threads))
	names(upper_models) <- traits
	
	lower_performance <- map(lower_models, "performance") %>% bind_rows() %>% mutate(group = labels[1])
	upper_performance <- map(upper_models, "performance") %>% bind_rows() %>% mutate(group = labels[2])
	
	performance_metrics <- bind_rows(lower_performance, upper_performance)
	
	return(list(models = list(lower = lower_models, upper = upper_models), 
							performance = performance_metrics,
							stratified_data = stratified_data))
}

# Set Quantile values; num_threads; and quantile labels
quantiles_25 <- c(0.25, 0.75)

num_threads <- ifelse(node_name == "threadeast", 4, 1)

quantile_labels <- list(
	temp_pc = c("Cold Temperatures (Lower 25%)", "Warm Temperatures (Upper 25%)"),
	soil_pc = c("Sandy Soils (Lower 25%)", "Water-Retentive Soils (Upper 25%)"),
	rain_pc = c("Low Precipitation (Lower 25%)", "High Precipitation (Upper 25%)"),
	elevation = c("Low Elevation (Lower 25%)", "High Elevation (Upper 25%)"),
	soil_ph = c("Low Soil pH (Lower 25%)", "High Soil pH (Upper 25%)"))

# Perform stratification, model training, and performance extraction for each variable
temp <- fit_models_on_strata(data, "temp_pc", traits, covariates, hyper_grid, num_threads, quantile_labels$temp_pc)
soil <- fit_models_on_strata(data, "soil_pc", traits, covariates, hyper_grid, num_threads, quantile_labels$soil_pc)
prcp <- fit_models_on_strata(data, "rain_pc", traits, covariates, hyper_grid, num_threads, quantile_labels$rain_pc)
elev <- fit_models_on_strata(data, "elevation", traits, covariates, hyper_grid, num_threads, quantile_labels$elevation)
ph <- fit_models_on_strata(data, "soil_ph", traits, covariates, hyper_grid, num_threads, quantile_labels$soil_ph)

# Get r squared values
performance <- bind_rows(temp$performance, soil$performance, prcp$performance,
												 elev$performance, ph$performance) %>%
	dplyr::select(trait, group, rsq)

# Function to calculate partial dependence and confidence intervals
calculate_partial_dependence <- function(model, data, feature, conf.level = 0.95) {
	
	partial_results <- partial(model, pred.var = feature, train = data, parallel = FALSE)
	
	# Create a template for predictions
	pred_data <- data
	
	# Calculate predictions for each grid point
	predictions <- foreach(i = seq_len(nrow(partial_results)), .combine = rbind, .packages = "ranger") %dopar% {
		pred_data[[feature]] <- partial_results[i, feature]
		predict(model, data = pred_data)$predictions
	}
	
	# Calculate mean and confidence intervals
	partial_results$yhat_mean <- rowMeans(predictions)
	partial_results$yhat_lower <- apply(predictions, 1, function(x) quantile(x, probs = (1 - conf.level) / 2))
	partial_results$yhat_upper <- apply(predictions, 1, function(x) quantile(x, probs = 1 - (1 - conf.level) / 2))
	
	return(partial_results)
}

# Calculate partial dependence for each trait and scenario with group labels
calculate_pdp_for_scenario <- function(scenario_data, traits, feature, labels) {
	results <- foreach(trait = traits, .combine = rbind, .packages = c('ranger', 'pdp', 'dplyr', 'foreach'), 
										 .export = c('calculate_partial_dependence')) %dopar% {
										 	
										 	model_lower <- scenario_data$models$lower[[trait]][["trait_mod"]]
										 	model_upper <- scenario_data$models$upper[[trait]][["trait_mod"]]
										 	
										 	data_lower <- scenario_data$stratified_data$lower
										 	data_upper <- scenario_data$stratified_data$upper
										 	
										 	# Diagnostic prints
										 	print(paste("Trait:", trait))
										 	print("Model lower class:")
										 	print(class(model_lower))
										 	print("Model upper class:")
										 	print(class(model_upper))
										 	print("Data lower class:")
										 	print(class(data_lower))
										 	print("Data upper class:")
										 	print(class(data_upper))
										 	
										 	pdp_lower <- calculate_partial_dependence(model_lower, data_lower, feature)
										 	pdp_upper <- calculate_partial_dependence(model_upper, data_upper, feature)
										 	
										 	pdp_lower$group <- labels[1]
										 	pdp_upper$group <- labels[2]
										 	pdp_lower$trait <- trait
										 	pdp_upper$trait <- trait
										 	
										 	bind_rows(pdp_lower, pdp_upper)
										 }
	return(results)
}

# Quantile labels
quantile_labels <- list(
	temp_pc = c("Cold Temperatures (Lower 25%)", "Warm Temperatures (Upper 25%)"),
	soil_pc = c("Sandy Soils (Lower 25%)", "Water-Retentive Soils (Upper 25%)"),
	rain_pc = c("Low Precipitation (Lower 25%)", "High Precipitation (Upper 25%)"),
	elevation = c("Low Elevation (Lower 25%)", "High Elevation (Upper 25%)"),
	soil_ph = c("Low Soil pH (Lower 25%)", "High Soil pH (Upper 25%)"))

# Perform PDP calculation for each scenario
pdp_temp <- calculate_pdp_for_scenario(temp, traits, "standage", quantile_labels$temp_pc) %>% mutate(variable = "Temperature (PC)")
pdp_soil <- calculate_pdp_for_scenario(soil, traits, "standage", quantile_labels$soil_pc) %>% mutate(variable = "Soil - Water Retention (PC)")
pdp_prcp <- calculate_pdp_for_scenario(prcp, traits, "standage", quantile_labels$rain_pc) %>% mutate(variable = "Precipitation (PC)")
pdp_elev <- calculate_pdp_for_scenario(elev, traits, "standage", quantile_labels$elevation) %>% mutate(variable = "Elevation")
pdp_ph <- calculate_pdp_for_scenario(ph, traits, "standage", quantile_labels$soil_ph) %>% mutate(variable = "Soil pH")

# Combine all results
pdp_data <- bind_rows(pdp_temp, pdp_soil, pdp_prcp, pdp_elev, pdp_ph) %>%
	left_join(performance, by = c("trait", "group"))

if (export) {
	write_csv(pdp_data, file = paste0(path_out, "/tables/pdp_data.csv"))
}

rm(pdp_temp, pdp_soil, pdp_prcp, pdp_elev, pdp_ph, elev, performance, ph, prcp, soil, temp); gc()

## >>> Plot pdp curves <<<

# Function for pdp plotting
create_pdp_plot <- function(data, levels, colors, labels, x_lab, y_lab, show_y_strip_labels) {
	data %>%
		mutate(group = factor(group, levels = levels)) %>%
		ggplot(aes(x = standage, y = yhat, color = group, shape = group)) +
		geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, linewidth = .5) +
		scale_color_manual(values = setNames(colors, labels)) +
		scale_fill_manual(values = setNames(colors, labels)) +
		labs(x = x_lab, y = y_lab, color = "Quantile", shape = "Quantile", fill = "Quantile") +
		theme_bw(base_line_size = .3, base_rect_size = .5) +
		theme(
			legend.position = "none",
			text = element_text(family = "sans", size = 8),
			strip.background.y = if (show_y_strip_labels) element_rect(fill = "white", color = "black", linewidth = .75) else element_blank(),
			strip.text.y = if (show_y_strip_labels) element_text() else element_blank()) +
		facet_grid(trait ~ group, scales = "free", labeller = labeller(trait = trait_labels))
}

# Define plot parameters and generate plots
plot_params <- list(
	list(data = pdp_data %>% filter(variable == "Temperature (PC)"),
			 levels = quantile_labels$temp_pc,
			 colors = c("darkslateblue", "darkred"),
			 labels = quantile_labels$temp_pc,
			 x_lab = NULL,
			 y_lab = "Predicted Trait Value (log-scaled)",
			 show_y_strip_labels = FALSE),
	
	list(data = pdp_data %>% filter(variable == "Precipitation (PC)"),
			 levels = quantile_labels$rain_pc,
			 colors = c("tan4", "cyan4"),
			 labels = quantile_labels$rain_pc,
			 x_lab = NULL,
			 y_lab = NULL,
			 show_y_strip_labels = TRUE),
	
	list(data = pdp_data %>% filter(variable == "Soil - Water Retention (PC)"),
			 levels = quantile_labels$soil_pc,
			 colors = c("khaki4", "springgreen3"),
			 labels = quantile_labels$soil_pc,
			 x_lab = NULL,
			 y_lab = "Predicted Trait Value (log-scaled)",
			 show_y_strip_labels = FALSE),
	
	list(data = pdp_data %>% filter(variable == "Elevation"),
			 levels = quantile_labels$elevation,
			 colors = c("bisque4", "royalblue4"),
			 labels = quantile_labels$elevation,
			 x_lab = "Standage",
			 y_lab = NULL,
			 show_y_strip_labels = TRUE),
	
	list(data = pdp_data %>% filter(variable == "Soil pH"),
			 levels = quantile_labels$soil_ph,
			 colors = c("orangered1", "darkorchid4"),
			 labels = quantile_labels$soil_ph,
			 x_lab = "Standage",
			 y_lab = "Predicted Trait Value (log-scaled)",
			 show_y_strip_labels = TRUE))

# Clean labels for traits
trait_labels <- c(
	wood_density = "Wood\nDens.",
	bark_thickness = "Bark\nThick.",
	conduit_diam = "Conduit\nDiam.",
	leaf_n = "Leaf\nNitro.",
	specific_leaf_area = "SLA",
	seed_dry_mass = "Seed\nDry\nMass",
	shade_tolerance = "Shade\nTol.",
	height = "Tree\nHeight")

# Quantile levels for table 
quantile_levels <- c("Cold Temperatures", "Warm Temperatures",
										 "Low Precipitation", "High Precipitation",
										 "Sandy Soils", "Water-Retentive Soils",
										 "Low Elevation", "High Elevation",
										 "Low Soil pH", "High Soil pH")

# Generate PDP plots
pdp_plots <- lapply(plot_params, function(params) {
	create_pdp_plot(params$data, params$levels, params$colors, params$labels, params$x_lab, params$y_lab, params$show_y_strip_labels)
})

# Create summary table with R² values
r2_summary <- pdp_data %>%
	group_by(trait, group) %>%
	summarize(rsq = round(mean(rsq), 3)) %>%
	mutate(trait = case_when(
		trait == "wood_density" ~ "Wood\nDensity",
		trait == "bark_thickness" ~ "Bark\nThickness",
		trait == "conduit_diam" ~ "Conduit\nDiameter",
		trait == "leaf_n" ~ "Leaf\nNitrogen",
		trait == "specific_leaf_area" ~ "SLA",
		trait == "seed_dry_mass" ~ "Seed\nDry Mass",
		trait == "shade_tolerance" ~ "Shade\nTolerance",
		trait == "height" ~ "Tree\nHeight",
		TRUE ~ NA_character_)) %>%
	spread(trait, rsq) %>%
	mutate(group = factor(group, levels = unique(group))) %>%
	mutate(group = factor(trimws(gsub("\\s*\\([^\\)]+\\)", "", group)), levels = quantile_levels)) %>%
	arrange(group) %>%
	rename(Quantile = group)

# Split r2_summary into two sub-tables
r2_summary_1 <- r2_summary[, c("Quantile", "Wood\nDensity", "Bark\nThickness", "Conduit\nDiameter", "Leaf\nNitrogen")]
r2_summary_2 <- r2_summary[, c("Quantile", "SLA", "Seed\nDry Mass", "Shade\nTolerance", "Tree\nHeight")]

# Define a custom theme for the table grob with left alignment and minimal padding
custom_theme <- ttheme_minimal(
	core = list(fg_params = list(hjust = 0, x = 0, fontsize = 5, just = "left"), bg_params = list(col = NA)),
	colhead = list(fg_params = list(hjust = 0, x = 0, fontsize = 5, just = "left")),
	padding = unit(c(2, 2), "mm"))

# Convert summary tables to compressed table grobs with custom themes
table_grob_1 <- tableGrob(r2_summary_1, rows = NULL, theme = custom_theme)
table_grob_2 <- tableGrob(r2_summary_2, rows = NULL, theme = custom_theme)

# Add horizontal lines to the tables
add_lines <- function(grob) {
	grob <- gtable_add_grob(grob,
													grobs = segmentsGrob( # Top horizontal line
														x0 = unit(0, "npc"),
														y0 = unit(1, "npc"),
														x1 = unit(1, "npc"),
														y1 = unit(1, "npc"),
														gp = gpar(lwd = 0.5)),
													t = 1, b = 1, l = 1, r = ncol(grob))
	grob <- gtable_add_grob(grob,
													grobs = segmentsGrob( # Bottom horizontal line
														x0 = unit(0, "npc"),
														y0 = unit(0, "npc"),
														x1 = unit(1, "npc"),
														y1 = unit(0, "npc"),
														gp = gpar(lwd = 0.5)),
													t = nrow(grob), b = nrow(grob), l = 1, r = ncol(grob))
	grob
}

table_grob_1 <- add_lines(table_grob_1)
table_grob_2 <- add_lines(table_grob_2)

# Create a title grob
title_grob <- textGrob("R² of Quantile Models", gp = gpar(fontsize = 7, fontface = "bold", hjust = 0, x = 0))

# Combine title and table grobs with decreased spacing
table_grob_combined <- arrangeGrob(
	title_grob,
	table_grob_1,
	table_grob_2,
	ncol = 1,
	heights = unit.c(unit(1.5, "lines"), unit(1, "null"), unit(1.2, "null")))

# Convert table grob with title to ggplot object
table_plot <- ggplot() +
	annotation_custom(table_grob_combined) +
	theme_void()

# Add the table to the list of PDP plots
pdp_plots <- c(pdp_plots, list(table_plot))

# Build compound figure with tables in the empty spots
pdp_25 <- plot_grid(
	plotlist = pdp_plots,
	ncol = 2,
	nrow = 3,
	rel_widths = c(1, 1),
	align = "h",
	axis = c("l"),
	greedy = TRUE,
	labels = "auto",
	label_fontfamily = "sans",
	label_size = 8)


# Save the plot if export is TRUE
if (export) {
	ggsave(filename = paste0(path_out, "/plots/pdp_quantiles.png"),
				 plot = pdp_25, 
				 bg = "white",
				 width = 180, 
				 height = 290, 
				 units = "mm", 
				 dpi = 600)
	write.csv(pdp_data, file = paste0(path_out, "/tables/pdp_25.csv"))
}

## Done with parallel processing - stop the cluster
stopCluster(cl)
rm(plot_params, pdp_plots, variables); gc()

## ---------- Conclusive analysis ----------

# Function to process PDP data
get_pdp_delta <- function(pdp_data) {
	delta <- pdp_data %>%
		group_by(trait, group) %>%
		do({
			fit <- lm(yhat ~ standage, data = .)
			intercept <- coef(fit)[1]
			slope <- coef(fit)[2]
			max_standage <- max(.$standage)
			y_max_standage <- intercept + slope * max_standage
			delta <- y_max_standage - intercept
			tibble(intercept = intercept, slope = slope, sd_yhat = sd(.$yhat), delta = delta, r_squared = summary(fit)$r.squared)
		}) %>%
		ungroup()
	return(delta)
}

# Summarise delta for traits
delta <- get_pdp_delta(pdp_data) 

# ------ Succession Figures -----

# For manuscript

succession <- delta %>%
	group_by(trait) %>%
	summarize(
		mean_slope = mean(slope),
		sd_slope = sd(slope),
		mean_delta = mean(delta),
		sd_delta = sd(delta),
		mean_rsq = mean(r_squared),
		sd_rsq = sd(r_squared)) %>%
	ungroup() %>%
	mutate(trait = case_when(
		trait == "wood_density" ~ "Wood Density",
		trait == "bark_thickness" ~ "Bark Thickness",
		trait == "conduit_diam" ~ "Conduit Diameter",
		trait == "leaf_n" ~ "Leaf Nitrogen",
		trait == "specific_leaf_area" ~ "Specific Leaf Area",
		trait == "seed_dry_mass" ~ "Seed Dry Mass",
		trait == "shade_tolerance" ~ "Shade Tolerance",
		trait == "height" ~ "Tree Height",
		TRUE ~ NA_character_)) 

succession <- succession %>%
	ggplot(aes(x = mean_delta, y = mean_rsq, color = trait, shape = trait, label = trait)) +
	geom_point(size = 3, fill = "white") +
	geom_errorbar(aes(ymin = mean_rsq - sd_rsq, ymax = mean_rsq + sd_rsq), linewidth = .2 ,width = 0.05) +
	geom_errorbar(aes(xmin = mean_delta - sd_delta, xmax = mean_delta + sd_delta), linewidth = .2, width = 0.05) +
	scale_colour_viridis_d() +
	scale_shape_manual(values = c(15, 16, 17, 18, 19, 15, 16, 17)) +
	labs(x = "Mean Difference Across Succession (log-scaled)", y = "Mean R²",
			 color = "Trait", shape = "Trait") +
	geom_vline(xintercept = median(succession$mean_delta), linetype = "dashed", color = "black", linewidth = .5) +
	geom_hline(yintercept = median(succession$mean_rsq), linetype = "dashed", color = "black", linewidth = .5) +
	theme_bw() +
	theme(legend.position = "right",
				text = element_text(family = "sans", size = 8))

print(succession) # check plot

# For appendix 
succession_append <- delta %>%
	mutate(group2 = factor(group, levels = c("Cold Temperatures (Lower 25%)", "Warm Temperatures (Upper 25%)",
																					"Low Precipitation (Lower 25%)", "High Precipitation (Upper 25%)",
																					"Sandy Soils (Lower 25%)", "Water-Retentive Soils (Upper 25%)",
																					"Low Elevation (Lower 25%)", "High Elevation (Upper 25%)",
																					"Low Soil pH (Lower 25%)", "High soil pH (Upper 25%)"))) %>%
	mutate(trait = case_when(
		trait == "wood_density" ~ "Wood Density",
		trait == "bark_thickness" ~ "Bark Thickness",
		trait == "conduit_diam" ~ "Conduit Diameter",
		trait == "leaf_n" ~ "Leaf Nitrogen",
		trait == "specific_leaf_area" ~ "Specific Leaf Area",
		trait == "seed_dry_mass" ~ "Seed Dry Mass",
		trait == "shade_tolerance" ~ "Shade Tolerance",
		trait == "height" ~ "Tree Height",
		TRUE ~ NA_character_))  %>%
	ggplot(aes(x = delta, y = r_squared, color = trait, shape = trait, label = trait)) +
	geom_point(size = 3, fill = "white") +
	scale_colour_viridis_d() +
	scale_shape_manual(values = c(15, 16, 17, 18, 19, 15, 16, 17)) +
	labs(x = "Difference Across Succession (log-scaled)", y = "R²",
			 color = "Trait", shape = "Trait") +
	theme_bw() +
	theme(legend.position = c(.75,.15),
				legend.direction = "horizontal",
				legend.title = element_blank(),
				legend.background = element_rect(fill = "transparent", colour = "transparent"),
				text = element_text(family = "sans", size = 8),
				strip.background = element_rect(fill = "white", color = "black", linewidth = .75)) +
	guides(color = guide_legend(ncol = 2, byrow = TRUE), shape = guide_legend(ncol = 2, byrow = TRUE)) +
	facet_wrap(~group, ncol = 4, nrow = 3)

if (export) {
	ggsave(filename = paste0(path_out, "/plots/succession_plot.png"),
				 plot = succession, 
				 bg = "white",
				 width = 130, 
				 height = 100, 
				 units = "mm", 
				 dpi = 600)
	
	ggsave(filename = paste0(path_out, "/plots/appendix/quantiles_succession.png"),
				 plot = succession_append, 
				 bg = "white",
				 width = 200, 
				 height = 180, 
				 units = "mm", 
				 dpi = 600)
	
	write.csv(delta, file = paste0(path_out, "/tables/succession.csv") )
}

rm(pdp_25, trait_labels, custom_theme, delta, hyper_grid, performance_metrics, quantile_labels, title_grob,
	 r2_summary, r2_summary_1, r2_summary_2, table_grob_1, table_grob_2, table_grob_combined, table_plot, cl); gc()

## Done! 
