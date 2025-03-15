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
downsample = TRUE 
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

if (parallel) { # set cluster 
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

# Visualize shapley values; first get a tibble for plotting
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
	"rain pc" = "PRCP (PC)",
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
													breaks = c(shap_min + 0.31 * shap_range,
																		 shap_max - 0.08 * shap_range),
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
		plotlist = shap_plots, 
		ncol = 4, nrow = 2, align = "v"),
	legend, ncol = 2, rel_widths = c(1, .07))

# Clean memory before parallel processing
rm(feature_importance, legend, model, shap_long, shap_plots, shap_max, shap_min, shap_range,
	 shap_values, shap_values_list, feature_labels, i, trait)

## -------- Calculate partial dependence curves ----------

## Calculate partial dependence for global trait models
pdp_traits <- map_dfr(names(best_models), function(trait) {
	model <- best_models[[trait]]
	result <- calculate_partial_dependence(model, data = data, feature = "standage")
	result <- result %>%
		mutate(trait = trait) 
	return(result)
})

# Plot trait pdp curves
trait_curves <- pdp_traits %>%
	group_by(trait) %>%
	mutate(
		signal = last(yhat) - first(yhat),
		trait_label = case_when(
			trait == "wood_density" ~ paste("Wood Density", "(Signal =", round(signal, 2), ")"),
			trait == "bark_thickness" ~ paste("Bark Thickness", "(Signal =", round(signal, 2), ")"),
			trait == "conduit_diam" ~ paste("Conduit Diameter", "(Signal =", round(signal, 2), ")"),
			trait == "leaf_n" ~ paste("Leaf Nitrogen", "(Signal =", round(signal, 2), ")"),
			trait == "specific_leaf_area" ~ paste("Specific Leaf Area", "(Signal =", round(signal, 2), ")"),
			trait == "seed_dry_mass" ~ paste("Seed Dry Mass", "(Signal =", round(signal, 2), ")"),
			trait == "shade_tolerance" ~ paste("Shade Tolerance", "(Signal =", round(signal, 2), ")"),
			trait == "height" ~ paste("Tree Height", "(Signal =", round(signal, 2), ")"),
			TRUE ~ NA_character_
		)
	) %>%
	ungroup() %>%
	ggplot(aes(x = standage, y = yhat, color = trait_label)) + 
	geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, linewidth = .7) +
	scale_fill_viridis_d() +
	scale_colour_viridis_d() +
	#scale_x_continuous(expand = c(0, 0)) +
	labs(x = "Standage (years)", y = "Predicted Trait Value (log-scaled)", color = "Trait", fill = "Trait") +
	theme_bw(base_line_size = .3, base_rect_size = .5) +
	facet_wrap(~trait_label, scales = "fixed", nrow = 2, ncol = 4) + 
	theme(
		legend.position = "none",
		text = element_text(family = "sans", size = 8),
		strip.background.x = element_rect(fill = "white", color = "black", linewidth = .75),
		strip.text.x = element_text(face = "bold"),
		plot.margin = unit(c(1,1.4,1,.3), "cm"))

# Build Fig2 compound figure
trait_plot <- plot_grid(
	shap_plot,
	trait_curves,
	nrow=2, ncol = 1,
	rel_heights = c(1, .6),
	labels = c("a.", "b."))

if (export) {
	ggsave(paste0(path_out, "/output/plots/fig2.png"),
				 plot = trait_plot,
				 bg = "white",
				 width = 200,  
				 height = 220, 
				 units = "mm",
				 dpi = 600)
	write_csv(shap_values, file = paste0(path_out, "/output/data/shap_values.csv"))
	write_csv(pdp_traits,  file = paste0(path_out, "/output/data/pdp_traits.csv"))
	
}

rm(best_models, test_data, train_data, trait_signal, signal_plot, pdp_append)

## >>>>> Scenario modelling <<<<<

# Next, we stratify the data based on 25% quantiles of our environmental variables and fit new 
# models for every scenario.Thereby we make sure to only consider scenarios observed in the data,
# and implicitly consider interaction between covariates and biomes.

# Here, we stratify the data based on 25% quantiles of our environmental variables and fit new 
# models for every scenario.Thereby we make sure to only consider scenarios observed in the data,
# and implicitly consider interaction between covariates and biomes.

## >>>>> Stratify data based on variable quantiles <<<<<

# Set Quantile values; covariates, num_threads; and quantile labels
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

## ----- Fit models calculate R squared again but without standage -----

# We also check for a change in model performance when the variable standage is not considered. Thus, 
# if trait expressions are ONLY explained by the environment (and biomes). 

covariates_ns <- covariates[covariates != "standage"]

# Fit models in parallel using foreach
temp_ns <- fit_models_on_strata(data, "temp_pc", traits, covariates_ns, hyper_grid, num_threads, quantile_labels$temp_pc)
soil_ns <- fit_models_on_strata(data, "soil_pc", traits, covariates_ns, hyper_grid, num_threads, quantile_labels$soil_pc)
prcp_ns <- fit_models_on_strata(data, "rain_pc", traits, covariates_ns, hyper_grid, num_threads, quantile_labels$rain_pc)
elev_ns <- fit_models_on_strata(data, "elevation", traits, covariates_ns, hyper_grid, num_threads, quantile_labels$elevation)
ph_ns <- fit_models_on_strata(data, "soil_ph", traits, covariates_ns, hyper_grid, num_threads, quantile_labels$soil_ph)

# Get r squared values and calculate difference due to standage
performance_ns <- bind_rows(temp_ns$performance, soil_ns$performance, prcp_ns$performance, elev_ns$performance, ph_ns$performance) %>%
	dplyr::select(trait, group, rsq) %>%
	rename(rsq_ns = rsq) %>%
	left_join(performance, by = c("trait", "group")) %>%
	mutate(rsq_diff = rsq_ns-rsq)

## >>>>>>>>>> Calculate Partial Dependence <<<<<<<<<<

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

rm(pdp_temp, pdp_soil, pdp_prcp, pdp_elev, pdp_ph, elev, ph, prcp, soil, temp)

## >>> Plot pdp curves <<<

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
			 show_y_strip_labels = FALSE),
	
	list(data = pdp_data %>% filter(variable == "Soil - Water Retention (PC)"),
			 levels = quantile_labels$soil_pc,
			 colors = c("goldenrod4", "springgreen3"),
			 labels = quantile_labels$soil_pc,
			 x_lab = NULL,
			 y_lab = "Predicted Trait Value (log-scaled)",
			 show_y_strip_labels = TRUE),
	
	list(data = pdp_data %>% filter(variable == "Elevation"),
			 levels = quantile_labels$elevation,
			 colors = c("darkolivegreen", "cadetblue3"),
			 labels = quantile_labels$elevation,
			 x_lab = "Standage",
			 y_lab = NULL,
			 show_y_strip_labels = FALSE),
	
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

# Quantile levels for barplot 
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
		trait == "wood_density" ~ "Wood Density",
		trait == "bark_thickness" ~ "Bark Thickness",
		trait == "conduit_diam" ~ "Conduit Diameter",
		trait == "leaf_n" ~ "Leaf Nitrogen",
		trait == "specific_leaf_area" ~ "Specific Leaf Area",
		trait == "seed_dry_mass" ~ "Seed Dry Mass",
		trait == "shade_tolerance" ~ "Shade Tolerance",
		trait == "height" ~ "Tree Height",
		TRUE ~ NA_character_)) %>%
	mutate(group = factor(group, levels = unique(group))) %>%
	mutate(group = factor(trimws(gsub("\\s*\\([^\\)]+\\)", "", group)), levels = quantile_levels)) %>%
	arrange(group) %>%
	rename(Quantile = group)

# Extract colors from plot parameters for bar plot
color_mapping <- unlist(lapply(plot_params, function(x) setNames(x$colors, x$quantile_levels)))
color_mapping <- color_mapping[unique(r2_summary$Quantile)]

# Create the bar plot for R² values
r2_bar_plot <- r2_summary %>%
	mutate(Quantile = factor(trimws(gsub("\\s*\\([^\\)]+\\)", "", Quantile)), levels = quantile_levels)) %>%
	ggplot(aes(x = Quantile, y = rsq, fill = Quantile)) +
	geom_bar(stat = "identity", position = "dodge", alpha = .5, colour = "black") +
	scale_fill_manual(values = color_mapping) +
	facet_wrap(~trait, nrow = 2, ncol = 4) +
	labs(x = NULL, y = "R²") +
	theme_bw() +
	theme(text = element_text(family = "sans", size = 7),
				strip.background =  element_rect(fill = "white", color = "black", linewidth = .75),
				strip.text.x = element_text(face = "bold"),
				axis.text.x = element_text(angle = 45, hjust = 1, size = 4),
				legend.position = "none")


# Add the bar plot to the list of PDP plots
pdp_plots <- c(pdp_plots, list(r2_bar_plot))

# Build compound figure with the bar plot in the last position
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
	ggsave(filename = paste0(path_out, "/output/plots/supplementary/s6.png"),
				 plot = pdp_25, 
				 bg = "white",
				 width = 210, 
				 height = 290, 
				 units = "mm", 
				 dpi = 600)
}

## >>> Plot pdp summary <<<

# Here, we summarise the difference between feature scenarios for the trait intercepts and signal to understand 
# environmental filtering versus modulation. 

pdp_summary <- pdp_data %>%
	select(trait, standage, yhat, variable, group) %>%
	mutate(group = ifelse(grepl("Upper 25%", group), "Upper", "Lower")) %>%
	
	# Calculate the intercepts (standage == 0) and highest standage for each trait, variable, and group
	group_by(trait, variable, group) %>%
	summarize(
		intercept = yhat[which.min(standage)],
		delta = yhat[which.max(standage)] - yhat[which.min(standage)]) %>%
	
	# Calculate the difference between upper and lower groups for intercepts and deltas
	pivot_wider(names_from = group, values_from = c(intercept, delta)) %>%
	mutate(intercept_diff = intercept_Upper - intercept_Lower,
				 delta_diff = delta_Upper - delta_Lower) %>%
	select(trait, variable, intercept_diff, delta_diff) %>% 
	
	# Make sure values are absolute
	mutate_if(is.numeric, abs) %>%
	
	# Set some labels
	mutate(trait = case_when(
		trait == "wood_density" ~ "Wood Density",
		trait == "bark_thickness" ~ "Bark Thickness",
		trait == "conduit_diam" ~ "Conduit Diameter",
		trait == "leaf_n" ~ "Leaf Nitrogen",
		trait == "specific_leaf_area" ~ "Specific Leaf Area",
		trait == "seed_dry_mass" ~ "Seed Dry Mass",
		trait == "shade_tolerance" ~ "Shade Tolerance",
		trait == "height" ~ "Tree Height",
		TRUE ~ NA_character_)) %>% 
	
	mutate(variable = case_when(
		variable == "Elevation" ~ "Elevation\n ",
		variable == "Precipitation (PC)" ~ "Precipitation\n (PC)",
		variable == "Soil - Water Retention (PC)" ~ "Soil Water\n Retention (PC)",
		variable == "Soil pH" ~ "Soil pH\n ",
		variable == "Temperature (PC)" ~ "Temperature\n (PC)",
		TRUE ~ NA_character_)) %>%
	
	# Set factor levels
	mutate(variable = factor(variable, levels = unique(variable)))

# Calculate summarise for features 
feature_summary <- pdp_summary %>%
	group_by(variable) %>%
	mutate_if(is.numeric, abs) %>%
	summarize(
		mean_intercept_diff = mean(intercept_diff, na.rm = TRUE),
		se_intercept_diff = sd(intercept_diff, na.rm = TRUE) / sqrt(n()),
		mean_delta_diff = mean(delta_diff, na.rm = TRUE),
		se_delta_diff = sd(delta_diff, na.rm = TRUE) / sqrt(n())) %>%
	mutate_if(is.numeric, abs) %>%
	mutate(variable = factor(variable, levels = unique(variable)))

# Calculate summarise for traits 
trait_summary <- pdp_summary %>%
	group_by(trait) %>%
	summarise(
		mean_delta_signal = mean(delta_diff),
		se_delta_signal = sd(delta_diff) / sqrt(n()),
		
		mean_delta_intercept = mean(intercept_diff),
		se_delta_intercept = sd(intercept_diff) / sqrt(n()))

## --- Build compound figure ---

# Feature labels
g.mid <- ggplot(pdp_summary, aes(x=1, y=variable)) +
	geom_text(aes(label=variable), hjust = 0.5, size = 4.5, family = "sans") +
	theme_void() +
	theme(axis.title=element_blank(),
				panel.grid=element_blank(),
				axis.text=element_blank(),
				axis.ticks=element_blank(),
				plot.margin = unit(c(1,-1,1.5,-1), "cm"))

# Difference in Intercept
p1 <- ggplot(pdp_summary, aes(x = variable, y = intercept_diff, fill = trait)) +
	geom_bar(stat = "identity", position = "dodge", color = "black", alpha = .8) +
	scale_fill_viridis_d() +
	coord_flip() +
	ylim(c(0, .9)) +
	theme_bw(base_rect_size = 1) +
	theme(text = element_text(family = "sans", size = 14),
				axis.title.y = element_blank(),
				axis.text.y = element_blank(),
				axis.ticks.y = element_blank(),
				plot.margin = unit(c(1, 0, 1, 1), "cm"),
				legend.direction = "vertical",
				legend.position = c(.225, .77),
				legend.background = element_rect(fill = "white", colour = "black", linewidth = .3)) +
	labs(y = expression(Delta~"Intercept (absolute & log-scaled)"), fill = NULL) +  
	scale_y_reverse() +
	coord_flip()

# Difference in Delta
p2 <- ggplot(pdp_summary, aes(x = variable, y = delta_diff, fill = trait)) +
	geom_bar(stat = "identity", position = "dodge", color = "black", alpha = .8) +
	scale_fill_viridis_d() +
	coord_flip() +
	ylim(c(0, .9)) +
	theme_bw(base_rect_size = 1) +
	theme(text = element_text(family = "sans", size = 14),
				axis.title.y = element_blank(),
				axis.text.y = element_blank(),
				axis.ticks.y = element_blank(),
				plot.margin = unit(c(1, 1, 1, 0), "cm"),
				legend.position = "none") +
	labs(y = expression(Delta~"Signal (absolute & log-scaled)"))

# Summary plot
variable_colors <- c("Temperature\n (PC)" = "darkred", 
										 "Soil pH\n " = "purple", 
										 "Soil Water\n Retention (PC)" = "orange", 
										 "Precipitation\n (PC)" = "navyblue", 
										 "Elevation\n " = "skyblue")

# Feature summary
p_all <- ggplot(feature_summary, aes(x = variable)) +
	geom_point(aes(y = mean_intercept_diff, color = variable, shape = "Intercept Difference"), 
						 size = 5, 
						 position = position_nudge(x = 0.2)) +
	geom_errorbar(aes(ymin = mean_intercept_diff - se_intercept_diff, 
										ymax = mean_intercept_diff + se_intercept_diff, 
										color = variable), 
								width = 0.2, 
								position = position_nudge(x = 0.2)) +
	geom_point(aes(y = mean_delta_diff, color = variable, shape = "Signal Difference"), 
						 size = 5, 
						 position = position_nudge(x = -0.2)) +
	geom_errorbar(aes(ymin = mean_delta_diff - se_delta_diff, 
										ymax = mean_delta_diff + se_delta_diff, 
										color = variable), 
								width = 0.2, 
								position = position_nudge(x = -0.2)) +
	
	coord_flip() +
	labs(x = NULL, y = "Absolute Difference (log-scale)") + 
	theme_bw(base_rect_size = 1) +
	scale_color_manual(values = variable_colors, guide = "none") + 
	scale_shape_manual(
		values = c("Intercept Difference" = 16, "Signal Difference" = 17),
		labels = c(expression(Delta~"Intercept"), expression(Delta~"Signal"))) + 
	theme(legend.position = c(.92, .87),
				legend.direction = "vertical",
				legend.background = element_rect(fill = "white", colour = "black", linewidth = .3),
				legend.title = element_blank(),
				text = element_text(family = "sans", size = 14),
				plot.margin = unit(c(.2,1,0,1), "cm"))

# Trait summary
p_traits <- trait_summary %>%	ggplot() +
	geom_point(aes(x = mean_delta_intercept, y = mean_delta_signal, group = trait, color = trait, shape = trait), 
						 size = 5, fill = "white") +
	
	geom_errorbar(aes(x = mean_delta_intercept, y = mean_delta_signal,
										xmin = mean_delta_intercept - se_delta_intercept, xmax = mean_delta_intercept + se_delta_intercept,
										group = trait, color = trait)) +
	geom_errorbar(aes(x = mean_delta_intercept, y = mean_delta_signal, 
										ymin = mean_delta_signal - se_delta_signal, ymax = mean_delta_signal + se_delta_signal,
										group = trait, color = trait)) +
	
	scale_colour_viridis_d() +
	
	scale_shape_manual(values = c(15, 16, 17, 18, 19, 15, 16, 17)) +
	labs(y = expression(Delta~"Signal (Abiotic Modulation)"),
			 x = expression(Delta~"Intercept (Abiotic Filtering)"),
			 fill = "Traits") + 
	
	theme_bw(base_rect_size = 1) +
	theme(legend.position = c(.9, .7),
				legend.direction = "vertical",
				legend.background = element_rect(fill = "white", colour = "black", linewidth = .3),
				legend.title = element_blank(),
				text = element_text(family = "sans", size = 14),
				plot.margin = unit(c(.2,1,0,1), "cm"))

# Build figure 
pdp_summary_plot <- plot_grid(p_all, 
															plot_grid(p1, g.mid, p2, 
																				ncol = 3, 
																				rel_widths = c(4/9, 1.2/9, 4/9),
																				labels = c("b.", " ", "c.")),
															p_traits,
															ncol = 1, nrow = 3, rel_heights = c(.65,1,.7), 
															labels = c("a.", " ", "d."))

# Save the plot if export is TRUE
if (export) {
	ggsave(filename = paste0(path_out, "/output/plots/fig3.png"),
				 plot = pdp_summary_plot, 
				 bg = "white",
				 width = 290, 
				 height = 390, 
				 units = "mm", 
				 dpi = 800)
}

## Done with parallel processing - stop the cluster
stopCluster(cl)
rm(plot_params, pdp_plots, g.mid, p1, p2, variable_colors, p_all); gc()

## ---------- Conclusive analysis ----------

# Calculate signal and intercept difference and display again r square
delta <- pdp_data %>% 
	
	get_pdp_delta() %>%
	
	select(trait, variable, group, intercept, r_squared) %>%
	group_by(trait, variable) %>%
	summarise(mean_rsq = mean(r_squared)) %>%
	ungroup() %>%
	
	# Set some labels
	mutate(trait = case_when(
		trait == "wood_density" ~ "Wood Density",
		trait == "bark_thickness" ~ "Bark Thickness",
		trait == "conduit_diam" ~ "Conduit Diameter",
		trait == "leaf_n" ~ "Leaf Nitrogen",
		trait == "specific_leaf_area" ~ "Specific Leaf Area",
		trait == "seed_dry_mass" ~ "Seed Dry Mass",
		trait == "shade_tolerance" ~ "Shade Tolerance",
		trait == "height" ~ "Tree Height",
		TRUE ~ NA_character_)) %>% 
	
	mutate(variable = case_when(
		variable == "Elevation" ~ "Elevation\n ",
		variable == "Precipitation (PC)" ~ "Precipitation\n (PC)",
		variable == "Soil - Water Retention (PC)" ~ "Soil Water\n Retention (PC)",
		variable == "Soil pH" ~ "Soil pH\n ",
		variable == "Temperature (PC)" ~ "Temperature\n (PC)",
		TRUE ~ NA_character_)) %>%
	
	left_join(pdp_summary, by = c("trait", "variable")) %>%
	rename(signal_diff = delta_diff) %>%
	
	group_by(trait) %>%
	summarise(mean_delta_intercept = mean(intercept_diff),
						se_delta_intercept = sd(intercept_diff) / sqrt(n()),
						mean_delta_signal = mean(signal_diff),
						se_delta_signal = sd(signal_diff) / sqrt(n()),
						trait_mean_rsq = mean(mean_rsq),
						trait_se_rsq = sd(mean_rsq) / sqrt(n())) %>%
	ungroup()


predictability <- delta %>%
	pivot_longer(cols = c(mean_delta_intercept, mean_delta_signal),
							 names_to = "metric", values_to = "value") %>%
	mutate(metric = recode(metric,
												 "mean_delta_intercept" = "Successional Intercept (log-scaled)",
												 "mean_delta_signal" = "Successional Signal (log-scaled)")) %>%
	ggplot(aes(x = value, y = trait_mean_rsq, color = trait, shape = trait, label = trait)) +
	geom_point(size = 3, fill = "white") +
	geom_errorbar(aes(ymin = trait_mean_rsq - trait_se_rsq, ymax = trait_mean_rsq + trait_se_rsq)) +
	geom_errorbar(aes(xmin = value - if_else(metric == "Successional Intercept (log-scaled)", se_delta_intercept, se_delta_signal),
										xmax = value + if_else(metric == "Successional Intercept (log-scaled)", se_delta_intercept, se_delta_signal))) +
	scale_colour_viridis_d() +
	scale_shape_manual(values = c(15, 16, 17, 18, 19, 15, 16, 17)) +
	labs(y = "Predictability (Mean R²)", x = "Mean Difference Between Feature Scenarios", color = "Trait", shape = "Trait") +
	theme_bw(base_line_size = .3) +
	xlim(c(0, .51)) +
	theme(
		legend.position = "top",
		legend.direction = "horizontal",
		legend.title = element_blank(),
		legend.spacing.x = unit(0.5, 'cm'),
		text = element_text(family = "sans", size = 12),
		strip.background = element_rect(fill = "white", color = "black", linewidth = .75)) +
	facet_grid(. ~ metric)

if (export) {
	ggsave(filename = paste0(path_out, "/output/plots/fig4.png"),
				 plot = predictability, 
				 bg = "white",
				 width = 220, 
				 height = 140, 
				 units = "mm", 
				 dpi = 600)
	
	write_csv(delta, file = paste0(path_out, "/output/data/delta.csv"))
	write_csv(performance_ns, file = paste0(path_out, "/output/data/performance_no_standage.csv"))
}

## ----- Other Supplementary Material -----


delta_intercept_append <- pdp_data %>% 
	
	get_pdp_delta() %>%
	
	select(trait, variable, group, intercept, r_squared) %>%
	group_by(trait, variable) %>%
	summarise(mean_rsq = mean(r_squared)) %>%
	ungroup() %>%
	
	# Set some labels
	mutate(trait = case_when(
		trait == "wood_density" ~ "Wood Density",
		trait == "bark_thickness" ~ "Bark Thickness",
		trait == "conduit_diam" ~ "Conduit Diameter",
		trait == "leaf_n" ~ "Leaf Nitrogen",
		trait == "specific_leaf_area" ~ "Specific Leaf Area",
		trait == "seed_dry_mass" ~ "Seed Dry Mass",
		trait == "shade_tolerance" ~ "Shade Tolerance",
		trait == "height" ~ "Tree Height",
		TRUE ~ NA_character_)) %>% 
	
	mutate(variable = case_when(
		variable == "Elevation" ~ "Elevation\n ",
		variable == "Precipitation (PC)" ~ "Precipitation\n (PC)",
		variable == "Soil - Water Retention (PC)" ~ "Soil Water\n Retention (PC)",
		variable == "Soil pH" ~ "Soil pH\n ",
		variable == "Temperature (PC)" ~ "Temperature\n (PC)",
		TRUE ~ NA_character_)) %>%
	
	left_join(pdp_summary, by = c("trait", "variable")) %>%
	rename(signal_diff = delta_diff) %>%
	
	ggplot(aes(x = intercept_diff, y = mean_rsq, color = trait, shape = trait, label = trait)) +
	geom_point(size = 3, fill = "white") +
	scale_colour_viridis_d() +
	scale_shape_manual(values = c(15, 16, 17, 18, 19, 15, 16, 17)) +
	labs(y = "Predictability (Mean R²)", x = "Intercept Difference in Feature Scenario", color = "Trait", shape = "Trait") +
	theme_bw(base_line_size = .3) +
	theme(legend.position = c(.85,.2),
				legend.direction = "horizontal",
				legend.title = element_blank(),
				legend.background = element_rect(fill = "transparent", colour = "transparent"),
				text = element_text(family = "sans", size = 12),
				strip.background = element_rect(fill = "white", color = "black", linewidth = .75)) +
	guides(color = guide_legend(ncol = 1, byrow = TRUE), shape = guide_legend(ncol = 2, byrow = TRUE)) +
	facet_wrap( ~ variable)

delta_signal_append <- pdp_data %>% 
	
	get_pdp_delta() %>%
	
	select(trait, variable, group, intercept, r_squared) %>%
	group_by(trait, variable) %>%
	summarise(mean_rsq = mean(r_squared)) %>%
	ungroup() %>%
	
	# Set some labels
	mutate(trait = case_when(
		trait == "wood_density" ~ "Wood Density",
		trait == "bark_thickness" ~ "Bark Thickness",
		trait == "conduit_diam" ~ "Conduit Diameter",
		trait == "leaf_n" ~ "Leaf Nitrogen",
		trait == "specific_leaf_area" ~ "Specific Leaf Area",
		trait == "seed_dry_mass" ~ "Seed Dry Mass",
		trait == "shade_tolerance" ~ "Shade Tolerance",
		trait == "height" ~ "Tree Height",
		TRUE ~ NA_character_)) %>% 
	
	mutate(variable = case_when(
		variable == "Elevation" ~ "Elevation\n ",
		variable == "Precipitation (PC)" ~ "Precipitation\n (PC)",
		variable == "Soil - Water Retention (PC)" ~ "Soil Water\n Retention (PC)",
		variable == "Soil pH" ~ "Soil pH\n ",
		variable == "Temperature (PC)" ~ "Temperature\n (PC)",
		TRUE ~ NA_character_)) %>%
	
	left_join(pdp_summary, by = c("trait", "variable")) %>%
	rename(signal_diff = delta_diff) %>%
	
	ggplot(aes(x = signal_diff, y = mean_rsq, color = trait, shape = trait, label = trait)) +
	geom_point(size = 3, fill = "white") +
	scale_colour_viridis_d() +
	scale_shape_manual(values = c(15, 16, 17, 18, 19, 15, 16, 17)) +
	labs(y = "Predictability (Mean R²)", x = "Signal Difference in Feature Scenario", color = "Trait", shape = "Trait") +
	theme_bw(base_line_size = .3) +
	theme(legend.position = c(.85,.2),
				legend.direction = "horizontal",
				legend.title = element_blank(),
				legend.background = element_rect(fill = "transparent", colour = "transparent"),
				text = element_text(family = "sans", size = 12),
				strip.background = element_rect(fill = "white", color = "black", linewidth = .75)) +
	guides(color = guide_legend(ncol = 1, byrow = TRUE), shape = guide_legend(ncol = 2, byrow = TRUE)) +
	facet_wrap( ~ variable)

if (export) {
	
	ggsave(filename = paste0(path_out, "/output/plots/supplementary/s7.png"),
				 plot = delta_intercept_append, 
				 bg = "white",
				 width = 200, 
				 height = 160, 
				 units = "mm", 
				 dpi = 600)
	
	ggsave(filename = paste0(path_out, "/output/plots/supplementary/s8.png"),
				 plot = delta_signal_append, 
				 bg = "white",
				 width = 200, 
				 height = 160, 
				 units = "mm", 
				 dpi = 600)
}

rm(pdp_25, trait_labels, custom_theme, delta, hyper_grid, performance_metrics, quantile_labels, title_grob,
	 r2_summary, r2_summary_1, r2_summary_2, table_grob_1, table_grob_2, table_grob_combined, table_plot, cl); gc()

## Done! 