############################################################################################################################
########################################  MSc Diss. Forest Succession Data Analysis ########################################  
############################################################################################################################

rm(list = ls()) # make sure environment is clean 
set.seed(42)    # set seed for reproducibility

# Load necessary libraries
library(caret)
library(ranger)
library(pdp)
library(mgcv)
library(rsample)
library(doParallel)
library(ggbeeswarm)
library(patchwork)
library(dggridR)
library(grid)
library(viridis)
library(tidyverse)

# For parallel processing: Register cores - 32 if on threadripper; 10 if local
num_cores <- 10
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Set path to run either on local device 
path_in <- "/Volumes/ritd-ag-project-rd01pr-dmayn10/merlin/data/fia_traits" 
path_out <- "/Users/serpent/Documents/MSc/Thesis/Code/analysis/traits" 

# Set paths to run on threadripper
path_in <- "/home/merlin/RDS_drive/merlin/data/fia_traits"
path_out <- "/home/merlin/traits_output" 

# ---------- Read and prepare input data ----------

# Read data and make some minor adjustments
data <- read_csv(paste0(path_in, "/plotlevel_data_2024-07-15.csv")) %>%
	
	# Keep only columns relevant for the traits analysis 
	select(starts_with("wmean_"), standage, biome, managed,
				 annual_mean_temperature, annual_precipitation, temperature_seasonality, mean_diurnal_range,
				 min_temperature_of_coldest_month, max_temperature_of_warmest_month,
				 elevation, pop_density, sand_content_015cm, soil_ph_015cm, water_capacity_015cm, LAT, LON) %>%
	filter(complete.cases(.)) %>%
	
	# Adjust some column names 
	rename_with(~ gsub("wmean_", "", .), starts_with("wmean_")) %>%
	rename_with(~ gsub("_015cm", "", .), ends_with("_015cm")) %>%
	rename(min_temperature = min_temperature_of_coldest_month,
				 max_temperature = max_temperature_of_warmest_month) %>%
	
	# Filter standage by the upper 10% quantiles
	filter(standage < quantile(standage, 0.9)) %>%
	# Disregard plots that have been managed 
	filter(managed == 0) %>%
	select(-managed)

# Create binary dummy variables based on biome levels with >100 observations; drop data of remaining biomes 
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

# Use spatial downsampling to speed up runtime and assess spatial autocorrelation? 
grid_res <- 12
dggs <- dgconstruct(res = grid_res, metric = TRUE, resround = 'nearest') 
data <- data %>% 
	mutate(cell = dgGEO_to_SEQNUM(dggs, in_lon_deg = LON, in_lat_deg = LAT)$seqnum) %>%
	dplyr::select(-LAT, -LON) %>%
	group_by(cell) %>%
	dplyr::summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = 'drop')

## ---------- Model training and validation ----------

# Split into training and test sets
split <- initial_split(data, prop = 0.8)

# Extract training and testing datasets
train_data <- training(split)
test_data  <- testing(split)

# Define traits
traits <- c("wood_density", "bark_thickness", "conduit_diam", "leaf_n", 
						"specific_leaf_area", "seed_dry_mass", "shade_tolerance", "height")

# Define covariates 
covariates <- c("standage", "annual_mean_temperature", "annual_precipitation", 
								"temperature_seasonality", "mean_diurnal_range", "min_temperature", 
								"max_temperature", "elevation", "pop_density", "sand_content", "soil_ph", 
								"water_capacity", "biome_boreal_forests_or_taiga", "biome_flooded_grasslands",        
								"biome_mediterranean_woodlands", "biome_temperate_broadleaf_forests", "biome_tundra",
								"biome_temperate_conifer_forests", "biome_temperate_grasslands", "biome_xeric_shrublands")

# Define a grid of hyperparameters including mtry and min.node.size
hyper_grid <- expand.grid(
	mtry = seq(2, length(covariates), by = 2),
	min.node.size = c(1, 5, 10),
	splitrule = "variance"
)

# Define function to perform cross-validation, find the best hyperparameters, and capture performance metrics
tune_rf_model <- function(trait, data, covariates, hyper_grid) {
	
	# Define model formula
	formula <- as.formula(paste(trait, "~", paste(covariates, collapse = " + ")))
	
	# Set training framework
	train_control <- trainControl(
		method = "cv", 
		number = 5, 
		savePredictions = "final", 
		summaryFunction = defaultSummary)
	
	# Train models
	model <- train(
		formula, 
		data = data,
		method = "ranger",
		trControl = train_control,
		tuneGrid = hyper_grid,
		importance = "permutation")
	
	# Get performance metrics
	results <- model$results
	best_tune <- model$bestTune
	results <- results %>%
		mutate(
			trait = trait, 
			best = rowSums(sapply(names(best_tune), 
														function(param) results[[param]] == best_tune[[param]])) == length(best_tune))
	
	return(list(model = model, results = results))
}

# Perform hyperparameter tuning for each trait and get performance metrics
tuned_models <- traits %>%
	map(~ tune_rf_model(.x, train_data, covariates, hyper_grid))

# Extract the best model for every trait
best_models <- tuned_models %>%
	map(~ .x$model$finalModel)

# Extract and write performance metrics
global_performance_metrics <- tuned_models %>%
	map_df(~ .x$results) %>%
	mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
	arrange(trait, best) %>%
	write_csv(file = paste0(path_out, "/tables/global_rf_performance_metrics.csv"))

## ---------- Calculate Shapley values for each trait model ----------

# Prediction function for fastshap
predict_fn <- function(object, newdata) {
	predict(object, data = newdata)$predictions
}

shap_values <- foreach(i = seq_along(best_models), .combine = 'rbind', .packages = c('fastshap', 'dplyr')) %dopar% {
	model <- best_models[[i]]
	trait <- traits[i]
	shap_values <- fastshap::explain(model, X = test_data %>% select(all_of(covariates)) %>% as.matrix(),
																	 pred_wrapper = predict_fn, nsim = 100, parallel = TRUE)
	shap_values <- as.data.frame(shap_values)
	shap_values$trait <- trait
	shap_values
}

## Visualize shapley values; first get a tibble for plotting
shap_long <- shap_values %>%
	pivot_longer(cols = -c(trait), names_to = "feature", values_to = "shap") %>%
	as_tibble() %>%
	mutate(feature = gsub("_", " ", feature),
				 trait = case_when(
				 	trait == "wood_density" ~ "Wood Density",
				 	trait == "bark_thickness" ~ "Bark Thickness",
				 	trait == "conduit_diam" ~ "Conduit Diameter",
				 	trait == "leaf_n" ~ "Leaf Nitrogen",
				 	trait == "specific_leaf_area" ~ "Specific Leaf Area",
				 	trait == "seed_dry_mass" ~ "Seed Dry Mass",
				 	trait == "shade_tolerance" ~ "Shade Tolernce",
				 	trait == "height" ~ "Tree Height",
				 	TRUE ~ NA_character_)) 

shap_plot <- shap_long %>%
	ggplot(aes(x = feature, y = shap, color = shap)) +
	geom_quasirandom(alpha = 0.5) +
	facet_wrap(~trait, ncol = 4, nrow = 2, scale = "fixed") +
	scale_color_viridis_c(option = "viridis",
												name = "Feature\nValue",
												breaks = c(min(shap_long$shap), max(shap_long$shap)),
												labels = c("low", "high")) +
	coord_flip() +
	labs(x = NULL, y = "Shapley Value") +
	theme_bw() +
	theme(text = element_text(family = "Arial"),
				legend.position = "right", 
				legend.key.width = unit(0.5, "cm"),
				legend.key.height = unit(2, "cm"), 
				legend.box.background = element_rect(color = "black", linewidth = .75), 
				strip.background = element_rect(fill = "white", color = "black", linewidth = .75),
				strip.text = element_text(color = "black")) 

ggsave(filename = paste0(path_out, "/plots/global_shap_plots.png"),
			 plot = shap_plot,
			 bg = "white",
			 width = 280,
			 height = 180,
			 units = "mm",
			 dpi = 1457)

# Sum absolute Shapley values to determine overall importance; create plot 
feature_importance <- shap_long %>%
	group_by(trait, feature) %>%
	summarize(importance = sum(abs(shap)), .groups = "drop") %>%
	arrange(trait, importance) %>%
	group_by(trait) %>%
	mutate(avg_importance = mean(importance)) %>%
	ungroup() 
feature_importance <- feature_importance %>%
	left_join(feature_importance %>%
		 	group_by(feature) %>%
		 	summarize(global_importance = sum(importance), .groups = "drop") %>%
		 	arrange(desc(global_importance)), by = "feature")


# Plot with reordered features based on global importance
importance_plot <- ggplot(feature_importance, aes(x = reorder(feature, global_importance), y = importance, fill = trait)) +
	geom_bar(stat = "identity", colour = "black", alpha = .7) +
	geom_hline(aes(yintercept = avg_importance), linetype = "dashed", colour = "red") +
	facet_wrap(~trait, scales = "fixed", ncol = 2, nrow = 4) +
	coord_flip() +
	scale_fill_viridis_d() +
	labs(x = NULL,
			 y = "Overall Feature Importance",
			 fill = "Trait",
			 color = "Average Importance") +
	theme_bw(base_size = 15) +
	theme(legend.position = "none",
				text = element_text(family = "Arial"),
				axis.text.y = element_text(size = 9))

ggsave(filename = paste0(path_out, "/plots/global_importance_plot.png"),
			 plot = importance_plot,
			 bg = "white",
			 width = 250,
			 height = 270,
			 units = "mm",
			 dpi = 1457)

## -------- Calculate partial dependence curves ----------

# General function to create partial dependence data for each trait based on quantiles
get_partial_dependence <- function(best_models, train_data, traits, pred.vars, quantiles, quantile_names = NULL) {
	
	# Get quantile values for covariate
	quantile_values <- quantile(train_data[[pred.vars[2]]], quantiles)
	
	if (is.null(quantile_names)) {
		quantile_names <- paste0("Quantile ", seq_along(quantiles))
	}
	
	# Print the threshold values for the quantiles
	for (i in seq_along(quantiles)) {
		cat(quantile_names[i], "threshold value:", quantile_values[i], "\n")
	}
	
	# Separate data for each quantile
	data_list <- lapply(seq_along(quantile_values), function(i) {
		q_value <- quantile_values[i]
		if (i == 1) {
			train_data %>% filter(!!sym(pred.vars[2]) <= q_value)
		} else {
			train_data %>% filter(!!sym(pred.vars[2]) >= q_value)
		}
	})
	
	pdp_data_combined <- traits %>%
		map_df(function(trait) {
			model <- best_models[[which(traits == trait)]]
			
			# Combine partial dependence data for each quantile
			pdp_data_list <- lapply(seq_along(quantile_values), function(i) {
				q_value <- quantile_values[i]
				data_subset <- data_list[[i]]
				
				# Build grid for partial dependence
				pred_grid <- expand.grid(
					standage = unique(data_subset$standage))
				pred_grid[[pred.vars[2]]] <- q_value
				
				# Calculate partial dependence
				pdp_data <- pdp::partial(
					object = model,
					pred.var = pred.vars,
					pred.grid = pred_grid,
					train = data_subset,
					plot = FALSE)
				
				# Label the data with the quantile name
				pdp_data %>% mutate(group = quantile_names[i])
			})
			
			pdp_data_combined <- bind_rows(pdp_data_list) %>%
				mutate(trait = trait,
							 lower = yhat - 1.96 * sd(yhat),  # Assuming normal distribution for 95% CI
							 upper = yhat + 1.96 * sd(yhat))
			
			return(pdp_data_combined)
		})
	
	return(pdp_data_combined)
}

##### Climate quantiles #####

# Clean labels for traits
trait_labels <- c(
	wood_density = "Wood Density",
	bark_thickness = "Bark Thickness",
	conduit_diam = "Conduit Diameter",
	leaf_n = "Leaf Nitrogen",
	specific_leaf_area = "Specific Leaf Area",
	seed_dry_mass = "Seed Dry Mass",
	shade_tolerance = "Shade Tolerance",
	height = "Tree Height"
)

# >>> 10% Quantiles for Temperature <<<

# Set Quantile values
quantiles_10 <- c(0.1, 0.9)
quantile_names_10 <- c("Coldest 10%", "Hottest 10%")
plot_title <- "Partial Dependence Plot for Annual Mean Temperature"

# Calculate partial dependence towards standage for all traits for the quantiles
pdp_temp_10 <- get_partial_dependence(
	best_models = best_models,
	train_data = train_data,
	traits = traits,
	pred.vars = c("standage", "annual_mean_temperature"),
	quantiles = quantiles_10,
	quantile_names = quantile_names_10)

# Plot the results with uncertainty ribbons and facet_grid
temp_10_plot <- ggplot(pdp_temp_10, aes(x = standage, y = yhat, color = group, shape = group)) +
	geom_ribbon(aes(ymin = lower, ymax = upper, fill = group), alpha = 0.2) +
	geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, linewidth = .7) +
	scale_color_manual(values = setNames(c("#0800af", "#c82300"), quantile_names_10)) +
	scale_fill_manual(values = setNames(c("#0800af", "#c82300"), quantile_names_10)) +
	labs(x = "Standage", y = "Predicted Trait Value (log-scaled)", color = "Quantile", shape = "Quantile", fill = "Quantile", title = plot_title) + 
	theme_bw() +
	theme(
		legend.position = "top",
		text = element_text(family = "Arial"),
		strip.text.y = element_text(size = 7),
		strip.background = element_rect(fill = "white", color = "black", linewidth = .75)) +
	facet_grid(trait ~ group, scales = "free", labeller = labeller(trait = trait_labels))

print(temp_10_plot)
ggsave(filename = paste0(path_out, "/plots/temp_10.png"),
			 plot = temp_10_plot, 
			 bg = "white",
			 width = 150, 
			 height = 250, 
			 units = "mm", 
			 dpi = 1457)

# >>> 25% Quantiles for Temperature <<<

# Set Quantile values
quantiles_25 <- c(0.25, 0.75)
quantile_names_25 <- c("Coldest 25%", "Hottest 25%")

# Calculate partial dependence towards standage for all traits for the quantiles
pdp_temp_25 <- get_partial_dependence(
	best_models = best_models,
	train_data = train_data,
	traits = traits,
	pred.vars = c("standage", "annual_mean_temperature"),
	quantiles = quantiles_25,
	quantile_names = quantile_names_25)

# Plot the results with uncertainty ribbons and facet_grid
temp_25_plot <- ggplot(pdp_temp_25, aes(x = standage, y = yhat, color = group, shape = group)) +
	geom_ribbon(aes(ymin = lower, ymax = upper, fill = group), alpha = 0.2) +
	geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, linewidth = .7) +
	scale_color_manual(values = setNames(c("#0800af", "#c82300"), quantile_names_25)) +
	scale_fill_manual(values = setNames(c("#0800af", "#c82300"), quantile_names_25)) +
	labs(x = "Standage", y = "Predicted Trait Value (log-scaled)", color = "Quantile", shape = "Quantile", fill = "Quantile", title = plot_title_25) + 
	theme_bw() +
	theme(
		legend.position = "top",
		text = element_text(family = "Arial"),
		strip.text.y = element_text(size = 7),
		strip.background = element_rect(fill = "white", color = "black", linewidth = .75)) +
	facet_grid(trait ~ group, scales = "free", labeller = labeller(trait = trait_labels))

print(temp_25_plot)
ggsave(filename = paste0(path_out, "/plots/temp_25.png"),
			 plot = temp_25_plot, 
			 bg = "white",
			 width = 150, 
			 height = 250, 
			 units = "mm", 
			 dpi = 1457)

# >>> 10% Quantiles for Precipitation <<<

# Set Quantile values
quantiles_10 <- c(0.1, 0.9)
quantile_names_10 <- c("Driest 10%", "Wettest 10%")
plot_title <- "Partial Dependence Plot for Annual Precipitation"

# Calculate partial dependence towards standage for all traits for the quantiles
pdp_prcp_10 <- get_partial_dependence(
	best_models = best_models,
	train_data = train_data,
	traits = traits,
	pred.vars = c("standage", "annual_precipitation"),
	quantiles = quantiles_10,
	quantile_names = quantile_names_10)

# Plot the results with uncertainty ribbons and facet_grid
prcp_10_plot <- ggplot(pdp_prcp_10, aes(x = standage, y = yhat, color = group, shape = group)) +
	geom_ribbon(aes(ymin = lower, ymax = upper, fill = group), alpha = 0.2) +
	geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, linewidth = .7) +
	scale_color_manual(values = setNames(c("#0800af", "#c82300"), quantile_names_10)) +
	scale_fill_manual(values = setNames(c("#0800af", "#c82300"), quantile_names_10)) +
	labs(x = "Standage", y = "Predicted Trait Value (log-scaled)", color = "Quantile", shape = "Quantile", fill = "Quantile", title = plot_title) + 
	theme_bw() +
	theme(
		legend.position = "top",
		text = element_text(family = "Arial"),
		strip.text.y = element_text(size = 7),
		strip.background = element_rect(fill = "white", color = "black", linewidth = .75)) +
	facet_grid(trait ~ group, scales = "free", labeller = labeller(trait = trait_labels))

print(prcp_10_plot)
ggsave(filename = paste0(path_out, "/plots/prcp_10.png"),
			 plot = temp_10_plot, 
			 bg = "white",
			 width = 150, 
			 height = 250, 
			 units = "mm", 
			 dpi = 1457)

# >>> 25% Quantiles for Precipitation <<<

# Set Quantile values
quantiles_25 <- c(0.25, 0.75)
quantile_names_25 <- c("Driest 25%", "Wettest 25%")

# Calculate partial dependence towards standage for all traits for the quantiles
pdp_prcp_25 <- get_partial_dependence(
	best_models = best_models,
	train_data = train_data,
	traits = traits,
	pred.vars = c("standage", "annual_precipitation"),
	quantiles = quantiles_25,
	quantile_names = quantile_names_25)

# Plot the results with uncertainty ribbons and facet_grid
prcp_25_plot <- ggplot(pdp_prcp_25, aes(x = standage, y = yhat, color = group, shape = group)) +
	geom_ribbon(aes(ymin = lower, ymax = upper, fill = group), alpha = 0.2) +
	geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, linewidth = .7) +
	scale_color_manual(values = setNames(c("#0800af", "#c82300"), quantile_names_25)) +
	scale_fill_manual(values = setNames(c("#0800af", "#c82300"), quantile_names_25)) +
	labs(x = "Standage", y = "Predicted Trait Value (log-scaled)", color = "Quantile", shape = "Quantile", fill = "Quantile", title = plot_title) + 
	theme_bw() +
	theme(
		legend.position = "top",
		text = element_text(family = "Arial"),
		strip.text.y = element_text(size = 7),
		strip.background = element_rect(fill = "white", color = "black", linewidth = .75)) +
	facet_grid(trait ~ group, scales = "free", labeller = labeller(trait = trait_labels))

print(prcp_25_plot)
ggsave(filename = paste0(path_out, "/plots/prcp_25.png"),
			 plot = temp_25_plot, 
			 bg = "white",
			 width = 150, 
			 height = 250, 
			 units = "mm", 
			 dpi = 1457)

## >>>>> Now look at different biomes <<<<<

# General function to create partial dependence data for each trait based on biomes
get_partial_dependence_biomes <- function(best_models, train_data, traits, pred.var, biome_vars, biome_names = NULL) {
	
	if (is.null(biome_names)) {
		biome_names <- biome_vars
	}
	
	# Print the biome names
	for (biome in biome_vars) {
		cat("Processing:", biome, "\n")
	}
	
	# Separate data for each biome
	data_list <- lapply(biome_vars, function(biome) {
		train_data %>% filter(!!sym(biome) == 1)
	})
	
	pdp_data_combined <- traits %>%
		map_df(function(trait) {
			model <- best_models[[which(traits == trait)]]
			
			# Combine partial dependence data for each biome
			pdp_data_list <- lapply(seq_along(biome_vars), function(i) {
				data_subset <- data_list[[i]]
				
				# Build grid for partial dependence
				pred_grid <- expand.grid(
					standage = unique(data_subset$standage))
				
				# Calculate partial dependence
				pdp_data <- pdp::partial(
					object = model,
					pred.var = pred.var,
					pred.grid = pred_grid,
					train = data_subset,
					plot = FALSE)
				
				# Label the data with the biome name
				pdp_data %>% mutate(group = biome_names[i])
			})
			
			pdp_data_combined <- bind_rows(pdp_data_list) %>%
				mutate(trait = trait,
							 lower = yhat - 1.96 * sd(yhat),
							 upper = yhat + 1.96 * sd(yhat))
			
			return(pdp_data_combined)
		})
	
	return(pdp_data_combined)
}


# Set biome names and variable names; we also consider if a system is dominated by gymnosperms or angiosperms
biomes <- tibble(
	biome_vars = c("biome_boreal_forests_or_taiga", "biome_flooded_grasslands", "biome_mediterranean_woodlands", 
								 "biome_temperate_broadleaf_forests", "biome_temperate_conifer_forests", "biome_temperate_grasslands", 
								 "biome_tundra", "biome_xeric_shrublands"),
	biome_names = c("Boreal Forests or Taiga", "Flooded Grasslands", "Mediterranean Woodlands", 
									"Temperate Broadleaf Forests", "Temperate Conifer Forests", "Temperate Grasslands", 
									"Tundra", "Xeric Shrublands")) %>% 
	mutate(
	dominance = c("Gymnosperm", "Angiosperm", "Angiosperm", 
								"Angiosperm", "Gymnosperm", "Angiosperm", 
								"Gymnosperm", "Gymnosperm"))

# Calculate partial dependence towards standage for all traits for the biomes
pdp_biomes <- get_partial_dependence_biomes(
	best_models = best_models,
	train_data = train_data,
	traits = traits,
	pred.var = "standage",
	biome_vars = biomes$biome_vars,
	biome_names = biomes$biome_names) %>%
	# Add biome dominance
	left_join(
		biomes %>% 
			rename(group = biome_names) %>% 
			dplyr::select(-biome_vars),
		by = "group")

# Plot partial dependence across biomes 
biome_plot <- ggplot(pdp_biomes, aes(x = standage, y = yhat, color = dominance, fill = dominance, group = group)) +
	geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .1, linewidth = .1) +
	geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, linewidth = .7) +
	scale_color_manual(values = dominance_colors) +
	scale_fill_manual(values = dominance_colors) +
	labs(x = "Standage", y = "Predicted Trait Value (log-scaled)", color = "Dominance", fill = "Dominance", title = NULL) + 
	theme_bw() +
	theme(
		legend.position = "top",
		text = element_text(family = "Arial"),
		strip.text.y = element_text(size = 7),
		strip.background = element_rect(fill = "white", color = "black", linewidth = .75)) +
	facet_grid(trait ~ group, scales = "free_y", labeller = labeller(trait = trait_labels))

print(biome_plot)
ggsave(filename = paste0(path_out, "/plots/biomes_plot.png"),
			 plot = biome_plot, 
			 bg = "white",
			 width = 150, 
			 height = 250, 
			 units = "mm", 
			 dpi = 1457)

## Done
stopCluster(cl)