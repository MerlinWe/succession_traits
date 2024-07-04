############################################################################################################################
########################################  MSc Diss. Forest Succession Data Analysis ########################################  
############################################################################################################################

# Objectives:
# 	
# 	1.	Understand the impact of stand age on forest traits.
#   2.	Interpret models results using Shapley values to understand feature importance.
#   3.	Predict how this relationship varies across different climate conditions and ecoregions.

# Outline of the Approach:
# 	
# 	1.	Data Preparation:
# 			•	Ensure the dataset is ready with all necessary variables: traits, stand age, ecoregion, climate variables, and other covariates.
# 			•	Split the data into training and testing sets.
#   2.	Hyperparameter Tuning:
# 			•	Perform grid search combined with cross-validation to find the optimal hyperparameters (mtry, num.trees, min.node.size).
# 	3.	Model Fitting with Optimal Hyperparameters:
# 			•	Use the optimal hyperparameters to fit random forest models for each trait.
# 			•	Ensure stand age is included as a mandatory predictor in all models.
#   4.	Interpretation with Shapley Values:
# 			•	Calculate Shapley values for feature importance using the fastshap package.
# 			•	Visualize the importance of stand age and other covariates across different models.
#   5.	Analysis Across Ecoregions and Climate Conditions:
# 			•	Stratify data by ecoregion and temperature quantiles.
# 			•	Fit separate models for each stratum using the optimal hyperparameters and compare results.

rm(list = ls()) # make sure environment is clean 
set.seed(42)    # set seed for reproducibility
Sys.setenv(R_FORK_SUPPORTED = "false") # disable forking for memory management

suppressWarnings({
	library(doParallel)
	registerDoParallel(num_cores)
})

# Load necessary libraries
library(caret)
library(ranger)
library(pdp)
library(mgcv)
library(rsample)
library(fastshap)
library(doParallel)
library(ggbeeswarm)
library(patchwork)
library(grid)
library(tidyverse)

# Read data and make some minor adjustments
data <- read_csv("/Volumes/ritd-ag-project-rd01pr-dmayn10/merlin/data/fia_traits/sub/traits_rf_clean.csv",
								 col_types = c("d","d","d","d","d","d","d","d","d","d","f","f","d","d","d","d","d","d","d","d","d","d","d","d")) %>%
	rename_with(~ gsub("wmean_", "", .), starts_with("wmean_")) %>%
	rename_with(~ gsub("_015cm", "", .), ends_with("_015cm")) %>%
	rename(min_temperature = min_temperature_of_coldest_month,
				 max_temperature = max_temperature_of_warmest_month)

data <- data %>% 
	# We filter Standage by the upper 10% quantiles
	filter(standage < quantile(standage, 0.9)) %>%
	
	sample_n(n() * 0.1) # use 10% sample for code development 

# Split into training and test sets
split <- initial_split(data, prop = 0.8)

# Extract training and testing datasets
train_data <- training(split)
test_data  <- testing(split)

# Define the traits and covariates
traits <- c("wood_density", "bark_thickness", "conduit_diam", "leaf_n", 
						"specific_leaf_area", "seed_dry_mass", "shade_tolerance", "height")

# Define numeric covariates and all covariates
covariates <- c("standage", "annual_mean_temperature", "annual_precipitation", 
								"temperature_seasonality", "mean_diurnal_range", "min_temperature", 
								"max_temperature", "elevation", "pop_density", 
								"sand_content", "soil_ph", "water_capacity")

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
		summaryFunction = defaultSummary
	)
	
	# Train models
	model <- train(
		formula, 
		data = data,
		method = "ranger",
		trControl = train_control,
		tuneGrid = hyper_grid,
		importance = "permutation"
	)
	
	# Get performance metrics
	results <- model$results
	best_tune <- model$bestTune
	results <- results %>%
		mutate(
			trait = trait, 
			best = rowSums(sapply(names(best_tune), 
														function(param) results[[param]] == best_tune[[param]])) == length(best_tune)
		)
	
	return(list(model = model, results = results))
}

# Register parallel cores before tuning and calculating Shapley values
num_cores <- detectCores() - 2
registerDoParallel(num_cores)

# Perform hyperparameter tuning for each trait and get performance metrics
tuned_models <- traits %>%
	map(~ tune_rf_model(.x, train_data, covariates, hyper_grid))

# Extract the best model for every trait
best_models <- tuned_models %>%
	map(~ .x$model$finalModel)

# Extract performance metrics
global_performance_metrics <- tuned_models %>%
	map_df(~ .x$results) %>%
	mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
	arrange(trait, best)

# Prediction function for fastshap
predict_fn <- function(object, newdata) {
	predict(object, data = newdata)$predictions
}

# Calculate Shapley values for each trait model using parallel processing
# NOTE: The more *nsim* the better, increase this for final run 
shap_values <- foreach(i = seq_along(best_models), .combine = 'rbind', .packages = c('fastshap', 'dplyr')) %dopar% {
	model <- best_models[[i]]
	trait <- traits[i]
	shap_values <- fastshap::explain(model, X = test_data %>% select(all_of(covariates)) %>% as.matrix(),
																	 pred_wrapper = predict_fn, nsim = 10, parallel = TRUE)
	shap_values <- as.data.frame(shap_values)
	shap_values$trait <- trait
	shap_values
}

# Stop the parallel backend
stopImplicitCluster()

# Plot shapley values for every trait
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

print(shap_plot) # examine plot 
ggsave(filename = "/Users/serpent/Documents/MSc/Thesis/Code/analysis/traits/plots/shap_plots.png",
			 plot = shap_plot,
			 bg = "white",
			 width = 250,
			 height = 175,
			 units = "mm",
			 dpi = 1457)

##### Stratify data based on annual mean temperature and make predictions 

# Stratify the data based on the upper and lower 10% quantiles of annual mean temperature
stratify_climate <- function(data, variable) {
	lower_quantile <- quantile(data[[variable]], 0.1)
	upper_quantile <- quantile(data[[variable]], 0.9)
	
	lower_data <- data %>% filter(.[[variable]] <= lower_quantile)
	upper_data <- data %>% filter(.[[variable]] >= upper_quantile)
	
	return(list(lower_data = lower_data, upper_data = upper_data))
}

# Function to train and predict for both stratified datasets and retrieve R-squared values
train_and_predict <- function(data, traits, covariates, hyper_grid) {
	# Split into training and test sets
	set.seed(42)
	split <- initial_split(data, prop = 0.8)
	
	# Extract training and testing datasets
	train_data <- training(split)
	test_data  <- testing(split)
	
	# Perform hyperparameter tuning for each trait and get performance metrics
	tuned_models <- traits %>% map(~ tune_rf_model(.x, train_data, covariates, hyper_grid))
	
	performance_metrics <- tuned_models %>%
		map_df(~ .x$results) %>%
		mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
		arrange(trait, best) %>%
		as_tibble()
	
	# Extract the best model for every trait
	best_models <- tuned_models %>% map(~ .x$model$finalModel)
	
	# Calculate R-squared for each model
	r_squared_values <- best_models %>% map2(traits, ~ caret::R2(predict(.x, test_data)$predictions, test_data[[.y]]))
	
	return(list(best_models = best_models, test_data = test_data, r_squared_values = r_squared_values, quant_performance_metrics = performance_metrics))
}

# Perform stratification and model training
stratified_data_climate <- stratify_climate(data, "annual_mean_temperature")

# Register parallel cores before tuning
registerDoParallel(num_cores)

lower_results <- train_and_predict(stratified_data_climate$lower_data, traits, covariates, hyper_grid)
upper_results <- train_and_predict(stratified_data_climate$upper_data, traits, covariates, hyper_grid)

# Stop the parallel backend
stopImplicitCluster()

# Retrieve performance metrics for lower and upper quantiles
quant_performance_metrics <- bind_rows(
	lower_results$quant_performance_metrics %>% mutate(group = "Lower 10%"),
	upper_results$quant_performance_metrics %>% mutate(group = "Upper 10%"))

# Function to create individual partial dependence plots with R-squared values and custom colours

trait_title <- function(trait) { # helper function for plot titles
	case_when(
		trait == "wood_density" ~ "Wood Density",
		trait == "bark_thickness" ~ "Bark Thickness",
		trait == "conduit_diam" ~ "Conduit Diameter",
		trait == "leaf_n" ~ "Leaf Nitrogen",
		trait == "specific_leaf_area" ~ "Specific Leaf Area",
		trait == "seed_dry_mass" ~ "Seed Dry Mass",
		trait == "shade_tolerance" ~ "Shade Tolerance",
		trait == "height" ~ "Tree Height",
		TRUE ~ NA_character_
	)
}

plot_partial_dependence_climate <- function(lower_results, upper_results, trait) {
	
	# Extract models
	lower_model <- lower_results$best_models[[which(traits == trait)]]
	upper_model <- upper_results$best_models[[which(traits == trait)]]
	
	# Generate partial dependence data
	pdp_lower <- pdp::partial(lower_model, pred.var = "standage", train = lower_results$test_data)
	pdp_upper <- pdp::partial(upper_model, pred.var = "standage", train = upper_results$test_data)
	
	# Add quantile identifiers
	pdp_lower$group <- "Lower 10%"
	pdp_upper$group <- "Upper 10%"
	
	# Combine partial dependence data
	pdp_data <- bind_rows(pdp_lower, pdp_upper)
	
	# Extract R-squared values
	lower_r2 <- lower_results$r_squared_values[[which(traits == trait)]]
	upper_r2 <- upper_results$r_squared_values[[which(traits == trait)]]
	
	# Raw data from lower and upper quantiles to show how deterministic it is?
	# lower_raw <- lower_results$test_data %>% mutate(group = "Lower 10%")
	# upper_raw <- upper_results$test_data %>% mutate(group = "Upper 10%")
	# raw_data <- bind_rows(lower_raw, upper_raw)
	
	ggplot(pdp_data, aes(x = standage, y = yhat, color = group, shape = group)) +
		
		# Smoothed curve using spline smoothing
		geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, linewidth = .7, aes(fill = group)) + 
		# Raw partial dependencies 
		geom_line(linewidth = .4) +
		
		scale_color_manual(values = c("Lower 10%" = "#0800af", "Upper 10%" = "#c82300")) +
		scale_shape_manual(values = c(16, 18)) + 
		labs(title = trait_title(trait),
				 subtitle = paste("Lower 10% R\u00B2:", round(lower_r2, 2), "| Upper 10% R\u00B2:", round(upper_r2, 2)),
				 x = "Standage", y = "Trait value (log scaled)") + 
		theme_bw() +
		theme(
			legend.position = "none",
			legend.title = element_blank(),
			legend.text = element_text(size = 10),
			text = element_text(family = "Palatino"),
			plot.title = element_text(face = "bold"),
			plot.subtitle = element_text(hjust = 0.5))
}

# Generate plots for all traits
pdp_plots_climate <- list()
for (trait in traits) {
	plot <- plot_partial_dependence_climate(lower_results, upper_results, trait)
	pdp_plots_climate[[trait]] <- plot
}

# Combine the plots using patchwork
climate_plot <- ((pdp_plots_climate$wood_density + theme(axis.title.x = element_text(colour = "transparent"))) +
										(pdp_plots_climate$bark_thickness+  theme(axis.title = element_text(colour = "transparent")))) /
	
	((pdp_plots_climate$conduit_diam + theme(axis.title.x = element_text(colour = "transparent"))) +
	 	(pdp_plots_climate$leaf_n + theme(axis.title = element_text(colour = "transparent")))) /
	
	((pdp_plots_climate$specific_leaf_area + theme(axis.title.x = element_text(colour = "transparent"))) +
	 	(pdp_plots_climate$seed_dry_mass + theme(axis.title = element_text(colour = "transparent")))) /
	
	(pdp_plots_climate$shade_tolerance + 
	 	(pdp_plots_climate$height + theme(axis.title.y = element_text(colour = "transparent")))) +
	
	plot_annotation(
		theme = theme(
			text = element_text(family = "Palatino"),
			plot.title = element_text(face = "bold")))

print(climate_plot) # examine figure
ggsave(filename = "/Users/serpent/Documents/MSc/Thesis/Code/analysis/traits/plots/clim_plots.png",
			 plot = climate_plot,
			 bg = "white",
			 width = 200,
			 height = 250,
			 units = "mm",
			 dpi = 1457)

##### Stratify data based on management and make predictions 

# Note: managed == 0 data will be down sampled to avoid bias 

# Stratify the data based on the managed variable
stratify_managed <- function(data, variable) {
	managed_data_0 <- data %>% filter(.[[variable]] == 0)
	managed_data_1 <- data %>% filter(.[[variable]] == 1)
	
	# Down sample the majority class (managed == 0) to match the minority class (managed == 1)
	set.seed(42)
	downsampled_data_0 <- managed_data_0 %>% sample_n(nrow(managed_data_1))
	
	return(list(managed_data_0 = downsampled_data_0, managed_data_1 = managed_data_1))
}

# Perform stratification and model training for managed variable
stratified_data_managed <- stratify_managed(data, "managed")

# Register parallel cores before tuning
registerDoParallel(num_cores)

managed_results_0 <- train_and_predict(stratified_data_managed$managed_data_0, traits, covariates, hyper_grid)
managed_results_1 <- train_and_predict(stratified_data_managed$managed_data_1, traits, covariates, hyper_grid)

# Stop the parallel backend
stopImplicitCluster()

# Retrieve performance metrics for managed levels
managed_performance_metrics <- bind_rows(
	managed_results_0$quant_performance_metrics %>% mutate(group = "Managed 0"),
	managed_results_1$quant_performance_metrics %>% mutate(group = "Managed 1")
)

# Function to create partial dependence plots for managed variable with R-squared values and custom colors
plot_partial_dependence_managed <- function(managed_results_0, managed_results_1, trait) {
	managed_model_0 <- managed_results_0$best_models[[which(traits == trait)]]
	managed_model_1 <- managed_results_1$best_models[[which(traits == trait)]]
	
	pdp_managed_0 <- pdp::partial(managed_model_0, pred.var = "standage", train = managed_results_0$test_data)
	pdp_managed_1 <- pdp::partial(managed_model_1, pred.var = "standage", train = managed_results_1$test_data)
	
	pdp_managed_0$group <- "Managed 0"
	pdp_managed_1$group <- "Managed 1"
	
	pdp_data <- bind_rows(pdp_managed_0, pdp_managed_1)
	
	managed_r2_0 <- managed_results_0$r_squared_values[[which(traits == trait)]]
	managed_r2_1 <- managed_results_1$r_squared_values[[which(traits == trait)]]
	
	ggplot(pdp_data, aes(x = standage, y = yhat, color = group, shape = group)) +
		geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, linewidth = .7, aes(fill = group)) + 
		geom_line(linewidth = .4) +
		scale_color_manual(values = c("Managed 0" = "#0800af", "Managed 1" = "#c82300")) +
		scale_shape_manual(values = c(16, 18)) + 
		labs(title = trait_title(trait),
				 subtitle = paste("Managed 0 R\u00B2:", round(managed_r2_0, 2), "| Managed 1 R\u00B2:", round(managed_r2_1, 2)),
				 x = "Standage", y = "Trait value (log scaled)") + 
		theme_bw() +
		theme(
			legend.position = "none",
			legend.title = element_blank(),
			legend.text = element_text(size = 10),
			text = element_text(family = "Palatino"),
			plot.title = element_text(face = "bold"),
			plot.subtitle = element_text(hjust = 0.5)
		)
}

# Generate plots for all traits
pdp_plots_managed <- list()
for (trait in traits) {
	plot <- plot_partial_dependence_managed(managed_results_0, managed_results_1, trait)
	pdp_plots_managed[[trait]] <- plot
}

# Combine the plots using patchwork
managed_plot <- ((pdp_plots_managed$wood_density + theme(axis.title.x = element_text(colour = "transparent"))) +
								 	(pdp_plots_managed$bark_thickness + theme(axis.title = element_text(colour = "transparent")))) /
	((pdp_plots_managed$conduit_diam + theme(axis.title.x = element_text(colour = "transparent"))) +
	 	(pdp_plots_managed$leaf_n + theme(axis.title = element_text(colour = "transparent")))) /
	((pdp_plots_managed$specific_leaf_area + theme(axis.title.x = element_text(colour = "transparent"))) +
	 	(pdp_plots_managed$seed_dry_mass + theme(axis.title = element_text(colour = "transparent")))) /
	(pdp_plots_managed$shade_tolerance +
	 	(pdp_plots_managed$height + theme(axis.title.y = element_text(colour = "transparent")))) +
	plot_annotation(
		theme = theme(
			text = element_text(family = "Palatino"),
			plot.title = element_text(face = "bold")))

print(managed_plot) # examine figure

# Save the combined plot
ggsave(filename = "/Users/serpent/Documents/MSc/Thesis/Code/analysis/traits/plots/managed_plots.png",
			 plot = managed_plot,
			 bg = "white",
			 width = 200,
			 height = 250,
			 units = "mm",
			 dpi = 1457)
