############################################################################################################################
########################################  MSc Diss. Forest Succession Data Analysis ########################################  
############################################################################################################################

# Objectives:
# 	
# 	1.	Understand the impact of stand age on forest traits.
#   2.	Examine how this relationship varies across different ecoregions and climate conditions.
#   3.	Interpret model results using Shapley values to understand feature importance.
# 
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

rm(list = ls()) # clean environment 

# Load necessary libraries
library(DT)
library(caret)
library(ranger)
library(pdp)
library(rsample)
library(fastshap)
library(doParallel)
library(ggbeeswarm)
library(tidyverse)

# Read data and make some minor adjustments
data <- read_csv("/Volumes/ritd-ag-project-rd01pr-dmayn10/merlin/data/fia_traits/traits_rf_clean.csv",
								 col_types = c("d","d","d","d","d","d","d","d","d","d","f","f","d","d","d","d","d","d","d","d","d","d","d","d")) %>%
	rename_with(~ gsub("wmean_", "", .), starts_with("wmean_")) %>%
	rename_with(~ gsub("_015cm", "", .), ends_with("_015cm")) %>%
	rename(min_temperature = min_temperature_of_coldest_month,
				 max_temperature = max_temperature_of_warmest_month)

data <- data %>% 
	filter(standage < 150) %>% # filter standage? 
	sample_n(1000) # only for code development 

# Split into training and test sets
set.seed(42)
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

# Perform cross-validation, find the best hyperparameters, capture performance metrics
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
performance_metrics <- tuned_models %>%
	map_df(~ .x$results) %>%
	mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
	arrange(trait, best)

# Prediction function for fastshap
predict_fn <- function(object, newdata) {
	predict(object, data = newdata)$predictions
}

# Calculate Shapley values for each trait model using parallel processing
shap_values <- foreach(i = seq_along(best_models), .combine = 'rbind', .packages = c('fastshap', 'dplyr')) %dopar% {
	model <- best_models[[i]]
	trait <- traits[i]
	shap_values <- fastshap::explain(model, X = test_data %>% select(all_of(covariates)) %>% as.matrix(),
																	 pred_wrapper = predict_fn, nsim = 100, parallel = TRUE)
	shap_values <- as.data.frame(shap_values)
	shap_values$trait <- trait
	shap_values
}

# Stop the parallel backend
stopImplicitCluster()

# Reshape the combined data frame to long format
shap_values_long <- shap_values %>%
	pivot_longer(cols = -c(trait), names_to = "feature", values_to = "shap") %>%
	as_tibble()

shap_values_long %>%
	ggplot( aes(x = feature, y = shap, color = shap)) +
	geom_quasirandom(alpha = 0.5) +
	facet_wrap(~trait, ncol = 4, nrow = 2) +
	scale_color_viridis_c(option = "viridis") +
	coord_flip() +
	labs(x = "Feature", y = "Shapley Value") +
	theme_bw() +
	theme(text = element_text(family = "Palatino"),
				legend.position = "none")
