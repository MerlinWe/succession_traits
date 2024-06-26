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


library(ggbeeswarm)
library(ranger)
library(fastshap)
library(doParallel)
library(caret)
library(tidyverse)

# Clear the workspace
rm(list = ls())

# Load the data
traits <- read_csv("/Volumes/ritd-ag-project-rd01pr-dmayn10/merlin/data/fia_traits/traits_rf_clean.csv")

# Split into training and test sets
set.seed(42)
data <- traits %>% sample_n(1000)

train_indices <- sample(1:nrow(data), size = 0.8 * nrow(data))
train_data <- data[train_indices, ]
test_data <- data[-train_indices, ]

# Define the traits and covariates
traits <- c("wmean_wood_density", "wmean_bark_thickness", "wmean_conduit_diam", "wmean_leaf_n", 
						"wmean_specific_leaf_area", "wmean_seed_dry_mass", "wmean_shade_tolerance", "wmean_dia", "wmean_height")

covariates <- c("standage", "annual_mean_temperature", "annual_precipitation", 
								"temperature_seasonality", "mean_diurnal_range", 
								"min_temperature_of_coldest_month", "max_temperature_of_warmest_month", 
								"elevation", "pop_density", "sand_content_015cm", 
								"soil_ph_015cm", "water_capacity_015cm")

# Define a grid of hyperparameters including splitrule
hyper_grid <- expand.grid(
	mtry = seq(2, length(covariates), by = 2),
	splitrule = "variance",
	min.node.size = c(1, 5, 10)
)

# Function to perform cross-validation and find the best hyperparameters
tune_rf_model <- function(trait, data, covariates, hyper_grid) {
	formula <- as.formula(paste(trait, "~", paste(covariates, collapse = " + ")))
	
	train_control <- trainControl(method = "cv", number = 5)
	
	train(
		formula, data = data,
		method = "ranger",
		trControl = train_control,
		tuneGrid = hyper_grid,
		importance = "permutation"
	)
}

# Perform hyperparameter tuning for each trait
tuned_models <- lapply(traits, function(trait) tune_rf_model(trait, train_data, covariates, hyper_grid))

# Extract the best models
best_models <- lapply(tuned_models, function(model) model$finalModel)

# Define prediction function for fastshap
predict_fn <- function(object, newdata) {
	predict(object, data = newdata)$predictions
}

# Calculate Shapley values
registerDoParallel(detectCores() - 1)
shap_values <- lapply(best_models, function(model) {
	fastshap::explain(model, X = test_data[, covariates, drop = FALSE], pred_wrapper = predict_fn, nsim = 100, parallel = TRUE)
})


# Define a function to plot Shapley values
plot_shap_values <- function(shap_values, covariates, trait) {
	# Reshape Shapley values to a long format
	shap_df <- as.data.frame(shap_values)
	shap_df$id <- 1:nrow(shap_df)
	shap_long <- tidyr::gather(shap_df, key = "feature", value = "shap", -id)
	
	# Ensure the features match the covariates
	shap_long$feature <- rep(covariates, each = nrow(shap_df))
	
	# Plotting
	ggplot(shap_long, aes(x = feature, y = shap)) +
		geom_quasirandom(alpha = 0.5) +
		coord_flip() +
		labs(title = paste(trait), x = "Feature", y = "Shapley Value") +
		theme_minimal()
}

# Plot Shapley values for each trait
plots <- lapply(seq_along(traits), function(i) plot_shap_values(shap_values[[i]], covariates, traits[i]))

# Arrange plots
gridExtra::grid.arrange(grobs = plots, ncol = 4, nrow = 3)

