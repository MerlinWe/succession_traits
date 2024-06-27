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

library(caret)
library(ranger)
library(pdp)
library(fastshap)
library(doParallel)
library(ggbeeswarm)
library(tidyverse)

# Clear the workspace
rm(list = ls())

# Load the data
traits <- read_csv("/Volumes/ritd-ag-project-rd01pr-dmayn10/merlin/data/fia_traits/traits_rf_clean.csv")

# Split into training and test sets
set.seed(42)
data <- traits %>% 
	filter(standage < 150) %>%
	sample_n(1000) # only for code development 

train_indices <- sample(1:nrow(data), size = 0.8 * nrow(data))
train_data <- data[train_indices, ]
test_data <- data[-train_indices, ]

# Define the traits and covariates
traits <- c("wmean_wood_density", "wmean_bark_thickness", "wmean_conduit_diam", "wmean_leaf_n", 
						"wmean_specific_leaf_area", "wmean_seed_dry_mass", "wmean_shade_tolerance", "wmean_height")

covariates <- c("standage", "annual_mean_temperature", "annual_precipitation", 
								"temperature_seasonality", "mean_diurnal_range", 
								"min_temperature_of_coldest_month", "max_temperature_of_warmest_month", 
								"elevation", "pop_density", "sand_content_015cm", 
								"soil_ph_015cm", "water_capacity_015cm")

# Define a grid of hyperparameters including mtry, min.niode.size, and split rule
hyper_grid <- expand.grid(
	mtry = seq(2, length(covariates), by = 2),
	splitrule = "variance",
	min.node.size = c(1, 5, 10)
)

# Function to perform cross-validation, find the best hyperparameters, and capture performance metrics
tune_rf_model <- function(trait, data, covariates, hyper_grid) {
	formula <- as.formula(paste(trait, "~", paste(covariates, collapse = " + ")))
	
	train_control <- trainControl(method = "cv", number = 5, savePredictions = "final", summaryFunction = defaultSummary)
	
	model <- train(
		formula, data = data,
		method = "ranger",
		trControl = train_control,
		tuneGrid = hyper_grid,
		importance = "permutation"
	)
	
	# Extract performance metrics
	results <- model$results
	best_tune <- model$bestTune
	results <- results %>%
		mutate(trait = trait, 
					 best = rowSums(sapply(names(best_tune), function(param) results[[param]] == best_tune[[param]])) == length(best_tune))
	
	return(list(model = model, results = results))
}

# Perform hyperparameter tuning for each trait and capture performance metrics
tuned_models <- lapply(traits, function(trait) tune_rf_model(trait, train_data, covariates, hyper_grid))

# Extract the best models and performance metrics
best_models <- lapply(tuned_models, function(x) x$model$finalModel)
performance_metrics <- do.call(rbind, lapply(tuned_models, function(x) x$results))

# •	RMSE: Root Mean Squared Error. This is a measure of the average deviation of the predicted values from the actual values. It is calculated as the square root of the average of squared differences between predicted and actual values. Lower values indicate better model performance.
# •	Rsquared: R-squared. This is a measure of how well the model’s predictions approximate the actual data. It represents the proportion of variance in the dependent variable that is predictable from the independent variables. Values range from 0 to 1, with higher values indicating better model performance.
# •	MAE: Mean Absolute Error. This is the average of the absolute differences between predicted and actual values. It provides a measure of how far the predictions are from the actual values on average. Lower values indicate better model performance.
# •	RMSESD: Standard Deviation of Root Mean Squared Error. This measures the variability of the RMSE across different cross-validation folds. Lower values indicate more consistent model performance.
# •	RsquaredSD: Standard Deviation of R-squared. This measures the variability of the R-squared values across different cross-validation folds. Lower values indicate more consistent model performance.
# •	MAESD: Standard Deviation of Mean Absolute Error. This measures the variability of the MAE across different cross-validation folds. Lower values indicate more consistent model performance.

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
	ggplot(shap_long, aes(x = feature, y = shap, color = shap)) +
		geom_quasirandom(alpha = 0.5) +
		scale_color_viridis_c(option = "viridis") +
		coord_flip() +
		labs(title = paste(trait), x = "Feature", y = "Shapley Value") +
		theme_bw() +
		theme(text = element_text(family = "Palatino"),
					legend.position = "none")
}

plots <- lapply(seq_along(traits), function(i) plot_shap_values(shap_values[[i]], covariates, traits[i]))
combined_plot <- cowplot::plot_grid(plotlist = plots, ncol = 4, nrow = 2)
print(combined_plot)


## Plot predictions using partial residual curves
plot_partial_dependence <- function(model, data, trait, covariate) {
	
	# Generate partial dependence data
	pd_data <- partial(model, pred.var = covariate, train = data, grid.resolution = 100)
	
	# Convert to data frame for ggplot
	pd_df <- as.data.frame(pd_data)
	
	# Plot partial dependence data
	ggplot(pd_df, aes(x = !!sym(covariate), y = yhat)) +
		geom_line(color = "blue") +
		labs(title = paste("Partial Dependence of", trait, "on", covariate),
				 x = covariate, y = "Predicted Value") +
		theme_bw() +
		theme(text = element_text(family = "Palatino"))
}

# Generate and plot partial dependence for each trait
partial_plots <- lapply(seq_along(traits), function(i) {
	plot_partial_dependence(best_models[[i]], test_data, traits[i], "standage")
})

# Combine all plots into a single plot
combined_partial_plot <- cowplot::plot_grid(plotlist = partial_plots, ncol = 4, nrow = 2)
print(combined_partial_plot)
