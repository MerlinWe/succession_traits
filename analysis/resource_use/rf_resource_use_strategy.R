rm(list = ls()) # make sure environment is clean 
set.seed(42)    # set seed for reproducibility

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
library(tidyverse)

# Read data 
data <- read_csv("/Volumes/ritd-ag-project-rd01pr-dmayn10/merlin/data/fia_traits/sub/rus_rf_clean.csv",
								 col_types = c("d","d","f","f","d","d","d","d","d","d","d","d","d","d","d","d")) %>%
	rename_with(~ gsub("wmean_", "", .), starts_with("wmean_")) %>%
	rename_with(~ gsub("_015cm", "", .), ends_with("_015cm")) %>%
	rename(min_temperature = min_temperature_of_coldest_month,
				 max_temperature = max_temperature_of_warmest_month)

data <- data %>% 
	# We filter Standage by the upper 10% quantiles
	filter(standage < quantile(standage, 0.9)) 

# Split into training and test sets
split <- initial_split(data, prop = 0.8)

# Extract training and testing datasets
train_data <- training(split)
test_data  <- testing(split)

# Define the dependent variable and covariates
dependent_variable <- "resource_use_score"

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
tune_rf_model <- function(dependent_variable, data, covariates, hyper_grid) {
	
	# Define model formula
	formula <- as.formula(paste(dependent_variable, "~", paste(covariates, collapse = " + ")))
	
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
			best = rowSums(sapply(names(best_tune), 
														function(param) results[[param]] == best_tune[[param]])) == length(best_tune)
		)
	
	return(list(model = model, results = results))
}

# Register parallel cores before tuning and calculating Shapley values
num_cores <- detectCores() - 2
registerDoParallel(num_cores)

# Perform hyperparameter tuning for the dependent variable and get performance metrics
tuned_model <- tune_rf_model(dependent_variable, train_data, covariates, hyper_grid)

# Extract the best model
best_model <- tuned_model$model$finalModel

# Extract performance metrics
global_performance_metrics <- tuned_model$results %>%
	mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
	arrange(best)

# Prediction function for fastshap
predict_fn <- function(object, newdata) {
	predict(object, data = newdata)$predictions
}

# Calculate Shapley values for the model using parallel processing
shap_values <- fastshap::explain(best_model, X = test_data %>% select(all_of(covariates)) %>% as.matrix(),
																 pred_wrapper = predict_fn, nsim = 100, parallel = TRUE)

# Convert Shapley values to a dataframe
shap_values <- as.data.frame(shap_values)

# Plot shapley values
shap_long <- shap_values %>%
	pivot_longer(cols = standage:water_capacity, names_to = "feature", values_to = "shap") %>%
	as_tibble() 

shap_plot <- shap_long %>%
	ggplot(aes(x = feature, y = shap, color = shap)) +
	geom_quasirandom(alpha = 0.5) +
	scale_color_viridis_c(option = "viridis") +
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
ggsave(filename = "/Users/serpent/Documents/MSc/Thesis/Code/analysis/resource_use/plots/shap_plots.png",
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
train_and_predict <- function(data, dependent_variable, covariates, hyper_grid) {
	# Split into training and test sets
	set.seed(42)
	split <- initial_split(data, prop = 0.8)
	
	# Extract training and testing datasets
	train_data <- training(split)
	test_data  <- testing(split)
	
	# Perform hyperparameter tuning for the dependent variable and get performance metrics
	tuned_model <- tune_rf_model(dependent_variable, train_data, covariates, hyper_grid)
	
	performance_metrics <- tuned_model$results %>%
		mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
		arrange(best) %>%
		as_tibble()
	
	# Extract the best model
	best_model <- tuned_model$model$finalModel
	
	# Calculate R-squared for the model
	r_squared_value <- caret::R2(predict(best_model, test_data)$predictions, test_data[[dependent_variable]])
	
	return(list(best_model = best_model, test_data = test_data, r_squared_value = r_squared_value, quant_performance_metrics = performance_metrics))
}

# Perform stratification and model training
stratified_data_climate <- stratify_climate(data, "annual_mean_temperature")

# Register parallel cores before tuning
registerDoParallel(num_cores)

lower_results <- train_and_predict(stratified_data_climate$lower_data, dependent_variable, covariates, hyper_grid)
upper_results <- train_and_predict(stratified_data_climate$upper_data, dependent_variable, covariates, hyper_grid)

# Stop the parallel backend
stopImplicitCluster()

# Retrieve performance metrics for lower and upper quantiles
quant_performance_metrics <- bind_rows(
	lower_results$quant_performance_metrics %>% mutate(group = "Lower 10%"),
	upper_results$quant_performance_metrics %>% mutate(group = "Upper 10%"))

# Create partial dependence plot for stratified data with R-squared values
plot_partial_dependence_climate <- function(lower_results, upper_results) {
	lower_model <- lower_results$best_model
	upper_model <- upper_results$best_model
	
	pdp_lower <- pdp::partial(lower_model, pred.var = "standage", train = lower_results$test_data)
	pdp_upper <- pdp::partial(upper_model, pred.var = "standage", train = upper_results$test_data)
	
	pdp_lower$group <- "Lower 10%"
	pdp_upper$group <- "Upper 10%"
	
	pdp_data <- bind_rows(pdp_lower, pdp_upper)
	
	lower_r2 <- lower_results$r_squared_value
	upper_r2 <- upper_results$r_squared_value
	
	ggplot(pdp_data, aes(x = standage, y = yhat, color = group, shape = group)) +
		geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, linewidth = .7, aes(fill = group)) + 
		geom_line() +
		scale_color_manual(values = c("Lower 10%" = "#0800af", "Upper 10%" = "#c82300")) +
		scale_shape_manual(values = c(16, 18)) + 
		labs(title = NULL,
				 subtitle = paste("Lower 10% R\u00B2:", round(lower_r2, 2), "| Upper 10% R\u00B2:", round(upper_r2, 2)),
				 x = "Standage", y = "Predicted Resource Use Score") + 
		theme_bw() +
		theme(
			legend.position = "bottom",
			legend.title = element_blank(),
			legend.text = element_text(size = 10),
			text = element_text(family = "Palatino"),
			plot.title = element_text(face = "bold", hjust = 0.5),
			plot.subtitle = element_text(hjust = 0.5)
		)
}

# Generate plot for partial dependence
climate_plot <- plot_partial_dependence_climate(lower_results, upper_results)
print(climate_plot) # examine figure


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

managed_results_0 <- train_and_predict(stratified_data_managed$managed_data_0, dependent_variable, covariates, hyper_grid)
managed_results_1 <- train_and_predict(stratified_data_managed$managed_data_1, dependent_variable, covariates, hyper_grid)

# Stop the parallel backend
stopImplicitCluster()

# Retrieve performance metrics for managed levels
managed_performance_metrics <- bind_rows(
	managed_results_0$quant_performance_metrics %>% mutate(group = "Managed 0"),
	managed_results_1$quant_performance_metrics %>% mutate(group = "Managed 1")
)

# Create partial dependence plot for managed variable with R-squared values
plot_partial_dependence_managed <- function(managed_results_0, managed_results_1) {
	managed_model_0 <- managed_results_0$best_model
	managed_model_1 <- managed_results_1$best_model
	
	pdp_managed_0 <- pdp::partial(managed_model_0, pred.var = "standage", train = managed_results_0$test_data)
	pdp_managed_1 <- pdp::partial(managed_model_1, pred.var = "standage", train = managed_results_1$test_data)
	
	pdp_managed_0$group <- "Managed 0"
	pdp_managed_1$group <- "Managed 1"
	
	pdp_data <- bind_rows(pdp_managed_0, pdp_managed_1)
	
	managed_r2_0 <- managed_results_0$r_squared_value
	managed_r2_1 <- managed_results_1$r_squared_value
	
	ggplot(pdp_data, aes(x = standage, y = yhat, color = group, shape = group)) +
		geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, linewidth = .7, aes(fill = group)) + 
		geom_line(linewidth = .4) +
		scale_color_manual(values = c("Managed 0" = "#0800af", "Managed 1" = "#c82300")) +
		scale_shape_manual(values = c(16, 18)) + 
		labs(title = NULL,
				 subtitle = paste("Managed 0 R\u00B2:", round(managed_r2_0, 2), "| Managed 1 R\u00B2:", round(managed_r2_1, 2)),
				 x = "Standage", y = NULL) + 
		theme_bw() +
		theme(
			legend.position = "bottom",
			legend.title = element_blank(),
			legend.text = element_text(size = 10),
			text = element_text(family = "Palatino"),
			plot.title = element_text(face = "bold", hjust = 0.5),
			plot.subtitle = element_text(hjust = 0.5)
		)
}

# Generate plot for partial dependence
managed_plot <- plot_partial_dependence_managed(managed_results_0, managed_results_1)
print(managed_plot) # examine figure

# Build Compound Figure
rus_pred <- climate_plot + managed_plot + 
	plot_annotation(
		title = "Partial Dependence of Standage on Resource Use Score",
		theme = theme(
			plot.title = element_text(face = "bold", hjust = 0.5, family = "Palatino")))

print(rus_pred)
ggsave(filename = "/Users/serpent/Documents/MSc/Thesis/Code/analysis/resource_use/plots/rus_pred_plot.png",
			 plot = rus_pred,
			 bg = "white",
			 width = 250,
			 height = 120,
			 units = "mm",
			 dpi = 1457)
