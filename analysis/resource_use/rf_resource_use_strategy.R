############################################################################################################################
########################################  MSc Diss. Forest Succession Data Analysis ########################################  
############################################################################################################################

# Objectives:
# 	1.	Understand the impact of stand age (and other features) on resource use strategy using Shapley values.
#   2.	Predict how patterns across successional time vary across different climate conditions, management, and ecoregions.

# Outline of the Approach:
# 	-	Data Preparation
#   -	Hyperparameter Tuning using grid search combined with cross-validation
# 	-	Model Fitting with Optimal Hyperparameters
#   -	Interpretation with Shapley Values
#   -	Analysis Across Ecoregions, Management and Climate Conditions using strata based on quantiles and factors

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

# Set path to run either on local device 
path_in <- "/Volumes/ritd-ag-project-rd01pr-dmayn10/merlin/data/fia_traits/sub" 
path_out <- "/Users/serpent/Documents/MSc/Thesis/Code/analysis/resource_use" 

# Set paths to run on threadripper
path_in <- "/home/RSD_storage/merlin/data/fia_traits/sub"
path_out <- "/home/RSD_storage/merlin/Code/analysis/resource_use" 

# Read data 
data <- read_csv(paste0(path_in, "/rus_rf_clean.csv"),
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
# If on threadripper set to 32, if local, set to 10 !!
num_cores <- 10
cl <- makeCluster(num_cores)
registerDoParallel(cl)
# Perform hyperparameter tuning for the dependent variable and get performance metrics
tuned_model <- tune_rf_model(dependent_variable, train_data, covariates, hyper_grid)

# Extract the best model
best_model <- tuned_model$model$finalModel

# Extract performance metrics
global_performance_metrics <- tuned_model$results %>%
	mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
	arrange(best) %>%
	write_csv(file = paste0(path_out, "/tables/global_rf_performance_metrics.csv"))

# Prediction function for fastshap
predict_fn <- function(object, newdata) {
	predict(object, data = newdata)$predictions
}

# Calculate Shapley values for the model using parallel processing
shap_values <- fastshap::explain(best_model, X = test_data %>% select(all_of(covariates)) %>% as.matrix(),
																 pred_wrapper = predict_fn, nsim = 100, parallel = TRUE)


## Visualize shapley values... first get a tibble for plotting
shap_long <- shap_values %>%
	as_tibble() %>%
	pivot_longer(cols = everything(), names_to = "feature", values_to = "shap")
	
# Sum absolute Shapley values to determine overall importance
feature_importance <- shap_long %>%
	group_by(feature) %>%
	summarize(importance = sum(abs(shap)), .groups = "drop") %>%
	arrange(importance)

shap_plot <- shap_long %>%
	mutate(feature = factor(feature, levels = feature_importance$feature)) %>%
	ggplot(aes(x = feature, y = shap, color = shap)) +
	geom_quasirandom(alpha = 0.5) +
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
ggsave(filename = paste0(path_out, "/plots/global_shap_plots.png"),
			 plot = shap_plot,
			 bg = "white",
			 width = 280,
			 height = 200,
			 units = "mm",
			 dpi = 1457)

## Plot feature importance 
importance_plot <- ggplot(feature_importance, aes(x = reorder(feature, importance), y = importance, fill = feature)) +
	geom_bar(stat = "identity", colour = "black", alpha = .7) +
	coord_flip() +
	scale_fill_viridis_d()+
	labs(title = "Feature Importance for Explaining Resource Use Strategy",
			 x = NULL,
			 y = "Overall Importance",
			 color = "Average Importance") +
	theme_bw(base_size = 15) +
	theme(legend.position = "none",
				text = element_text(family = "Arial"),
				plot.title = element_text(hjust = 0.5),
				plot.subtitle = element_text(hjust = 0.5))

print(importance_plot) # examine plot 
ggsave(filename = paste0(path_out, "/plots/global_importance_plot.png"),
			 plot = importance_plot,
			 bg = "white",
			 width = 250,
			 height = 200,
			 units = "mm",
			 dpi = 1457)

## Plot relationship of traits with standage 

feature_data <- shap_long %>% 
	filter(feature == "standage")
observed_values <- test_data %>% pull(standage)
	
shap_standage_plot <- ggplot(feature_data, aes(x = observed_values, y = shap, color = shap)) +
		geom_beeswarm(alpha = 0.5, colour = "darkgreen") +
		geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.5) +
		labs(x = "standage", y = "Influence (Shapley)") +
		theme_bw() +
		theme(text = element_text(family = "Arial"),
					legend.position = "right",
					legend.key.width = unit(0.5, "cm"),
					legend.key.height = unit(2, "cm"),
					legend.box.background = element_rect(color = "black", linewidth = .75),
					plot.title = element_text(hjust = 0.5)) 

print(shap_standage_plot)
ggsave(filename = paste0(path_out, "/plots/global_shap_standage_plots.png"),
			 plot = shap_standage_plot,
			 bg = "white",
			 width = 250,
			 height = 130,
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

lower_results <- train_and_predict(stratified_data_climate$lower_data, dependent_variable, covariates, hyper_grid)
upper_results <- train_and_predict(stratified_data_climate$upper_data, dependent_variable, covariates, hyper_grid)

# Retrieve performance metrics for lower and upper quantiles
quant_performance_metrics <- bind_rows(
	lower_results$quant_performance_metrics %>% mutate(group = "Lower 10%"),
	upper_results$quant_performance_metrics %>% mutate(group = "Upper 10%")) %>%
	
	write_csv(file = paste0(path_out, "/tables/climate_quantiles_rf_performance_metrics.csv"))

# Create partial dependence plot for stratified data with R-squared values
plot_partial_dependence_climate <- function(lower_results, upper_results, dependent_variable) {
	
	lower_model <- lower_results$best_model
	upper_model <- upper_results$best_model
	
	pdp_lower <- pdp::partial(lower_model, pred.var = "standage", train = lower_results$test_data)
	pdp_upper <- pdp::partial(upper_model, pred.var = "standage", train = upper_results$test_data)
	
	pdp_lower$group <- "Lower 10%"
	pdp_upper$group <- "Upper 10%"
	
	pdp_data <- bind_rows(pdp_lower, pdp_upper)
	
	lower_r2 <- lower_results$r_squared_value
	upper_r2 <- upper_results$r_squared_value
	
	lower_preds <- predict(lower_model, lower_results$test_data)
	upper_preds <- predict(upper_model, upper_results$test_data)
	
	lower_results$test_data$yhat <- lower_preds$predictions
	upper_results$test_data$yhat <- upper_preds$predictions
	
	lower_results$test_data <- lower_results$test_data %>%
		mutate(residuals = yhat - !!sym(dependent_variable),
					 group = "Lower 10%")
	
	upper_results$test_data <- upper_results$test_data %>%
		mutate(residuals = yhat - !!sym(dependent_variable),
					 group = "Upper 10%")
	
	residual_data <- bind_rows(lower_results$test_data, upper_results$test_data)
	
	plot <- ggplot(pdp_data, aes(x = standage, y = yhat, color = group, shape = group)) +
		geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, linewidth = .7, aes(fill = group)) + 
		geom_line(linewidth = .4) +
		geom_jitter(data = residual_data, aes(x = standage, y = yhat + residuals), width = 0.3, height = 0, alpha = 0.4) +
		scale_color_manual(values = c("Lower 10%" = "#0800af", "Upper 10%" = "#c82300")) +
		scale_shape_manual(values = c(16, 18)) + 
		labs(title = "Resource use strategy",
				 subtitle = paste("Lower 10% R\u00B2:", round(lower_r2, 2), "| Upper 10% R\u00B2:", round(upper_r2, 2)),
				 x = "Standage", y = "Resource use score") + 
		theme_bw() +
		ylim(0, 1) +
		facet_wrap(~group, ncol = 2, nrow = 1) +
		theme(
			legend.position = "none",
			legend.title = element_blank(),
			legend.text = element_text(size = 10),
			text = element_text(family = "Arial"),
			plot.title = element_text(face = "bold"),
			plot.subtitle = element_text(hjust = 0.5) 
		)
	
	r2_tibble <- tibble(group = c("Lower 10%", "Upper 10%"),
											r2 = c(lower_r2, upper_r2))
	
	return(list(plot = plot, r2_tibble = r2_tibble))
}

result <- plot_partial_dependence_climate(lower_results, upper_results, dependent_variable)
pdp_plot_climate <- result$plot
r2_values <- bind_rows(result$r2_tibble)

print(pdp_plot_climate) # examine figure
ggsave(filename = paste0(path_out, "/plots/climate_quantiles_pdp.png"),
			 plot = pdp_plot_climate,
			 bg = "white",
			 width = 250,
			 height = 150,
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

managed_results_0 <- train_and_predict(stratified_data_managed$managed_data_0, dependent_variable, covariates, hyper_grid)
managed_results_1 <- train_and_predict(stratified_data_managed$managed_data_1, dependent_variable, covariates, hyper_grid)

# Retrieve performance metrics for managed levels
managed_performance_metrics <- bind_rows(
	managed_results_0$quant_performance_metrics %>% mutate(group = "Managed 0"),
	managed_results_1$quant_performance_metrics %>% mutate(group = "Managed 1")
)

# Create partial dependence plot for managed variable with R-squared values
plot_partial_dependence_managed <- function(managed_results_0, managed_results_1, dependent_variable) {
	
	managed_model_0 <- managed_results_0$best_model
	managed_model_1 <- managed_results_1$best_model
	
	pdp_managed_0 <- pdp::partial(managed_model_0, pred.var = "standage", train = managed_results_0$test_data)
	pdp_managed_1 <- pdp::partial(managed_model_1, pred.var = "standage", train = managed_results_1$test_data)
	
	pdp_managed_0$group <- "Managed 0"
	pdp_managed_1$group <- "Managed 1"
	
	pdp_data <- bind_rows(pdp_managed_0, pdp_managed_1)
	
	managed_r2_0 <- managed_results_0$r_squared_value
	managed_r2_1 <- managed_results_1$r_squared_value
	
	managed_0_preds <- predict(managed_model_0, managed_results_0$test_data)
	managed_1_preds <- predict(managed_model_1, managed_results_1$test_data)
	
	managed_results_0$test_data$yhat <- managed_0_preds$predictions
	managed_results_1$test_data$yhat <- managed_1_preds$predictions
	
	managed_results_0$test_data <- managed_results_0$test_data %>%
		mutate(residuals = yhat - !!sym(dependent_variable),
					 group = "Managed 0")
	
	managed_results_1$test_data <- managed_results_1$test_data %>%
		mutate(residuals = yhat - !!sym(dependent_variable),
					 group = "Managed 1")
	
	residual_data <- bind_rows(managed_results_0$test_data, managed_results_1$test_data)
	
	plot <- ggplot(pdp_data, aes(x = standage, y = yhat, color = group, shape = group)) +
		geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, linewidth = .7, aes(fill = group)) + 
		#geom_line(linewidth = .4) +
		geom_jitter(data = residual_data, aes(x = standage, y = yhat + residuals), width = 0.3, height = 0, alpha = 0.4) +
		scale_color_manual(values = c("Managed 0" = "black", "Managed 1" = "firebrick4")) +
		scale_shape_manual(values = c("Managed 0" = 16, "Managed 1" = 18)) + 
		labs(title = "Resource use strategy",
				 subtitle = paste("Managed 0 R\u00B2:", round(managed_r2_0, 2), "| Managed 1 R\u00B2:", round(managed_r2_1, 2)),
				 x = "Standage", y = "Resource use score") + 
		theme_bw() +
		ylim(0, 1) +
		theme(
			legend.position = "none",
			legend.title = element_blank(),
			legend.text = element_text(size = 10),
			text = element_text(family = "Arial"),
			plot.title = element_text(face = "bold"),
			plot.subtitle = element_text(hjust = 0.5) 
		) +
		facet_wrap(~group)
	
	r2_tibble <- tibble(group = c("Managed 0", "Managed 1"),
											r2 = c(managed_r2_0, managed_r2_1))
	
	return(list(plot = plot, r2_tibble = r2_tibble))
}

# Generate plot for partial dependence
result <- plot_partial_dependence_managed(managed_results_0, managed_results_1, dependent_variable)
pdp_plot_managed <- result$plot

print(pdp_plot_managed)
ggsave(filename = paste0(path_out, "/plots/managed_pdp.png"),
			 plot = pdp_plot_managed,
			 bg = "white",
			 width = 250,
			 height = 200,
			 units = "mm",
			 dpi = 1457)

stopCluster(cl) # stop the parallel backend
