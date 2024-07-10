############################################################################################################################
########################################  MSc Diss. Forest Succession Data Analysis ########################################  
############################################################################################################################

# Objectives:
# 	1.	Understand the impact of stand age (and other features) on functional diversity indices using Shapley values.
#   2.	Predict how index patterns across successional time vary across different climate conditions, management, and ecoregions.

# Outline of the Approach:
# 	-	Data Preparation
#   -	Hyperparameter Tuning using grid search combined with cross-validation
# 	-	Model Fitting with Optimal Hyperparameters
#   -	Interpretation with Shapley Values
#   -	Analysis Across Ecoregions, Management and Climate Conditions using strata based on quantiles and factors

rm(list = ls()) # make sure environment is clean 
set.seed(42)    # set seed for reproducibility
Sys.setenv(R_FORK_SUPPORTED = "false") # disable forking for memory management

suppressWarnings({
	library(doParallel)
	registerDoParallel()
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
data <- read_csv("/Volumes/ritd-ag-project-rd01pr-dmayn10/merlin/data/fia_traits/sub/fun_rf_clean.csv",
								 col_types = c("d","d","d","d","d","d","d","d","d","d","f","f","d","d","d","d","d","d","d","d","d","d","d","d")) %>%
	rename_with(~ gsub("wmean_", "", .), starts_with("wmean_")) %>%
	rename_with(~ gsub("_015cm", "", .), ends_with("_015cm")) %>%
	rename(min_temperature = min_temperature_of_coldest_month,
				 max_temperature = max_temperature_of_warmest_month)

data <- data %>% 
	# We filter Standage by the upper 10% quantiles
	filter(standage < quantile(standage, 0.9)) %>%
	sample_n(n() * 0.2) # use a 20% sample for code development 

# Split into training and test sets
split <- initial_split(data, prop = 0.8)

# Extract training and testing datasets
train_data <- training(split)
test_data  <- testing(split)

# Define the indices and covariates
indices <- c("fun_div", "fun_disp", "fun_even")

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
tune_rf_model <- function(index, data, covariates, hyper_grid) {
	
	# Define model formula
	formula <- as.formula(paste(index, "~", paste(covariates, collapse = " + ")))
	
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
			index = index, 
			best = rowSums(sapply(names(best_tune), 
														function(param) results[[param]] == best_tune[[param]])) == length(best_tune)
		)
	
	return(list(model = model, results = results))
}

# Register parallel cores before tuning and calculating Shapley values
num_cores <- detectCores() - 2
registerDoParallel(num_cores)

# Perform hyperparameter tuning for each index and get performance metrics
tuned_models <- indices %>%
	map(~ tune_rf_model(.x, train_data, covariates, hyper_grid))

# Extract the best model for every index
best_models <- tuned_models %>%
	map(~ .x$model$finalModel)

# Extract performance metrics
global_performance_metrics <- tuned_models %>%
	map_df(~ .x$results) %>%
	mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
	arrange(index, best)

# Prediction function for fastshap
predict_fn <- function(object, newdata) {
	predict(object, data = newdata)$predictions
}

## ----- Calculate Shapley values for each index model using parallel processing -----

shap_values <- foreach(i = seq_along(best_models), .combine = 'rbind', .packages = c('fastshap', 'dplyr')) %dopar% {
	model <- best_models[[i]]
	index <- indices[i]
	shap_values <- fastshap::explain(model, X = test_data %>% select(all_of(covariates)) %>% as.matrix(),
																	 pred_wrapper = predict_fn, nsim = 10, parallel = TRUE)
	shap_values <- as.data.frame(shap_values)
	shap_values$index <- index
	shap_values
}

# Stop the parallel backend
stopImplicitCluster()

## Visualize shapley values... first get a tibble for plotting
shap_long <- shap_values %>%
	pivot_longer(cols = -c(index), names_to = "feature", values_to = "shap") %>%
	as_tibble() %>%
	mutate(feature = gsub("_", " ", feature),
				 index = case_when(
				 	index == "fun_disp" ~ "Functional Dispersion",
				 	index == "fun_div" ~ "Functional Divergence",
				 	index == "fun_even" ~ "Functional Evenness",
				 	TRUE ~ NA_character_)) 

shap_plot <- shap_long %>%
	ggplot(aes(x = feature, y = shap, color = shap)) +
	geom_quasirandom(alpha = 0.5) +
	facet_wrap(~index, ncol = 3, nrow = 1, scale = "free_x") +
	scale_color_viridis_c(option = "viridis",
												name = "Feature\nValue",
												breaks = c(min(shap_long$shap), max(shap_long$shap)),
												labels = c("low", "high")) +
	coord_flip() +
	labs(x = NULL, y = "Shapley Value") +
	theme_bw() +
	theme(text = element_text(family = "Arial"),
				legend.position = "top", 
				strip.background = element_rect(fill = "white", color = "black", linewidth = .75),
				strip.text = element_text(color = "black")) 

print(shap_plot) # examine plot 
ggsave(filename = "/Users/serpent/Documents/MSc/Thesis/Code/analysis/indices/plots/shap_plots.png",
			 plot = shap_plot,
			 bg = "white",
			 width = 250,
			 height = 175,
			 units = "mm",
			 dpi = 1457)

## Calculate and plot feature importance 

# Sum absolute Shapley values to determine overall importance
feature_importance <- shap_long %>%
	group_by(index, feature) %>%
	summarize(importance = sum(abs(shap)), .groups = "drop") %>%
	group_by(index) %>%
	mutate(avg_importance = mean(importance)) %>%
	ungroup()

# Create the plot
importance_plot <- ggplot(feature_importance, aes(x = reorder(feature, importance), y = importance, fill = index)) +
	geom_bar(stat = "identity", colour = "black", alpha = .7) +
	geom_hline(aes(yintercept = avg_importance), linetype = "dashed", colour = "red") +
	facet_wrap(~index, scales = "free_x", ncol = 3, nrow = 1) +
	coord_flip() +
	scale_fill_viridis_d()+
	labs(title = "Feature Importance by Index",
			 x = "Feature",
			 y = "Importance",
			 fill = "Index",
			 color = "Average Importance") +
	theme_bw(base_size = 15) +
	theme(legend.position = "none",
				text = element_text(family = "Arial"),
				plot.title = element_text(hjust = 0.5),
				plot.subtitle = element_text(hjust = 0.5))

print(importance_plot) # examine plot 
ggsave(filename = "/Users/serpent/Documents/MSc/Thesis/Code/analysis/indices/plots/importance_plot.png",
			 plot = importance_plot,
			 bg = "white",
			 width = 250,
			 height = 200,
			 units = "mm",
			 dpi = 1457)

## Plot relationship of traits with standage 

plots <- list()
for (index in unique(shap_long$index)) {
	
	index_data <- shap_long %>% filter(index == !!index)
	feature_data <- index_data %>% filter(feature == "standage")
	observed_values <- test_data %>% pull(standage)
	
	plot <- ggplot(feature_data, aes(x = observed_values, y = shap, color = shap)) +
		geom_beeswarm(alpha = 0.5, colour = "darkgreen") +
		geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.5) +
		labs(x = "standage", y = "Influence (Shapley)") +
		theme_bw() +
		theme(text = element_text(family = "Arial"),
					legend.position = "right",
					legend.key.width = unit(0.5, "cm"),
					legend.key.height = unit(2, "cm"),
					legend.box.background = element_rect(color = "black", linewidth = .75),
					plot.title = element_text(hjust = 0.5)) +
		ggtitle(paste(index))
	
	plots[[index]] <- plot
}

# Combine all plots into one figure with 2 columns and 4 rows
shap_standage_plots <- wrap_plots(plots, ncol = 3, nrow = 1)
print(shap_standage_plots)

ggsave(filename = "/Users/serpent/Documents/MSc/Thesis/Code/analysis/indices/plots/shap_standage_plot.png",
			 plot = shap_standage_plots,
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
train_and_predict <- function(data, indices, covariates, hyper_grid) {
	
	# Split into training and test sets
	set.seed(42)
	split <- initial_split(data, prop = 0.8)
	
	# Extract training and testing datasets
	train_data <- training(split)
	test_data  <- testing(split)
	
	# Perform hyperparameter tuning for each trait and get performance metrics
	tuned_models <- indices %>% map(~ tune_rf_model(.x, train_data, covariates, hyper_grid))
	
	performance_metrics <- tuned_models %>%
		map_df(~ .x$results) %>%
		mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
		arrange(index, best) %>%
		as_tibble()
	
	# Extract the best model for every trait
	best_models <- tuned_models %>% map(~ .x$model$finalModel)
	
	# Calculate R-squared for each model
	r_squared_values <- best_models %>% map2(indices, ~ caret::R2(predict(.x, test_data)$predictions, test_data[[.y]]))
	
	return(list(best_models = best_models, test_data = test_data, r_squared_values = r_squared_values, quant_performance_metrics = performance_metrics))
}

# Perform stratification and model training
stratified_data_climate <- stratify_climate(data, "annual_mean_temperature")

# Register parallel cores before tuning
registerDoParallel(num_cores)

lower_results <- train_and_predict(stratified_data_climate$lower_data, indices, covariates, hyper_grid)
upper_results <- train_and_predict(stratified_data_climate$upper_data, indices, covariates, hyper_grid)

# Stop the parallel backend
stopImplicitCluster()

# Retrieve performance metrics for lower and upper quantiles
quant_performance_metrics <- bind_rows(
	lower_results$quant_performance_metrics %>% mutate(group = "Lower 10%"),
	upper_results$quant_performance_metrics %>% mutate(group = "Upper 10%"))

# Function to create individual partial dependence plots with R-squared values and custom colours

index_title <- function(index) { # helper function for plot titles
	case_when(
			index == "fun_disp" ~ "Functional Dispersion",
			index == "fun_div" ~ "Functional Divergence",
			index == "fun_even" ~ "Functional Evenness",
			TRUE ~ NA_character_)
}

# Function to create individual partial dependence plots with R-squared values and custom colours
plot_partial_dependence_climate <- function(lower_results, upper_results, index) {
	
	# Extract models
	lower_model <- lower_results$best_models[[which(indices == index)]]
	upper_model <- upper_results$best_models[[which(indices == index)]]
	
	# Generate partial dependence data
	pdp_lower <- pdp::partial(lower_model, pred.var = "standage", train = lower_results$test_data)
	pdp_upper <- pdp::partial(upper_model, pred.var = "standage", train = upper_results$test_data)
	
	# Add quantile identifiers
	pdp_lower$group <- "Lower 10%"
	pdp_upper$group <- "Upper 10%"
	
	# Combine partial dependence data
	pdp_data <- bind_rows(pdp_lower, pdp_upper)
	
	# Extract R-squared values
	lower_r2 <- lower_results$r_squared_values[[which(indices == index)]]
	upper_r2 <- upper_results$r_squared_values[[which(indices == index)]]
	
	# Calculate residuals for jitter
	lower_preds <- predict(lower_model, lower_results$test_data)
	upper_preds <- predict(upper_model, upper_results$test_data)
	
	lower_results$test_data$yhat <- lower_preds$predictions
	upper_results$test_data$yhat <- upper_preds$predictions
	
	lower_results$test_data <- lower_results$test_data %>%
		mutate(residuals = yhat - !!sym(index),
					 group = "Lower 10%")
	
	upper_results$test_data <- upper_results$test_data %>%
		mutate(residuals = yhat - !!sym(index),
					 group = "Upper 10%")
	
	residual_data <- bind_rows(lower_results$test_data, upper_results$test_data)
	
	plot <- ggplot(pdp_data, aes(x = standage, y = yhat, color = group, shape = group)) +
		geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, linewidth = .7, aes(fill = group)) + 
		geom_line(linewidth = .4) +
		geom_jitter(data = residual_data, aes(x = standage, y = yhat + residuals), width = 0.3, height = 0, alpha = 0.4) +
		scale_color_manual(values = c("Lower 10%" = "#0800af", "Upper 10%" = "#c82300")) +
		scale_shape_manual(values = c(16, 18)) + 
		labs(title = index_title(index),
				 subtitle = paste("Lower 10% R\u00B2:", round(lower_r2, 2), "| Upper 10% R\u00B2:", round(upper_r2, 2)),
				 x = "Standage", y = "Index value") + 
		theme_bw() +
		theme(
			legend.position = "none",
			legend.title = element_blank(),
			legend.text = element_text(size = 10),
			text = element_text(family = "Palatino"),
			plot.title = element_text(face = "bold"),
			plot.subtitle = element_text(hjust = 0.5)
		) +
		facet_wrap(~group, scales = "fixed")
	
	r2_tibble <- tibble(index = index,
											group = c("Lower 10%", "Upper 10%"),
											r2 = c(lower_r2, upper_r2))
	
	return(list(plot = plot, r2_tibble = r2_tibble))
}

# Generate plot data for all traits
pdp_plots_climate <- list()
r2_values <- tibble(index = character(), group = character(), r2 = numeric())
for (index in indices) {
	result <- plot_partial_dependence_climate(lower_results, upper_results, index)
	pdp_plots_climate[[index]] <- result$plot
	r2_values <- bind_rows(r2_values, result$r2_tibble)
}

# Combine the pdp plots using patchwork
pdp_plots_climate <- wrap_plots(pdp_plots_climate, ncol = 3, nrow = 1) +
	plot_annotation(
		theme = theme(
			text = element_text(family = "Palatino"),
			plot.title = element_text(face = "bold")
		)
	)

# Display the pdp plot
print(pdp_plots_climate) # examine figure
ggsave(filename = "/Users/serpent/Documents/MSc/Thesis/Code/analysis/indices/plots/clim_plots.png",
			 plot = pdp_plots_climate,
			 bg = "white",
			 width = 250,
			 height = 150,
			 units = "mm",
			 dpi = 1457)

# Create bar chart for R-squared values
r2_plot <- ggplot(r2_values, aes(x = index, y = r2, fill = group)) +
	geom_bar(stat = "identity", position = position_dodge(width = 0.5), colour = "black") +
	scale_fill_manual(values = c("Lower 10%" = "#0800af", "Upper 10%" = "#c82300")) +
	labs(x = NULL, y = "R²") +
	guides(fill = guide_legend(title = "Annual mean temperature quantiles:")) +
	theme_bw() +
	theme(
		text = element_text(family = "Palatino"),
		plot.title = element_text(face = "bold"),
		axis.text.x = element_text(angle = 45, hjust = 1),
		legend.position = "top",
		legend.justification = "left",
		legend.background = element_rect(fill = "transparent"))

# Display the bar chart
print(r2_plot) # examine R-squared values
ggsave(filename = "/Users/serpent/Documents/MSc/Thesis/Code/analysis/indices/plots/r2_plot.png",
			 plot = r2_plot,
			 bg = "white",
			 width = 250,
			 height = 200,
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

managed_results_0 <- train_and_predict(stratified_data_managed$managed_data_0, indices, covariates, hyper_grid)
managed_results_1 <- train_and_predict(stratified_data_managed$managed_data_1, indices, covariates, hyper_grid)

# Stop the parallel backend
stopImplicitCluster()

# Retrieve performance metrics for managed levels
managed_performance_metrics <- bind_rows(
	managed_results_0$quant_performance_metrics %>% mutate(group = "Managed 0"),
	managed_results_1$quant_performance_metrics %>% mutate(group = "Managed 1")
)

# Function to create partial dependence plots for managed variable with R-squared values and custom colors
plot_partial_dependence_managed <- function(managed_results_0, managed_results_1, index) {
	managed_model_0 <- managed_results_0$best_models[[which(indices == index)]]
	managed_model_1 <- managed_results_1$best_models[[which(indices == index)]]
	
	pdp_managed_0 <- pdp::partial(managed_model_0, pred.var = "standage", train = managed_results_0$test_data)
	pdp_managed_1 <- pdp::partial(managed_model_1, pred.var = "standage", train = managed_results_1$test_data)
	
	pdp_managed_0$group <- "Managed 0"
	pdp_managed_1$group <- "Managed 1"
	
	pdp_data <- bind_rows(pdp_managed_0, pdp_managed_1)
	
	managed_r2_0 <- managed_results_0$r_squared_values[[which(indices == index)]]
	managed_r2_1 <- managed_results_1$r_squared_values[[which(indices == index)]]
	
	# Calculate residuals for jitter
	managed_preds_0 <- predict(managed_model_0, managed_results_0$test_data)
	managed_preds_1 <- predict(managed_model_1, managed_results_1$test_data)
	
	managed_results_0$test_data$yhat <- managed_preds_0$predictions
	managed_results_1$test_data$yhat <- managed_preds_1$predictions
	
	managed_results_0$test_data <- managed_results_0$test_data %>%
		mutate(residuals = yhat - !!sym(index),
					 group = "Managed 0")
	
	managed_results_1$test_data <- managed_results_1$test_data %>%
		mutate(residuals = yhat - !!sym(index),
					 group = "Managed 1")
	
	residual_data <- bind_rows(managed_results_0$test_data, managed_results_1$test_data)
	
	plot <- ggplot(pdp_data, aes(x = standage, y = yhat, color = group, shape = group)) +
		geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, linewidth = .7, aes(fill = group)) + 
		geom_line(linewidth = .4) +
		geom_jitter(data = residual_data, aes(x = standage, y = yhat + residuals), width = 0.3, height = 0, alpha = 0.4) +  # Add jittered residuals
		scale_color_manual(values = c("Managed 0" = "#0800af", "Managed 1" = "#c82300")) +
		scale_shape_manual(values = c(16, 18)) + 
		labs(title = index_title(index),
				 subtitle = paste("Managed 0 R\u00B2:", round(managed_r2_0, 2), "| Managed 1 R\u00B2:", round(managed_r2_1, 2)),
				 x = "Standage", y = "Index value") + 
		theme_bw() +
		theme(
			legend.position = "none",
			legend.text = element_text(size = 10),
			text = element_text(family = "Palatino"),
			plot.title = element_text(face = "bold"),
			plot.subtitle = element_text(hjust = 0.5)
		) +
		facet_wrap(~group, scales = "fixed")
	
	r2_tibble <- tibble(index = index,
											group = c("Managed 0", "Managed 1"),
											r2 = c(managed_r2_0, managed_r2_1))
	
	return(list(plot = plot, r2_tibble = r2_tibble))
}

# Generate plots for all traits
pdp_plots_managed <- list()
r2_values_managed <- tibble(index = character(), group = character(), r2 = numeric())
for (index in indices) {
	result <- plot_partial_dependence_managed(managed_results_0, managed_results_1, index)
	pdp_plots_managed[[index]] <- result$plot
	r2_values_managed <- bind_rows(r2_values_managed, result$r2_tibble)
}

# Combine the plots using patchwork
managed_plot <- wrap_plots(pdp_plots_managed, ncol = 3, nrow = 1) +
	plot_annotation(
		theme = theme(
			text = element_text(family = "Palatino"),
			plot.title = element_text(face = "bold")
		)
	)

# Display the combined plot
print(managed_plot) # examine figure
ggsave(filename = "/Users/serpent/Documents/MSc/Thesis/Code/analysis/indices/plots/managed_plots.png",
			 plot = managed_plot,
			 bg = "white",
			 width = 150,
			 height = 250,
			 units = "mm",
			 dpi = 1457)

# Create bar chart for R-squared values
r2_plot_managed <- ggplot(r2_values_managed, aes(x = index, y = r2, fill = group)) +
	geom_bar(stat = "identity", position = position_dodge(width = 0.5), colour = "black") +
	scale_fill_manual(values = c("Managed 0" = "#0800af", "Managed 1" = "#c82300")) +
	labs(x = NULL, y = "R²") +
	guides(fill = guide_legend(title = "Management status:")) +
	theme_bw() +
	theme(
		text = element_text(family = "Palatino"),
		plot.title = element_text(face = "bold"),
		axis.text.x = element_text(angle = 45, hjust = 1),
		legend.position = "top",
		legend.justification = "left",
		legend.background = element_rect(fill = "transparent")
	)

# Display the bar chart
print(r2_plot_managed) # examine R-squared values
ggsave(filename = "/Users/serpent/Documents/MSc/Thesis/Code/analysis/traits/plots/r2_plot_managed.png",
			 plot = r2_plot_managed,
			 bg = "white",
			 width = 250,
			 height = 200,
			 units = "mm",
			 dpi = 1457)

##### Stratify data based on biome and make predictions 