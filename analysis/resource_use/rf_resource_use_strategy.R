############################################################################################################################
########################################  MSc Diss. Forest Succession Data Analysis ########################################  
############################################################################################################################

rm(list = ls())
options(scipen = 999)
set.seed(42)

library(caret)
library(ranger)
library(pdp)
library(fastshap)
library(doParallel)
library(ggbeeswarm)
library(tidyverse)

## Read most recent data
fia_clean <- read_csv("/Volumes/ritd-ag-project-rd01pr-dmayn10/merlin/data/fia_traits/plotlvl_data_2024-06-25.csv") %>% 
	mutate(
		PID = as.character(PID),
		PID_rep = as.character(PID_rep),
		rep_measure = as.logical(rep_measure),
		PID_measure = as.integer(PID_measure),
		state = as.factor(state),
		standage = as.integer(standage),
		INVYR = as.integer(INVYR),
		FORTYPCD = as.integer(FORTYPCD),
		foresttype = as.factor(foresttype),
		biome = as.factor(biome),
		ownership = as.factor(ownership),
		managed = as.factor(managed),
		ll_id = as.character(ll_id),
		across(starts_with("wmean_"), as.numeric),
		across(starts_with("total_"), as.numeric),
		across(starts_with("fun_"), as.numeric),
		resource_use_score = as.numeric(resource_use_score),
		annual_mean_temperature = as.numeric(annual_mean_temperature),
		annual_precipitation = as.numeric(annual_precipitation),
		isothermality = as.numeric(isothermality),
		max_temperature_of_warmest_month = as.numeric(max_temperature_of_warmest_month),
		mean_diurnal_range = as.numeric(mean_diurnal_range),
		mean_temperature_of_coldest_quarter = as.numeric(mean_temperature_of_coldest_quarter),
		mean_temperature_of_driest_quarter = as.numeric(mean_temperature_of_driest_quarter),
		mean_temperature_of_warmest_quarter = as.numeric(mean_temperature_of_warmest_quarter),
		mean_temperature_of_wettest_quarter = as.numeric(mean_temperature_of_wettest_quarter),
		min_temperature_of_coldest_month = as.numeric(min_temperature_of_coldest_month),
		precipitation_seasonality = as.numeric(precipitation_seasonality),
		precipitation_of_coldest_quarter = as.numeric(precipitation_of_coldest_quarter),
		precipitation_of_driest_month = as.numeric(precipitation_of_driest_month),
		precipitation_of_driest_quarter = as.numeric(precipitation_of_driest_quarter),
		precipitation_of_warmest_quarter = as.numeric(precipitation_of_warmest_quarter),
		precipitation_of_wettest_month = as.numeric(precipitation_of_wettest_month),
		precipitation_of_wettest_quarter = as.numeric(precipitation_of_wettest_quarter),
		temperature_annual_range = as.numeric(temperature_annual_range),
		temperature_seasonality = as.numeric(temperature_seasonality),
		elevation = as.numeric(elevation),
		pop_density = as.numeric(pop_density),
		ecoregion = as.integer(ecoregion),
		across(starts_with("sand_content_"), as.numeric),
		across(starts_with("soil_ph_"), as.numeric),
		across(starts_with("water_capacity_"), as.numeric))

# Add ecoregions with names 
ecoregions <- read_csv("/Volumes/ritd-ag-project-rd01pr-dmayn10/merlin/data/composite/resolve_ecoregions_legend.csv", locale = locale(encoding = "UTF-8")) %>%
	mutate(
		Name = iconv(Name, from = "UTF-8", to = "ASCII//TRANSLIT"), 
		Name = stri_replace_all_regex(tolower(Name), "\\s", "_")) %>%
	mutate(Code = as.integer(Code))
fia_clean <- fia_clean %>%
	left_join(ecoregions, by = c("ecoregion" = "Code")) %>%
	mutate(ecoregion = Name) %>%
	dplyr::select(-Name) #

## >>>>>>>>>>>>>>>>>>>> Predicting resource use strategy using random forest <<<<<<<<<<<<<<<<<<<<

# Get input data for fitting a random forest 
data <- fia_clean %>% 
	filter(standage < 500) %>%
	sample_n(1000) %>%
	
	dplyr::select(resource_use_score, standage, ecoregion, managed, 
				 annual_mean_temperature, annual_precipitation, temperature_seasonality, mean_diurnal_range, 
				 min_temperature_of_coldest_month, max_temperature_of_warmest_month,
				 elevation, pop_density, sand_content_015cm, soil_ph_015cm, water_capacity_015cm) %>%

	filter(resource_use_score < 1 & resource_use_score > 0) %>% 
	filter(complete.cases(.)) %>%
	as_tibble() 

# Split into training and test sets
set.seed(42)

train_indices <- sample(1:nrow(data), size = 0.8 * nrow(data))
train_data <- data[train_indices, ]
test_data <- data[-train_indices, ]

# Define the dependent variable and covariates
dependent_variable <- "resource_use_score"

covariates <- c("standage", "annual_mean_temperature", "annual_precipitation", 
								"temperature_seasonality", "mean_diurnal_range", 
								"min_temperature_of_coldest_month", "max_temperature_of_warmest_month", 
								"elevation", "pop_density", "sand_content_015cm", 
								"soil_ph_015cm", "water_capacity_015cm")

# Define a grid of hyperparameters including mtry, min.node.size, and split rule
hyper_grid <- expand.grid(
	mtry = seq(2, length(covariates), by = 2),
	splitrule = "variance",
	min.node.size = c(1, 5, 10)
)

# Function to perform cross-validation, find the best hyperparameters, and capture performance metrics
tune_rf_model <- function(data, dependent_variable, covariates, hyper_grid) {
	formula <- as.formula(paste(dependent_variable, "~", paste(covariates, collapse = " + ")))
	
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
		mutate(best = rowSums(sapply(names(best_tune), function(param) results[[param]] == best_tune[[param]])) == length(best_tune))
	
	return(list(model = model, results = results))
}

# Perform hyperparameter tuning and capture performance metrics
tuned_model <- tune_rf_model(train_data, dependent_variable, covariates, hyper_grid)

# Extract the best model and performance metrics
best_model <- tuned_model$model$finalModel
performance_metrics <- tuned_model$results

# Define prediction function for fastshap
predict_fn <- function(object, newdata) {
	predict(object, data = newdata)$predictions
}

# Calculate Shapley values
registerDoParallel(detectCores() - 1)
shap_values <- fastshap::explain(best_model, X = test_data[, covariates, drop = FALSE], pred_wrapper = predict_fn, nsim = 100, parallel = TRUE)

# Define a function to plot Shapley values
plot_shap_values <- function(shap_values, covariates, dependent_variable) {
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
		labs(title = paste(dependent_variable), x = "Feature", y = "Shapley Value") +
		theme_bw() +
		theme(text = element_text(family = "Palatino"),
					legend.position = "none")
}

plot_shap_values(shap_values, covariates, dependent_variable)

# Plot predictions using partial residual curves
plot_partial_dependence <- function(model, data, dependent_variable, covariate) {
	
	# Generate partial dependence data
	pd_data <- partial(model, pred.var = covariate, train = data, grid.resolution = 100)
	
	# Convert to data frame for ggplot
	pd_df <- as.data.frame(pd_data)
	
	# Plot partial dependence data
	ggplot(pd_df, aes(x = !!sym(covariate), y = yhat)) +
		geom_point(color = "blue") +
		labs(title = paste("Partial Dependence of", dependent_variable, "on", covariate),
				 x = covariate, y = "Predicted Value") +
		theme_bw() +
		theme(text = element_text(family = "Palatino"))
}

# Generate and plot partial dependence
partial_plots <- lapply(covariates, function(covariate) {
	plot_partial_dependence(best_model, test_data, dependent_variable, covariate)
})

# Combine all plots into a single plot
combined_partial_plot <- cowplot::plot_grid(plotlist = partial_plots, ncol = 4, nrow = 3)
print(combined_partial_plot)

