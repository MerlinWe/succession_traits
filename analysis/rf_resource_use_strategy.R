############################################################################################################################
########################################  MSc Diss. Forest Succession Data Analysis ########################################  
############################################################################################################################

rm(list = ls())
options(scipen = 999)
set.seed(42)

library(car)
library(caret)
library(MASS)
library(randomForest)
library(ranger)
library(stringi)
library(tidyverse)

## Read most recent data
fia_clean <- list.files("/Volumes/ritd-ag-project-rd01pr-dmayn10/merlin/data/fia_traits", full.names = TRUE) %>%
	file.info() %>%
	as_tibble(rownames = "file") %>%
	arrange(desc(mtime)) %>%
	slice(1) %>%
	pull(file) %>%
	read_csv() %>% 
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
	select(-Name) #

## >>>>>>>>>>>>>>>>>>>> Predicting resource use strategy using random forest <<<<<<<<<<<<<<<<<<<<

# Get input data for fitting a random forest 
strat <- fia_clean %>% 
	filter(standage < 500) %>%
	
	select(resource_use_score, standage, ecoregion, managed, 
				 annual_mean_temperature, annual_precipitation, temperature_seasonality, mean_diurnal_range, 
				 min_temperature_of_coldest_month, max_temperature_of_warmest_month,
				 elevation, pop_density, sand_content_015cm, soil_ph_015cm, water_capacity_015cm) %>%

	filter(resource_use_score < 1 & resource_use_score > 0) %>% 
	filter(complete.cases(.)) %>%
	as_tibble() 

table(strat$managed)
table(strat$ecoregion)

## Random forest regression workflow 
set.seed(42)

# withhold 10% for final testing; get the test ids at random
test_ids <- sample(1:nrow(strat), round(.10*nrow(strat)))
test_data <- strat %>% slice(test_ids)

# the rest of the data go into testing/training
train_data <- strat %>% slice(-test_ids)

## ----- 1. Initial Model Fit (as a Baseline) -----

# Initial model fit with default parameters
baseline_rf <- ranger(resource_use_score ~ standage + managed +
		annual_mean_temperature + annual_precipitation + temperature_seasonality +
		mean_diurnal_range + min_temperature_of_coldest_month + max_temperature_of_warmest_month +
		elevation + pop_density + sand_content_015cm + soil_ph_015cm + water_capacity_015cm,
		data = train_data, importance = "impurity")

# Evaluate the initial model on the test set
initial_predictions <- predict(baseline_rf, data = test_data)
initial_r2 <- cor(log(initial_predictions$predictions), log(test_data$resource_use_score))^2
print(initial_r2) # quite bad 

## ----- 2. Hyperparameter Tuning with Cross-Validation -----

# Define the grid of hyperparameters to search
mtry_values <- c(2,5,10)
ntree_values <- c(100, 200, 500)
nodesize_values <- c(1, 5, 10)

# Create a tibble to store the results
results <- expand.grid(mtry = mtry_values, ntree = ntree_values, nodesize = nodesize_values)
results$mean_r2 <- NA

# Specify the number of folds
nfolds <- 10

# Perform cross-validation for each combination of hyperparameters
for(i in 1:nrow(results)) {
	mtry <- results$mtry[i]
	ntree <- results$ntree[i]
	nodesize <- results$nodesize[i]
	
	fold_r2 <- rep(NA, nfolds)
	fold_id <- sample(1:nfolds, nrow(train_data), replace = TRUE)
	
	for(j in 1:nfolds) {
		val_data <- train_data[fold_id == j,]
		train_fold_data <- train_data[fold_id != j,]
		
		model <- ranger(resource_use_score ~ standage + managed +
											annual_mean_temperature + annual_precipitation + temperature_seasonality +
											mean_diurnal_range + min_temperature_of_coldest_month + max_temperature_of_warmest_month +
											elevation + pop_density + sand_content_015cm + soil_ph_015cm + water_capacity_015cm,
										data = train_fold_data,
										mtry = mtry,
										num.trees = ntree,
										min.node.size = nodesize)

		
		val_predict <- predict(model, data = val_data)
		fold_r2[j] <- cor(log(val_predict$predictions), log(val_data$resource_use_score))^2
	}
	
	results$mean_r2[i] <- mean(fold_r2)
}

# Find the best hyperparameters
best_params <- results[which.max(results$mean_r2),]

# Print the best hyperparameters
print(best_params)

## ---------- 3.Final Model Training with Best Hyperparameters ----------

# Fit the final model using the best hyperparameters
final_rf <- ranger(resource_use_score ~ standage + managed +
									 	annual_mean_temperature + annual_precipitation + temperature_seasonality +
									 	mean_diurnal_range + min_temperature_of_coldest_month + max_temperature_of_warmest_month +
									 	elevation + pop_density + sand_content_015cm + soil_ph_015cm + water_capacity_015cm,
	data = train_data,
	mtry = best_params$mtry,
	num.tree = best_params$ntree,
	min.node.size = best_params$nodesize,
	importance = "impurity"
)

# Plot variable importance (Percentage Increase in Mean Squared Error and Increase in Node Purity)
var_importance <-  importance(final_rf) %>%
	as_tibble(rownames = "variable") %>%
	arrange(desc(value))

# Evaluate the final model on the test set
final_predictions <- predict(final_rf, data = test_data)
cat("r2 (final model) =", cor(log(final_predictions$predictions), log(test_data$resource_use_score))^2) #
cat("accuracy (final model) =", mean((final_predictions$predictions - strat$resource_use_score)^2)) 

# Plot Predictions on the input data
tibble(Observed = strat$resource_use_score,
			 Predicted = predict(final_rf, data = strat)$predictions) %>%
	ggplot(aes(x = Observed, y = Predicted)) +
	geom_point(alpha = .5) +
	geom_abline(slope = 1, intercept = 0, color = "red") +
	theme_bw() +
	labs(title = 'Observed vs Predicted Resource Use Score',
			 x = 'Observed Resource Use Score',
			 y = 'Predicted Resource Use Score')
