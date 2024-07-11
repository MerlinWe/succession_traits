############################################################################################################################
########################################  MSc Diss. Forest Succession Data Analysis ########################################  
############################################################################################################################

# Objectives:
# 	1.	Understand the impact of stand age (and other features) on tree traits using Shapley values.
#   2.	Predict how trait patterns across successional time vary across different climate conditions and ecoregions.

# Outline of the Approach:
# 	-	Data Preparation
#   -	Hyperparameter Tuning using grid search combined with cross-validation
# 	-	Model Fitting with Optimal Hyperparameters
#   -	Interpretation with Shapley Values
#   -	Analysis Across Ecoregions and Climate Conditions using strata based on quantiles and factors

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
library(grid)
library(tidyverse)

# Set path to run either on local device 
path_in <- "/Volumes/ritd-ag-project-rd01pr-dmayn10/merlin/data/fia_traits/sub" 
path_out <- "/Users/serpent/Documents/MSc/Thesis/Code/analysis/traits" 

# Set paths to run on threadripper
path_in <- "/home/merlin/RDS_drive/merlin/data/fia_traits/sub"
path_out <- "/home/merlin/RDS_drive/merlin/Code/analysis/traits" 

# Read data and make some minor adjustments
data <- read_csv(paste0(path_in, "/traits_rf_clean.csv"),
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
# If on threadripper set to 32, if local, set to 10 !!
num_cores <- 10
cl <- makeCluster(num_cores)
registerDoParallel(cl)

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

# Prediction function for fastshap
predict_fn <- function(object, newdata) {
	predict(object, data = newdata)$predictions
}

## ----- Calculate Shapley values for each trait model using parallel processing -----

shap_values <- foreach(i = seq_along(best_models), .combine = 'rbind', .packages = c('fastshap', 'dplyr')) %dopar% {
	model <- best_models[[i]]
	trait <- traits[i]
	shap_values <- fastshap::explain(model, X = test_data %>% select(all_of(covariates)) %>% as.matrix(),
																	 pred_wrapper = predict_fn, nsim = 100, parallel = TRUE)
	shap_values <- as.data.frame(shap_values)
	shap_values$trait <- trait
	shap_values
}

## Visualize shapley values... first get a tibble for plotting
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
ggsave(filename = paste0(path_out, "/plots/global_shap_plots.png"),
			 plot = shap_plot,
			 bg = "white",
			 width = 280,
			 height = 200,
			 units = "mm",
			 dpi = 1457)

## Calculate and plot feature importance 

# Sum absolute Shapley values to determine overall importance
feature_importance <- shap_long %>%
	group_by(trait, feature) %>%
	summarize(importance = sum(abs(shap)), .groups = "drop") %>%
	arrange(trait, importance) %>%
	group_by(trait) %>%
	mutate(avg_importance = mean(importance)) %>%
	ungroup()

# Create the plot
importance_plot <- ggplot(feature_importance, aes(x = reorder(feature, importance), y = importance, fill = trait)) +
	geom_bar(stat = "identity", colour = "black", alpha = .7) +
	geom_hline(aes(yintercept = avg_importance), linetype = "dashed", colour = "red") +
	facet_wrap(~trait, scales = "fixed", ncol = 2, nrow = 4) +
	coord_flip() +
	scale_fill_viridis_d()+
	labs(x = NULL,
			 y = "Overall Feature Importance",
			 fill = "Trait",
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

plots <- list()
for (trait in unique(shap_long$trait)) {
	
	trait_data <- shap_long %>% filter(trait == !!trait)
	feature_data <- trait_data %>% filter(feature == "standage")
	observed_values <- test_data %>% pull(standage)
	
	plot <- ggplot(feature_data, aes(x = observed_values, y = shap, color = shap)) +
		geom_beeswarm(alpha = 0.5, colour = "darkgreen") +
		geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.5) +
		labs(x = "standage", y = "Influence (Shapley)") +
		theme_bw() +
		ylim(-1, 1) +
		theme(text = element_text(family = "Arial"),
					legend.position = "right",
					legend.key.width = unit(0.5, "cm"),
					legend.key.height = unit(2, "cm"),
					legend.box.background = element_rect(color = "black", linewidth = .75),
					plot.title = element_text(hjust = 0.5)) +
		ggtitle(paste(trait))
	
	plots[[trait]] <- plot
}

# Combine all plots into one figure with 2 columns and 4 rows
shap_standage_plots <- wrap_plots(plots, ncol = 2, nrow = 4)
print(shap_standage_plots)

ggsave(filename = paste0(path_out, "/plots/global_shap_standage_plots.png"),
			 plot = shap_standage_plots,
			 bg = "white",
			 width = 250,
			 height = 200,
			 units = "mm",
			 dpi = 1457)

##### Stratify data based on annual mean temperature and make predictions, using 10% quantiles, 50% quantiles, and halves: 

## ---------- Use 10% quantiles ----------

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

lower_results <- train_and_predict(stratified_data_climate$lower_data, traits, covariates, hyper_grid)
upper_results <- train_and_predict(stratified_data_climate$upper_data, traits, covariates, hyper_grid)

# Retrieve performance metrics for lower and upper quantiles
quant_performance_metrics <- bind_rows(
	lower_results$quant_performance_metrics %>% mutate(group = "Lower 10%"),
	upper_results$quant_performance_metrics %>% mutate(group = "Upper 10%")) %>%
	
	write_csv(file = paste0(path_out, "/tables/climate_10_quantiles_rf_performance_metrics.csv"))

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

# Function to create individual partial dependence plots with R-squared values and custom colours
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
	
	# Calculate residuals for jitter
	lower_preds <- predict(lower_model, lower_results$test_data)
	upper_preds <- predict(upper_model, upper_results$test_data)
	
	lower_results$test_data$yhat <- lower_preds$predictions
	upper_results$test_data$yhat <- upper_preds$predictions
	
	lower_results$test_data <- lower_results$test_data %>%
		mutate(residuals = yhat - !!sym(trait),
					 group = "Lower 10%")
	
	upper_results$test_data <- upper_results$test_data %>%
		mutate(residuals = yhat - !!sym(trait),
					 group = "Upper 10%")
	
	residual_data <- bind_rows(lower_results$test_data, upper_results$test_data)
	
	plot <- ggplot(pdp_data, aes(x = standage, y = yhat, color = group, shape = group)) +
		geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, linewidth = .7, aes(fill = group)) + 
		geom_line(linewidth = .4) +
		geom_point(data = residual_data, aes(x = standage, y = yhat + residuals), alpha = 0.4) +
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
			text = element_text(family = "Arial"),
			plot.title = element_text(face = "bold"),
			plot.subtitle = element_text(hjust = 0.5)
		) +
		facet_wrap(~group, scales = "fixed")
	
	r2_tibble <- tibble(trait = trait,
											group = c("Lower 10%", "Upper 10%"),
											r2 = c(lower_r2, upper_r2))
	
	return(list(plot = plot, r2_tibble = r2_tibble))
}

# Generate plot data for all traits
pdp_plots_climate <- list()
r2_values <- tibble(trait = character(), group = character(), r2 = numeric())
for (trait in traits) {
	result <- plot_partial_dependence_climate(lower_results, upper_results, trait)
	pdp_plots_climate[[trait]] <- result$plot
	r2_values <- bind_rows(r2_values, result$r2_tibble)
}

# Combine the pdp plots using patchwork
pdp_plots_climate <- wrap_plots(pdp_plots_climate, ncol = 2, nrow = 4) +
	plot_annotation(
		theme = theme(
			text = element_text(family = "Arial"),
			plot.title = element_text(face = "bold")
		)
	)

# Display the pdp plot
print(pdp_plots_climate) # examine figure
ggsave(filename = paste0(path_out, "/plots/climate_10_quantiles_pdp.png"),
			 plot = pdp_plots_climate,
			 bg = "white",
			 width = 200,
			 height = 250,
			 units = "mm",
			 dpi = 1457)

# Create bar chart for R-squared values
r2_plot <- ggplot(r2_values, aes(x = trait, y = r2, fill = group)) +
	geom_bar(stat = "identity", position = position_dodge(width = 0.5), colour = "black") +
	scale_fill_manual(values = c("Lower 10%" = "#0800af", "Upper 10%" = "#c82300")) +
	labs(x = NULL, y = "R²") +
	guides(fill = guide_legend(title = "Annual mean temperature quantiles:")) +
	theme_bw() +
	theme(
		text = element_text(family = "Arial"),
		plot.title = element_text(face = "bold"),
		axis.text.x = element_text(angle = 45, hjust = 1),
		legend.position = "top",
		legend.justification = "left",
		legend.background = element_rect(fill = "transparent"))

# Display the bar chart
print(r2_plot) # examine R-squared values
ggsave(filename = paste0(path_out, "/plots/climate_10_quantiles_rsquared.png"),
			 plot = r2_plot,
			 bg = "white",
			 width = 150,
			 height = 100,
			 units = "mm",
			 dpi = 1457)

## ---------- Use 50% quantiles ---------- 

# Stratify the data based on the upper and lower 50% quantiles of annual mean temperature

stratify_climate_50 <- function(data, variable) {
	lower_quantile_50 <- quantile(data[[variable]], 0.25)
	upper_quantile_50 <- quantile(data[[variable]], 0.75)
	
	lower_data_50 <- data %>% filter(.[[variable]] <= lower_quantile_50)
	upper_data_50 <- data %>% filter(.[[variable]] >= upper_quantile_50)
	
	return(list(lower_data_50 = lower_data_50, upper_data_50 = upper_data_50))
}

# Perform stratification for 50% quantiles and model training
stratified_data_climate_50 <- stratify_climate_50(data, "annual_mean_temperature")

lower_results_50 <- train_and_predict(stratified_data_climate_50$lower_data_50, traits, covariates, hyper_grid)
upper_results_50 <- train_and_predict(stratified_data_climate_50$upper_data_50, traits, covariates, hyper_grid)

# Retrieve performance metrics for lower and upper 50% quantiles
quant_performance_metrics_50 <- bind_rows(
	lower_results_50$quant_performance_metrics %>% mutate(group = "Lower 50%"),
	upper_results_50$quant_performance_metrics %>% mutate(group = "Upper 50%")
) %>%
	write_csv(file = paste0(path_out, "/tables/climate_50_quantiles_rf_performance_metrics.csv"))

# Function to create individual partial dependence plots with R-squared values and custom colours for 50% quantiles
plot_partial_dependence_climate_50 <- function(lower_results_50, upper_results_50, trait) {
	
	# Extract models
	lower_model_50 <- lower_results_50$best_models[[which(traits == trait)]]
	upper_model_50 <- upper_results_50$best_models[[which(traits == trait)]]
	
	# Generate partial dependence data
	pdp_lower_50 <- pdp::partial(lower_model_50, pred.var = "standage", train = lower_results_50$test_data)
	pdp_upper_50 <- pdp::partial(upper_model_50, pred.var = "standage", train = upper_results_50$test_data)
	
	# Add quantile identifiers
	pdp_lower_50$group <- "Lower 50%"
	pdp_upper_50$group <- "Upper 50%"
	
	# Combine partial dependence data
	pdp_data_50 <- bind_rows(pdp_lower_50, pdp_upper_50)
	
	# Extract R-squared values
	lower_r2_50 <- lower_results_50$r_squared_values[[which(traits == trait)]]
	upper_r2_50 <- upper_results_50$r_squared_values[[which(traits == trait)]]
	
	# Calculate residuals for jitter
	lower_preds_50 <- predict(lower_model_50, lower_results_50$test_data)
	upper_preds_50 <- predict(upper_model_50, upper_results_50$test_data)
	
	lower_results_50$test_data$yhat <- lower_preds_50$predictions
	upper_results_50$test_data$yhat <- upper_preds_50$predictions
	
	lower_results_50$test_data <- lower_results_50$test_data %>%
		mutate(residuals = yhat - !!sym(trait),
					 group = "Lower 50%")
	
	upper_results_50$test_data <- upper_results_50$test_data %>%
		mutate(residuals = yhat - !!sym(trait),
					 group = "Upper 50%")
	
	residual_data_50 <- bind_rows(lower_results_50$test_data, upper_results_50$test_data)
	
	plot_50 <- ggplot(pdp_data_50, aes(x = standage, y = yhat, color = group, shape = group)) +
		geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, linewidth = .7, aes(fill = group)) + 
		geom_line(linewidth = .4) +
		geom_point(data = residual_data_50, aes(x = standage, y = yhat + residuals), alpha = 0.4) +
		scale_color_manual(values = c("Lower 50%" = "#0800af", "Upper 50%" = "#c82300")) +
		scale_shape_manual(values = c(16, 18)) + 
		labs(title = trait_title(trait),
				 subtitle = paste("Lower 50% R\u00B2:", round(lower_r2_50, 2), "| Upper 50% R\u00B2:", round(upper_r2_50, 2)),
				 x = "Standage", y = "Trait value (log scaled)") + 
		theme_bw() +
		theme(
			legend.position = "none",
			legend.title = element_blank(),
			legend.text = element_text(size = 10),
			text = element_text(family = "Arial"),
			plot.title = element_text(face = "bold"),
			plot.subtitle = element_text(hjust = 0.5)
		) +
		facet_wrap(~group, scales = "fixed")
	
	r2_tibble_50 <- tibble(trait = trait,
												 group = c("Lower 50%", "Upper 50%"),
												 r2 = c(lower_r2_50, upper_r2_50))
	
	return(list(plot = plot_50, r2_tibble = r2_tibble_50))
}

# Generate plot data for all traits for 50% quantiles
pdp_plots_climate_50 <- list()
r2_values_50 <- tibble(trait = character(), group = character(), r2 = numeric())
for (trait in traits) {
	result_50 <- plot_partial_dependence_climate_50(lower_results_50, upper_results_50, trait)
	pdp_plots_climate_50[[trait]] <- result_50$plot
	r2_values_50 <- bind_rows(r2_values_50, result_50$r2_tibble)
}

# Combine the pdp plots using patchwork for 50% quantiles
pdp_plots_climate_50 <- wrap_plots(pdp_plots_climate_50, ncol = 2, nrow = 4) +
	plot_annotation(
		theme = theme(
			text = element_text(family = "Arial"),
			plot.title = element_text(face = "bold")
		)
	)

# Display the pdp plot for 50% quantiles
print(pdp_plots_climate_50) # examine figure
ggsave(filename = paste0(path_out, "/plots/climate_50_quantiles_pdp.png"),
			 plot = pdp_plots_climate_50,
			 bg = "white",
			 width = 200,
			 height = 250,
			 units = "mm",
			 dpi = 1457)

# Create bar chart for R-squared values for 50% quantiles
r2_plot_50 <- ggplot(r2_values_50, aes(x = trait, y = r2, fill = group)) +
	geom_bar(stat = "identity", position = position_dodge(width = 0.5), colour = "black") +
	scale_fill_manual(values = c("Lower 50%" = "#0800af", "Upper 50%" = "#c82300")) +
	labs(x = NULL, y = "R²") +
	guides(fill = guide_legend(title = "Annual mean temperature quantiles:")) +
	theme_bw() +
	theme(
		text = element_text(family = "Arial"),
		plot.title = element_text(face = "bold"),
		axis.text.x = element_text(angle = 45, hjust = 1),
		legend.position = "top",
		legend.justification = "left",
		legend.background = element_rect(fill = "transparent"))

# Display the bar chart for 50% quantiles
print(r2_plot_50) # examine R-squared values
ggsave(filename = paste0(path_out, "/plots/climate_50_quantiles_rsquared.png"),
			 plot = r2_plot_50,
			 bg = "white",
			 width = 150,
			 height = 100,
			 units = "mm",
			 dpi = 1457)

## ---------- Using upper and lower halves of the data ----------

# Stratify the data based on the median of annual mean temperature
stratify_climate_halves <- function(data, variable) {
	median_value <- median(data[[variable]])
	
	lower_data_halves <- data %>% filter(.[[variable]] <= median_value)
	upper_data_halves <- data %>% filter(.[[variable]] > median_value)
	
	return(list(lower_data_halves = lower_data_halves, upper_data_halves = upper_data_halves))
}

# Perform stratification for halves and model training
stratified_data_climate_halves <- stratify_climate_halves(data, "annual_mean_temperature")

lower_results_halves <- train_and_predict(stratified_data_climate_halves$lower_data_halves, traits, covariates, hyper_grid)
upper_results_halves <- train_and_predict(stratified_data_climate_halves$upper_data_halves, traits, covariates, hyper_grid)

# Retrieve performance metrics for lower and upper halves
quant_performance_metrics_halves <- bind_rows(
	lower_results_halves$quant_performance_metrics %>% mutate(group = "Lower Half"),
	upper_results_halves$quant_performance_metrics %>% mutate(group = "Upper Half")
) %>%
	write_csv(file = paste0(path_out, "/tables/climate_halves_rf_performance_metrics.csv"))

# Function to create individual partial dependence plots with R-squared values and custom colours for halves
plot_partial_dependence_climate_halves <- function(lower_results_halves, upper_results_halves, trait) {
	
	# Extract models
	lower_model_halves <- lower_results_halves$best_models[[which(traits == trait)]]
	upper_model_halves <- upper_results_halves$best_models[[which(traits == trait)]]
	
	# Generate partial dependence data
	pdp_lower_halves <- pdp::partial(lower_model_halves, pred.var = "standage", train = lower_results_halves$test_data)
	pdp_upper_halves <- pdp::partial(upper_model_halves, pred.var = "standage", train = upper_results_halves$test_data)
	
	# Add quantile identifiers
	pdp_lower_halves$group <- "Lower Half"
	pdp_upper_halves$group <- "Upper Half"
	
	# Combine partial dependence data
	pdp_data_halves <- bind_rows(pdp_lower_halves, pdp_upper_halves)
	
	# Extract R-squared values
	lower_r2_halves <- lower_results_halves$r_squared_values[[which(traits == trait)]]
	upper_r2_halves <- upper_results_halves$r_squared_values[[which(traits == trait)]]
	
	# Calculate residuals for jitter
	lower_preds_halves <- predict(lower_model_halves, lower_results_halves$test_data)
	upper_preds_halves <- predict(upper_model_halves, upper_results_halves$test_data)
	
	lower_results_halves$test_data$yhat <- lower_preds_halves$predictions
	upper_results_halves$test_data$yhat <- upper_preds_halves$predictions
	
	lower_results_halves$test_data <- lower_results_halves$test_data %>%
		mutate(residuals = yhat - !!sym(trait),
					 group = "Lower Half")
	
	upper_results_halves$test_data <- upper_results_halves$test_data %>%
		mutate(residuals = yhat - !!sym(trait),
					 group = "Upper Half")
	
	residual_data_halves <- bind_rows(lower_results_halves$test_data, upper_results_halves$test_data)
	
	plot_halves <- ggplot(pdp_data_halves, aes(x = standage, y = yhat, color = group, shape = group)) +
		geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, linewidth = .7, aes(fill = group)) + 
		geom_line(linewidth = .4) +
		geom_point(data = residual_data_halves, aes(x = standage, y = yhat + residuals), alpha = 0.4) +
		scale_color_manual(values = c("Lower Half" = "#0800af", "Upper Half" = "#c82300")) +
		scale_shape_manual(values = c(16, 18)) + 
		labs(title = trait_title(trait),
				 subtitle = paste("Lower Half R\u00B2:", round(lower_r2_halves, 2), "| Upper Half R\u00B2:", round(upper_r2_halves, 2)),
				 x = "Standage", y = "Trait value (log scaled)") + 
		theme_bw() +
		theme(
			legend.position = "none",
			legend.title = element_blank(),
			legend.text = element_text(size = 10),
			text = element_text(family = "Arial"),
			plot.title = element_text(face = "bold"),
			plot.subtitle = element_text(hjust = 0.5)
		) +
		facet_wrap(~group, scales = "fixed")
	
	r2_tibble_halves <- tibble(trait = trait,
														 group = c("Lower Half", "Upper Half"),
														 r2 = c(lower_r2_halves, upper_r2_halves))
	
	return(list(plot = plot_halves, r2_tibble = r2_tibble_halves))
}

# Generate plot data for all traits for halves
pdp_plots_climate_halves <- list()
r2_values_halves <- tibble(trait = character(), group = character(), r2 = numeric())
for (trait in traits) {
	result_halves <- plot_partial_dependence_climate_halves(lower_results_halves, upper_results_halves, trait)
	pdp_plots_climate_halves[[trait]] <- result_halves$plot
	r2_values_halves <- bind_rows(r2_values_halves, result_halves$r2_tibble)
}

# Combine the pdp plots using patchwork for halves
pdp_plots_climate_halves <- wrap_plots(pdp_plots_climate_halves, ncol = 2, nrow = 4) +
	plot_annotation(
		theme = theme(
			text = element_text(family = "Arial"),
			plot.title = element_text(face = "bold")
		)
	)

# Display the pdp plot for halves
print(pdp_plots_climate_halves) # examine figure
ggsave(filename = paste0(path_out, "/plots/climate_halves_pdp.png"),
			 plot = pdp_plots_climate_halves,
			 bg = "white",
			 width = 200,
			 height = 250,
			 units = "mm",
			 dpi = 1457)

# Create bar chart for R-squared values for halves
r2_plot_halves <- ggplot(r2_values_halves, aes(x = trait, y = r2, fill = group)) +
	geom_bar(stat = "identity", position = position_dodge(width = 0.5), colour = "black") +
	scale_fill_manual(values = c("Lower Half" = "#0800af", "Upper Half" = "#c82300")) +
	labs(x = NULL, y = "R²") +
	guides(fill = guide_legend(title = "Annual mean temperature halves:")) +
	theme_bw() +
	theme(
		text = element_text(family = "Arial"),
		plot.title = element_text(face = "bold"),
		axis.text.x = element_text(angle = 45, hjust = 1),
		legend.position = "top",
		legend.justification = "left",
		legend.background = element_rect(fill = "transparent"))

# Display the bar chart for halves
print(r2_plot_halves) # examine R-squared values
ggsave(filename = paste0(path_out, "/plots/climate_halves_rsquared.png"),
			 plot = r2_plot_halves,
			 bg = "white",
			 width = 150,
			 height = 100,
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

managed_results_0 <- train_and_predict(stratified_data_managed$managed_data_0, traits, covariates, hyper_grid)
managed_results_1 <- train_and_predict(stratified_data_managed$managed_data_1, traits, covariates, hyper_grid)

# Retrieve performance metrics for managed levels
managed_performance_metrics <- bind_rows(
	managed_results_0$quant_performance_metrics %>% mutate(group = "Managed 0"),
	managed_results_1$quant_performance_metrics %>% mutate(group = "Managed 1")) %>%
	write_csv(file = paste0(path_out, "/tables/managed_rf_performance_metrics.csv"))

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
	
	# Calculate residuals for jitter
	managed_preds_0 <- predict(managed_model_0, managed_results_0$test_data)
	managed_preds_1 <- predict(managed_model_1, managed_results_1$test_data)
	
	managed_results_0$test_data$yhat <- managed_preds_0$predictions
	managed_results_1$test_data$yhat <- managed_preds_1$predictions
	
	managed_results_0$test_data <- managed_results_0$test_data %>%
		mutate(residuals = yhat - !!sym(trait),
					 group = "Managed 0")
	
	managed_results_1$test_data <- managed_results_1$test_data %>%
		mutate(residuals = yhat - !!sym(trait),
					 group = "Managed 1")
	
	residual_data <- bind_rows(managed_results_0$test_data, managed_results_1$test_data)
	
	plot <- ggplot(pdp_data, aes(x = standage, y = yhat, color = group, shape = group)) +
		geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, linewidth = .7, aes(fill = group)) + 
		geom_line(linewidth = .4) +
		geom_jitter(data = residual_data, aes(x = standage, y = yhat + residuals), width = 0.3, height = 0, alpha = 0.4) +  # Add jittered residuals
		scale_color_manual(values = c("Managed 0" = "#0800af", "Managed 1" = "#c82300")) +
		scale_shape_manual(values = c(16, 18)) + 
		labs(title = trait_title(trait),
				 subtitle = paste("Managed 0 R\u00B2:", round(managed_r2_0, 2), "| Managed 1 R\u00B2:", round(managed_r2_1, 2)),
				 x = "Standage", y = "Trait value (log scaled)") + 
		theme_bw() +
		theme(
			legend.position = "none",
			legend.text = element_text(size = 10),
			text = element_text(family = "Arial"),
			plot.title = element_text(face = "bold"),
			plot.subtitle = element_text(hjust = 0.5)
		) +
		facet_wrap(~group, scales = "fixed")
	
	r2_tibble <- tibble(trait = trait,
											group = c("Managed 0", "Managed 1"),
											r2 = c(managed_r2_0, managed_r2_1))
	
	return(list(plot = plot, r2_tibble = r2_tibble))
}

# Generate plots for all traits
pdp_plots_managed <- list()
r2_values_managed <- tibble(trait = character(), group = character(), r2 = numeric())
for (trait in traits) {
	result <- plot_partial_dependence_managed(managed_results_0, managed_results_1, trait)
	pdp_plots_managed[[trait]] <- result$plot
	r2_values_managed <- bind_rows(r2_values_managed, result$r2_tibble)
}

# Combine the plots using patchwork
managed_plot <- wrap_plots(pdp_plots_managed, ncol = 2, nrow = 4) +
	plot_annotation(
		theme = theme(
			text = element_text(family = "Arial"),
			plot.title = element_text(face = "bold")
		)
	)

# Display the combined plot
print(managed_plot) # examine figure
ggsave(filename = paste0(path_out, "/plots/managed_pdp.png"),
			 plot = managed_plot,
			 bg = "white",
			 width = 200,
			 height = 250,
			 units = "mm",
			 dpi = 1457)

# Create bar chart for R-squared values
r2_plot_managed <- ggplot(r2_values_managed, aes(x = trait, y = r2, fill = group)) +
	geom_bar(stat = "identity", position = position_dodge(width = 0.5), colour = "black") +
	scale_fill_manual(values = c("Managed 0" = "#0800af", "Managed 1" = "#c82300")) +
	labs(x = NULL, y = "R²") +
	guides(fill = guide_legend(title = "Management status:")) +
	theme_bw() +
	theme(
		text = element_text(family = "Arial"),
		plot.title = element_text(face = "bold"),
		axis.text.x = element_text(angle = 45, hjust = 1),
		legend.position = "top",
		legend.justification = "left",
		legend.background = element_rect(fill = "transparent")
	)

# Display the bar chart
print(r2_plot_managed) # examine R-squared values
ggsave(filename = paste0(path_out, "/plots/managed_rsquared.png"),
			 plot = r2_plot_managed,
			 bg = "white",
			 width = 150,
			 height = 100,
			 units = "mm",
			 dpi = 1457)

##### Stratify data based on major biomes and make predictions 

# We'll look at the two major biomes; Temperate broadleaf forests and Temperate conifer forests (both have roughly the same number of observations)

# Stratify the data based on the biome
stratify_biome <- function(data, variable) {
	coniferous_data <- data %>% filter(.[[variable]] == "Temperate conifer forests")
	broadleaf_data <- data %>% filter(.[[variable]] == "Temperate broadleaf forests")
	
	return(list(coniferous_data = coniferous_data, broadleaf_data = broadleaf_data))
}

# Perform stratification and model training for biome variable
stratified_data_biome <- stratify_biome(data, "biome")

biome_results_coniferous <- train_and_predict(stratified_data_biome$coniferous_data, traits, covariates, hyper_grid)
biome_results_broadleaf  <- train_and_predict(stratified_data_biome$broadleaf_data, traits, covariates, hyper_grid)

# Retrieve performance metrics for managed levels
biome_performance_metrics <- bind_rows(
	biome_results_coniferous$quant_performance_metrics %>% mutate(group = "coniferous"),
	biome_results_broadleaf$quant_performance_metrics %>% mutate(group = "broadleaf")) %>%
	write_csv(file = paste0(path_out, "/tables/biome_rf_performance_metrics.csv"))

# Function to create partial dependence plots for biome variable with R-squared values
plot_partial_dependence_biome <- function(biome_results_coniferous, biome_results_broadleaf, trait) {
	
	biome_model_coniferous <- biome_results_coniferous$best_models[[which(traits == trait)]]
	biome_model_broadleaf  <- biome_results_broadleaf$best_models[[which(traits == trait)]]
	
	pdp_biome_coniferous <- pdp::partial(biome_model_coniferous, pred.var = "standage", train = biome_results_coniferous$test_data)
	pdp_biome_broadleaf  <- pdp::partial(biome_model_broadleaf, pred.var = "standage", train = biome_results_broadleaf$test_data)
	
	pdp_biome_coniferous$group <- "coniferous"
	pdp_biome_broadleaf$group <- "broadleaf"
	
	pdp_data <- bind_rows(pdp_biome_coniferous, pdp_biome_broadleaf)
	
	biome_r2_coniferous <- biome_results_coniferous$r_squared_values[[which(traits == trait)]]
	biome_r2_broadleaf  <- biome_results_broadleaf$r_squared_values[[which(traits == trait)]]
	
	# Calculate residuals for jitter
	biome_preds_coniferous <- predict(biome_model_coniferous, biome_results_coniferous$test_data)
	biome_preds_broadleaf  <- predict(biome_model_broadleaf, biome_results_broadleaf$test_data)
	
	biome_results_coniferous$test_data$yhat <- biome_preds_coniferous$predictions
	biome_results_broadleaf$test_data$yhat <- biome_preds_broadleaf$predictions
	
	biome_results_coniferous$test_data <- biome_results_coniferous$test_data %>%
		mutate(residuals = yhat - !!sym(trait),
					 group = "coniferous")
	
	biome_results_broadleaf$test_data <- biome_results_broadleaf$test_data %>%
		mutate(residuals = yhat - !!sym(trait),
					 group = "broadleaf")
	
	residual_data <- bind_rows(biome_results_coniferous$test_data, biome_results_broadleaf$test_data)
	
	plot <- ggplot(pdp_data, aes(x = standage, y = yhat, color = group, shape = group)) +
		geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, linewidth = .7, aes(fill = group)) + 
		geom_line(linewidth = .4) +
		geom_jitter(data = residual_data, aes(x = standage, y = yhat + residuals), width = 0.3, height = 0, alpha = 0.4) + 
		scale_color_manual(values = c("coniferous" = "#616536", "broadleaf" = "#954535")) +
		scale_shape_manual(values = c(16, 18)) + 
		labs(title = trait_title(trait),
				 subtitle = paste("Coniferous R\u00B2:", round(biome_r2_coniferous, 2), "| Broadleaf R\u00B2:", round(biome_r2_broadleaf, 2)),
				 x = "Standage", y = "Trait value (log scaled)") + 
		theme_bw() +
		theme(
			legend.position = "none",
			legend.text = element_text(size = 10),
			text = element_text(family = "Arial"),
			plot.title = element_text(face = "bold"),
			plot.subtitle = element_text(hjust = 0.5)
		) +
		facet_wrap(~group, scales = "fixed")
	
	r2_tibble <- tibble(trait = trait,
											group = c("coniferous", "broadleaf"),
											r2 = c(biome_r2_coniferous, biome_r2_broadleaf))
	
	return(list(plot = plot, r2_tibble = r2_tibble))
}

# Generate plots for all traits
pdp_plots_biome <- list()
r2_values_biome <- tibble(trait = character(), group = character(), r2 = numeric())
for (trait in traits) {
	result <- plot_partial_dependence_biome(biome_results_coniferous, biome_results_broadleaf, trait)
	pdp_plots_biome[[trait]] <- result$plot
	r2_values_biome <- bind_rows(r2_values_biome, result$r2_tibble)
}

# Combine the plots using patchwork
biome_plot <- wrap_plots(pdp_plots_biome, ncol = 2, nrow = 4) +
	plot_annotation(
		theme = theme(
			text = element_text(family = "Arial"),
			plot.title = element_text(face = "bold")
		)
	)

# Display the combined plot
print(biome_plot) # examine figure
ggsave(filename = paste0(path_out, "/plots/biomes_pdp.png"),
			 plot = biome_plot,
			 bg = "white",
			 width = 200,
			 height = 250,
			 units = "mm",
			 dpi = 1457)

# Create bar chart for R-squared values
r2_plot_biome <- ggplot(r2_values_biome, aes(x = trait, y = r2, fill = group)) +
	geom_bar(stat = "identity", position = position_dodge(width = 0.5), colour = "black") +
	scale_fill_manual(values = c("coniferous" = "#616536", "broadleaf" = "#954535")) +
	labs(x = NULL, y = "R²") +
	guides(fill = guide_legend(title = "Biome:")) +
	theme_bw() +
	theme(
		text = element_text(family = "Arial"),
		plot.title = element_text(face = "bold"),
		axis.text.x = element_text(angle = 45, hjust = 1),
		legend.position = "top",
		legend.justification = "left",
		legend.background = element_rect(fill = "transparent")
	)

# Display the bar chart
print(r2_plot_biome) # examine R-squared values
ggsave(filename = paste0(path_out, "/plots/biomes_rsquared.png"),
			 plot = r2_plot_biome,
			 bg = "white",
			 width = 150,
			 height = 100,
			 units = "mm",
			 dpi = 1457)

##########  End of analysis ##########

# End the parallel backend
stopCluster(cl)
