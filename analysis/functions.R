########################################################################
########## succession_traits: Scenario Modeling - Source Code ##########
########################################################################

## This script contains additional functions called by the main script

## ----- Composite function for PCA -----

# Function to filter scores by varimax rotation
filter_scores <- function(scores, loadings, threshold = threshold) {
	
	strong_loadings <- abs(loadings) >= threshold
	filtered_scores <- matrix(NA, nrow = nrow(scores), ncol = ncol(scores))
	for (i in 1:ncol(loadings)) {
		
		filtered_scores[, i] <- ifelse(strong_loadings[, i], scores[, i], NA)
		
	}
	rownames(filtered_scores) <- rownames(scores)
	colnames(filtered_scores) <- colnames(scores)
	return(filtered_scores)
}

# Function to create contribution plots with abbreviated x-axis labels
create_contrib_plot <- function(var_contrib_data, title, fill_color, abbreviations) {
	var_contrib_data %>%
		as_tibble() %>%
		mutate(name = abbreviations[match(name, names(abbreviations))]) %>%
		ggplot(aes(x = reorder(name, -contrib), y = contrib)) +
		geom_bar(stat = "identity", color = "black", fill = fill_color, alpha = 0.8) +
		labs(y = "Contributions (%)", x = NULL, title = title) +
		theme_bw() +
		theme(text = element_text(family = "sans", size = 8),
					axis.text.x = element_text(angle = 45, hjust = 1))
}

# Function to get varimax rotation contributions for plotting
varimax_contrib <- function(loadings) {
	loadings^2 / colSums(loadings^2) * 100
}

## ------ Random Forest Modelling -----

# Define function to perform cross-validation and find the best hyperparameters
tune_rf_model <- function(trait, data, covariates, hyper_grid, num_threads = 1) {
	
	# Define model formula
	formula <- as.formula(paste(trait, "~", paste(covariates, collapse = " + ")))
	
	# Perform tuning
	prediction_error <- foreach(
		num.trees = hyper_grid$num.trees,
		mtry = hyper_grid$mtry,
		min.node.size = hyper_grid$min.node.size,
		.combine = 'c', 
		.packages = "ranger"
		
	) %dopar% {
		
		# Fit models
		mod <- ranger::ranger(
			formula = formula,
			data = data,
			dependent.variable.name = trait,
			num.trees = num.trees,
			mtry = mtry,
			min.node.size = min.node.size,
			num.threads = num_threads)
		
		# Return prediction error
		return(mod$prediction.error)
	}
	
	# Plot prediction errors
	error_plot <- hyper_grid %>% 
		mutate(prediction_error = prediction_error) %>%
		ggplot(aes(x = mtry, y = as.factor(min.node.size), fill = prediction_error)) + 
		facet_wrap(~ num.trees) + 
		geom_tile() + 
		scale_y_discrete(breaks = c(1, 10, 20)) + 
		scale_fill_viridis_c() + 
		ylab("min.node.size") + 
		ggtitle(trait) + 
		theme(text = element_text(family = "sans", size = 6),
					legend.key.height = unit(3, "mm"),
					legend.key.width = unit(3, "mm"),
					plot.title = element_text(face = "bold", size = 6))
	
	# Get best set 
	best_hyperparameters <- hyper_grid %>% 
		mutate(trait = trait) %>%
		arrange(prediction_error) %>% 
		slice(1)
	
	return(list(error_plot = error_plot, best_hyperparameters = best_hyperparameters))
}

# Define function to fit a random forest with best hyperparameters per trait
fit_rf_model <- function(trait, data, covariates, hyper_parameters, num_threads = 1) {
	
	# Define model formula for each trait
	formula <- as.formula(paste(trait, "~", paste(covariates, collapse = " + ")))
	
	# Get trait-specific hyperparameters
	hyper_grid <- hyper_parameters[hyper_parameters$trait == trait, ]
	
	# Fit model
	trait_mod <- ranger::ranger(
		formula = formula,
		data = data,
		num.trees = hyper_grid$num.trees[1],
		mtry = hyper_grid$mtry[1],
		min.node.size = hyper_grid$min.node.size[1],
		num.threads = num_threads)
	
	# Get performance metrics 
	performance <- data.frame(
		trait = trait,
		mtry = trait_mod$mtry,
		num_trees = trait_mod$num.trees,
		min_node_size = trait_mod$min.node.size,
		rsq = trait_mod$r.squared,
		pred_error = trait_mod$prediction.error)
	
	# Return model and performance metrics 
	return(list(trait_mod = trait_mod, performance = performance))
}

## ----- Shapley calculations and bootstrapping -----

# Prediction function for fastshap
predict_fn <- function(object, newdata) {
	predict(object, data = newdata)$prediction
}

# Function to train models and compute SHAP values per bootstrap iteration
bootstrap_shap <- function(data, traits, covariates, hyper_grid, num_threads = num_cores, boot_id) {
	
	# Sample with replacement
	boot_split <- initial_split(data, prop = 0.8, strata = NULL)
	train_data <- training(boot_split)
	test_data <- testing(boot_split)
	
	# Fit models in parallel using foreach
	best_models <- foreach(trait = traits, .packages = c('ranger', 'dplyr')) %dopar% {
		fit_rf_model(trait, train_data, covariates, hyper_grid, num_threads)
	}
	
	# Extract models
	best_models <- map(best_models, "trait_mod")
	
	# Compute SHAP values for each trait
	shap_values_list <- foreach(i = seq_along(best_models), .combine = 'rbind', .packages = c('fastshap', 'dplyr', 'tidyverse')) %dopar% {
		
		model <- best_models[[i]]
		trait <- traits[i]
		
		shap_values <- fastshap::explain(model, 
																		 X = test_data %>% select(all_of(covariates)) %>% as.matrix(),
																		 pred_wrapper = predict_fn, 
																		 nsim = 10, 
																		 parallel = TRUE)  
		
		shap_values <- as.data.frame(shap_values)
		shap_values$trait <- trait
		shap_values$bootstrap_id <- boot_id  # Track bootstrap iteration
		
		# Reshape for external processing
		shap_values %>%
			pivot_longer(-c(trait, bootstrap_id), names_to = "variable", values_to = "shap_value")
	}
	
	return(shap_values_list)
}











## ----- Partial Dependence -----

# Define data stratification
stratify <- function(data, variable, quantiles) {
	
	lower_quantile <- quantile(data[[variable]], quantiles[1])
	upper_quantile <- quantile(data[[variable]], quantiles[2])
	
	lower_data <- data %>% filter(.[[variable]] <= lower_quantile)
	upper_data <- data %>% filter(.[[variable]] >= upper_quantile)
	
	return(list(lower = lower_data, upper = upper_data))
}

# Function to stratify, fit models, and extract performance metrics
fit_models_on_strata <- function(data, variable, traits, covariates, hyper_grid, num_threads, labels) {
	
	stratified_data <- stratify(data, variable, c(0.25, 0.75))
	
	lower_models <- map(traits, ~ fit_rf_model(.x, stratified_data$lower, covariates, hyper_grid, num_threads))
	names(lower_models) <- traits
	
	upper_models <- map(traits, ~ fit_rf_model(.x, stratified_data$upper, covariates, hyper_grid, num_threads))
	names(upper_models) <- traits
	
	lower_performance <- map(lower_models, "performance") %>% bind_rows() %>% mutate(group = labels[1])
	upper_performance <- map(upper_models, "performance") %>% bind_rows() %>% mutate(group = labels[2])
	
	performance_metrics <- bind_rows(lower_performance, upper_performance)
	
	return(list(models = list(lower = lower_models, upper = upper_models), 
							performance = performance_metrics,
							stratified_data = stratified_data))
}

# Function to calculate partial dependence and confidence intervals
calculate_partial_dependence <- function(model, data, feature, conf.level = 0.95) {
	
	partial_results <- partial(model, pred.var = feature, train = data, parallel = FALSE)
	
	# Create a template for predictions
	pred_data <- data
	
	# Calculate predictions for each grid point
	predictions <- foreach(i = seq_len(nrow(partial_results)), .combine = rbind, .packages = "ranger") %dopar% {
		pred_data[[feature]] <- partial_results[i, feature]
		predict(model, data = pred_data)$predictions
	}
	
	# Calculate mean and confidence intervals
	partial_results$yhat_mean <- rowMeans(predictions)
	partial_results$yhat_lower <- apply(predictions, 1, function(x) quantile(x, probs = (1 - conf.level) / 2))
	partial_results$yhat_upper <- apply(predictions, 1, function(x) quantile(x, probs = 1 - (1 - conf.level) / 2))
	
	# Calculate intercepts 
	partial_results$intercept <- partial_results$yhat_mean[1]
	
	return(partial_results)
}

# Calculate partial dependence for each trait and scenario with group labels
calculate_pdp_for_scenario <- function(scenario_data, traits, feature, labels) {
	
	results <- foreach(trait = traits, .combine = rbind, .packages = c('ranger', 'pdp', 'dplyr', 'foreach'), 
										 .export = c('calculate_partial_dependence')) %dopar% {
										 	
										 	model_lower <- scenario_data$models$lower[[trait]][["trait_mod"]]
										 	model_upper <- scenario_data$models$upper[[trait]][["trait_mod"]]
										 	
										 	data_lower <- scenario_data$stratified_data$lower
										 	data_upper <- scenario_data$stratified_data$upper
										 	
										 	# Diagnostic prints
										 	print(paste("Trait:", trait))
										 	print("Model lower class:")
										 	print(class(model_lower))
										 	print("Model upper class:")
										 	print(class(model_upper))
										 	print("Data lower class:")
										 	print(class(data_lower))
										 	print("Data upper class:")
										 	print(class(data_upper))
										 	
										 	pdp_lower <- calculate_partial_dependence(model_lower, data_lower, feature)
										 	pdp_upper <- calculate_partial_dependence(model_upper, data_upper, feature)
										 	
										 	pdp_lower$group <- labels[1]
										 	pdp_upper$group <- labels[2]
										 	pdp_lower$trait <- trait
										 	pdp_upper$trait <- trait
										 	
										 	bind_rows(pdp_lower, pdp_upper)
										 }
	return(results)
}


# Function to recode groups into "high" and "low"
recode_group <- function(group) {
	case_when(
		str_detect(group, "Upper 25%") ~ "high",
		str_detect(group, "Lower 25%") ~ "low"
	)
}


