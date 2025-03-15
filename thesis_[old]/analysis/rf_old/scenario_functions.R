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

## ----- Shapley Calculations -----

# Prediction function for fastshap
predict_fn <- function(object, newdata) {
	predict(object, data = newdata)$prediction
}

# Function to create shapley plots for a single trait
plot_shapley_for_trait <- function(trait_id, hide_x_axis_labels = FALSE) {
	
	feature_importance %>%
		filter(trait == trait_id) %>%
		{
			shap_range <- shap_max - shap_min
			y_lim_min <- min(.$shap) - 0.1 * shap_range
			y_lim_max <- max(.$shap) + 0.1 * shap_range
			
			# Plot shapley bee swarms
			p1 <- ggplot(., aes(x = reorder(feature, importance), y = shap, color = shap)) +
				geom_quasirandom(alpha = 0.5) +
				scale_color_viridis_c(option = "viridis",guide = "none") +
				scale_x_discrete(labels = feature_labels) +
				scale_y_continuous(n.breaks = 5, limits = c(y_lim_min, y_lim_max)) +
				coord_flip() +
				labs(x = NULL, y = "Shapley Value", title = glue("{distinct(., trait_name)$trait_name}")) +
				theme_bw(base_line_size = .3, base_rect_size = .5) +
				theme(text = element_text(family = "sans", size = 6),
							axis.title.x = element_text(family = "sans", size = 5),
							plot.title = element_text(face = "bold", size = 6),
							strip.text = element_text(face = "bold"),
							plot.margin = margin(t=7.5, b=5, r=1, l=5)) +
				theme(axis.title.x = if(hide_x_axis_labels) element_blank() else element_text())
			
			# Plot absolute shapley importance 
			p2 <- ggplot(distinct(., feature, importance, rsq), aes(x = reorder(feature, importance), y = importance)) +
				geom_bar(stat = "identity", fill = "#1B7E74", color = "black", alpha = 0.7, linewidth = 0.2) +
				scale_y_continuous(n.breaks = 3) +
				labs(x = NULL, y = "Feat. Importance", title = glue("R² = {distinct(., rsq)$rsq}")) +
				coord_flip() +
				theme_bw(base_line_size = .3, base_rect_size = .5) +
				theme(text = element_text(family = "sans", size = 6),
							plot.title = element_text(size = 6),
							axis.text.y = element_blank(),
							axis.ticks.y = element_blank(),
							axis.title.x = element_text(family = "sans", size = 5),
							plot.margin = margin(t = 15, b = 2.5, r = 5, l = 1)) +
				theme(axis.title.x = if(hide_x_axis_labels) element_blank() else element_text())
			
			# Set plot structure 
			plot_grid(p1, p2, 
								ncol = 2, 
								nrow = 1, 
								rel_widths = c(.7, .3),
								align = "h", 
								axis = c("l"), 
								greedy = TRUE)
		}
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
	
	stratified_data <- stratify(data, variable, quantiles_25)
	
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

# Function for pdp plotting
create_pdp_plot <- function(data, levels, colors, labels, x_lab, y_lab, show_y_strip_labels) {
	data %>%
		mutate(group = factor(group, levels = levels)) %>%
		ggplot(aes(x = standage, y = yhat, color = group, shape = group)) +
		geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, linewidth = .7) +
		scale_color_manual(values = setNames(colors, labels)) +
		scale_fill_manual(values = setNames(colors, labels)) +
		labs(x = x_lab, y = y_lab, color = "Quantile", shape = "Quantile", fill = "Quantile") +
		theme_bw(base_line_size = .3, base_rect_size = .5) +
		theme(
			legend.position = "none",
			text = element_text(family = "sans", size = 8),
			strip.background.x = element_rect(fill = "white", color = "black", linewidth = .75),
			strip.background.y = if (show_y_strip_labels) element_rect(fill = "white", color = "black", linewidth = .75) else element_blank(),
			strip.text.y = if (show_y_strip_labels) element_text() else element_blank(),
			strip.text.x = element_text(face = "bold")) +
		facet_grid(trait ~ group, scales = "free", labeller = labeller(trait = trait_labels))
}


# Calculate standard error
calculate_se <- function(x) {
	n <- length(x)
	sd(x) / sqrt(n)
}

# Summarise successional PDP data
get_pdp_delta <- function(pdp_data) {
	delta <- pdp_data %>%
		group_by(trait, group) %>%
		summarize(
			intercept = first(intercept),  
			se_yhat = sd(yhat_mean) / sqrt(n()), 
			delta = last(yhat_mean) - first(yhat_mean),
			r_squared = mean(rsq)  
		) %>%
		ungroup()
	return(delta)
}



