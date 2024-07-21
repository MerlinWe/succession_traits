########################################################################
############ succession_traits: Grid Modeling - Source Code ############
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
	loadings^2 / rowSums(loadings^2) * 100
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

# General function to create partial dependence data for each trait based on quantiles
get_partial_dependence <- function(best_models, train_data, traits, pred.vars, quantiles, quantile_names = NULL) {
	
	# Get quantile values for covariate
	quantile_values <- quantile(train_data[[pred.vars[2]]], quantiles)
	
	if (is.null(quantile_names)) {
		quantile_names <- paste0("Quantile ", seq_along(quantiles))
	}
	
	# Print the threshold values for the quantiles
	for (i in seq_along(quantiles)) {
		cat(quantile_names[i], "threshold value:", quantile_values[i], "\n")
	}
	
	# Separate data for each quantile
	data_list <- lapply(seq_along(quantile_values), function(i) {
		q_value <- quantile_values[i]
		if (i == 1) {
			train_data %>% filter(!!sym(pred.vars[2]) <= q_value)
		} else {
			train_data %>% filter(!!sym(pred.vars[2]) >= q_value)
		}
	})
	
	# Perform partial dependence calculation in parallel
	pdp_data_combined <- foreach(trait = traits, .combine = 'rbind', .packages = c('pdp', 'dplyr')) %dopar% {
		model <- best_models[[which(traits == trait)]]
		
		# Combine partial dependence data for each quantile
		pdp_data_list <- lapply(seq_along(quantile_values), function(i) {
			q_value <- quantile_values[i]
			data_subset <- data_list[[i]]
			
			# Build grid for partial dependence
			pred_grid <- expand.grid(
				standage = unique(data_subset$standage))
			pred_grid[[pred.vars[2]]] <- q_value
			
			# Calculate partial dependence
			pdp_data <- pdp::partial(
				object = model,
				pred.var = pred.vars,
				pred.grid = pred_grid,
				train = data_subset,
				plot = FALSE, 
				parallel = FALSE)
			
			# Label the data with the quantile name
			pdp_data %>% mutate(group = quantile_names[i])
		})
		
		# Calculate uncertainty as 95% CI assuming a normal distribution
		pdp_data_combined <- bind_rows(pdp_data_list) %>%
			mutate(trait = trait,
						 lower = yhat - 1.96 * sd(yhat),
						 upper = yhat + 1.96 * sd(yhat))
		
		return(pdp_data_combined)
	}
	
	return(pdp_data_combined)
}

# Define a function to get partial dependence for different variables
get_pdp <- function(variable, names) {
	get_partial_dependence(
		best_models = best_models,
		train_data = train_data,
		traits = traits,
		pred.vars = c("standage", variable),
		quantiles = quantiles_25,
		quantile_names = names)
}

# Get partial dependence for each covariate of
pdp_results <- lapply(names(variables), function(var) {
	get_pdp(var, variables[[var]])
})


# Define a function to create pdp plots
create_pdp_plot <- function(data, levels, colors, labels, x_lab, y_lab) {
	data %>%
		mutate(group = factor(group, levels = levels)) %>%
		ggplot(aes(x = standage, y = yhat, color = group, shape = group)) +
		geom_ribbon(aes(ymin = lower, ymax = upper, fill = group), alpha = 0.2, linewidth = .3) +
		geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, linewidth = .5) +
		scale_color_manual(values = setNames(colors, labels)) +
		scale_fill_manual(values = setNames(colors, labels)) +
		labs(x = x_lab, y = y_lab, color = "Quantile", shape = "Quantile", fill = "Quantile") +
		theme_bw(base_line_size = .3, base_rect_size = .5) +
		theme(
			legend.position = "none",
			text = element_text(family = "sans", size = 8),
			strip.background.y = element_blank(),
			strip.text.y = element_blank(),
			strip.background = element_rect(fill = "white", color = "black", linewidth = .75)) +
		facet_grid(trait ~ group, scales = "free", labeller = labeller(trait = trait_labels))
}

process_pdp_data <- function(pdp_data) {
	slopes_intercepts <- pdp_data %>%
		group_by(trait, group) %>%
		do({
			fit <- lm(yhat ~ standage, data = .)
			intercept <- coef(fit)[1]
			slope <- coef(fit)[2]
			max_standage <- max(.$standage)
			y_max_standage <- intercept + slope * max_standage
			delta <- y_max_standage - intercept
			tibble(intercept = intercept, slope = slope, sd_yhat = sd(.$yhat), delta = delta)
		}) %>%
		ungroup()
	return(slopes_intercepts)
}





