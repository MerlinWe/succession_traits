########################################################################
########## succession_traits: helper functions source code ##########
########################################################################

## This script contains helper functions called by the main script

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

# ---------- Helpers for biome-stratification ----------

# Derive a single biome_group factor from dummies
make_biome_group <- function(df, biome_dummy_cols) {
	# Pick the 1-coded biome; if multiple or none, mark "other"
	bm <- df[, biome_dummy_cols, drop = FALSE]
	idx <- max.col(bm, ties.method = "first")  # fastest, but assumes at least one 1
	has_one <- rowSums(bm == 1) == 1
	grp <- rep("other", nrow(df))
	grp[has_one] <- sub("^biome_", "", biome_dummy_cols[idx[has_one]])
	
	# Map to 3 broad groups
	map3 <- c(
		"boreal_forests_or_taiga"   = "forest",
		"temperate_broadleaf_forests" = "forest",
		"temperate_conifer_forests" = "forest",
		"temperate_grasslands"      = "grassland",
		"flooded_grasslands"        = "grassland",
		"xeric_shrublands"          = "xeric_med",
		"mediterranean_woodlands"   = "xeric_med",
		"tundra"                    = "other"     # exclude?
	)
	grp3 <- dplyr::recode(grp, !!!map3, .default = "other")
	factor(grp3, levels = c("forest","grassland","xeric_med","other"))
}

# --------- RF modelling helpers ---------
 
split_train_test <- function(df, ntest = NULL, test_frac = 0.2, seed = 42) {
	set.seed(seed)
	df <- dplyr::mutate(df, .id = dplyr::row_number())
	if (is.null(ntest)) {
		ntest <- max(1L, floor(nrow(df) * test_frac))
	} else {
		ntest <- min(ntest, nrow(df))
	}
	test_ids <- sample(df$.id, size = ntest)
	list(
		train = dplyr::filter(df, !.id %in% test_ids),
		test  = dplyr::filter(df,  .id %in% test_ids)
	)
}

# RF tuner: iterate the FULL Cartesian grid; use OOB error
tune_rf_model <- function(trait, data, covariates, hyper_grid, num_threads = 1) {
	formula <- as.formula(paste(trait, "~", paste(covariates, collapse = " + ")))
	hg <- hyper_grid
	
	# parallel over rows of grid; avoid overkill -> num_threads = 1 inside ranger
	pred_err <- foreach::foreach(i = seq_len(nrow(hg)), .combine = "c", .packages = "ranger") %dopar% {
		row <- hg[i, ]
		mod <- ranger::ranger(
			formula = formula, data = data,
			num.trees = row$num.trees, mtry = row$mtry, min.node.size = row$min.node.size,
			num.threads = 1L, importance = "none",
			respect.unordered.factors = "order", keep.inbag = FALSE, seed = 15
		)
		mod$prediction.error  # OOB error
	}
	
	error_plot <- dplyr::bind_cols(hg, prediction_error = pred_err) %>%
		ggplot2::ggplot(ggplot2::aes(x = mtry, y = as.factor(min.node.size), fill = prediction_error)) +
		ggplot2::facet_wrap(~ num.trees) +
		ggplot2::geom_tile() +
		ggplot2::scale_y_discrete(breaks = c(1,10,20)) +
		ggplot2::scale_fill_viridis_c() +
		ggplot2::ylab("min.node.size") +
		ggplot2::ggtitle(trait) +
		ggplot2::theme(
			text = ggplot2::element_text(family = "sans", size = 6),
			legend.key.height = grid::unit(3, "mm"),
			legend.key.width  = grid::unit(3, "mm"),
			plot.title = ggplot2::element_text(face = "bold", size = 6)
		)
	
	best <- dplyr::bind_cols(hg, prediction_error = pred_err) %>%
		dplyr::mutate(trait = trait) %>%
		dplyr::arrange(prediction_error) %>%
		dplyr::slice(1)
	
	list(error_plot = error_plot, best_hyperparameters = best)
}

# Ranger random forest fit function
fit_rf_model <- function(trait, df_train, covariates, hyper_parameters,
												 num_threads = 1, case_weights_col = NULL) {
	stopifnot(trait %in% names(df_train))
	formula <- as.formula(paste(trait, "~", paste(covariates, collapse = " + ")))
	row <- hyper_parameters[hyper_parameters$trait == trait, , drop = FALSE]
	if (nrow(row) != 1) stop("Hyper-grid missing/duplicated row for trait: ", trait)
	
	cw <- if (!is.null(case_weights_col) && case_weights_col %in% names(df_train))
		df_train[[case_weights_col]] else NULL
	
	mod <- ranger::ranger(
		formula = formula, data = df_train,
		num.trees = row$num_trees[1], mtry = row$mtry[1], min.node.size = row$min_node_size[1],
		num.threads = num_threads, case.weights = cw,
		importance = "none", respect.unordered.factors = "order", keep.inbag = FALSE, seed = 15
	)
	
	perf <- tibble::tibble(
		trait = trait,
		mtry = mod$mtry,
		num_trees = mod$num.trees,
		min_node_size = mod$min.node.size,
		rsq = mod$r.squared,
		pred_error = mod$prediction.error
	)
	list(trait_mod = mod, performance = perf)
}

## ----- Shapley -----

predict_fn <- function(object, newdata) {
	p <- predict(object, data = as.data.frame(newdata))
	if (!is.null(p$predictions)) return(p$predictions)
	if (!is.null(p$prediction))  return(p$prediction)
	stop("Unknown predict() return structure from ranger.")
}

# Robust symmetric scaling for colour
symmetric_scale <- function(x) {
	x <- as.numeric(x)
	if (!any(is.finite(x))) return(rep(0, length(x)))
	z <- as.numeric(scale(x))
	z[!is.finite(z)] <- 0
	pmax(pmin(z, 3), -3) / 3
}

plot_beeswarm_by_leaf <- function(df, trait_label, max_points = 5000) {
	
	df <- df %>%
		mutate(leaf_type = recode(leaf_type,
															"broadleaf" = "Broadleaf",
															"coniferous" = "Coniferous"
		))
	
	# ---- Define pretty labels for predictor variables ----
	var_labels <- c(
		"standage"   = "Stand age",
		"temp_pc"    = "Temperature PC",
		"rain_pc"    = "Precipitation PC",
		"soil_pc"    = "Soil PC",
		"elevation"  = "Elevation",
		"soil_ph"    = "Soil pH"
	)
	
	# ---- Filter, rescale, and sample ----
	df_sub <- df %>%
		filter(trait == trait_label) %>%
		group_by(variable) %>%
		mutate(value_col = scales::rescale(feature_value, to = c(0, 1))) %>%
		ungroup() %>%
		{ if (nrow(.) > max_points) slice_sample(., n = max_points) else . }
	
	# ---- Beeswarm plot ----
	ggplot(df_sub,
				 aes(x = reorder(variable, shap_value, FUN = median),
				 		y = shap_value,
				 		color = value_col)) +
		ggbeeswarm::geom_quasirandom(alpha = 0.35, size = 0.9) +
		facet_wrap(~leaf_type, ncol = 1, scales = "free_y") +
		coord_flip() +
		scale_x_discrete(labels = var_labels) +
		scale_color_viridis_c(name = "Feature value") +
		labs(y = "SHAP value", x = NULL) +
		theme_bw(base_size = 11) +
		theme(
			panel.grid.minor = element_blank(),
			panel.grid.major.y = element_blank(),
			strip.background = element_rect(fill = "grey95"),
			strip.text = element_text(face = "bold", size = 11),
			axis.text.y = element_text(size = 9, hjust = 1),
			axis.text.x = element_text(size = 9),
			axis.title.y = element_text(size = 11),
			legend.position = "right",
			legend.key.height = unit(3.5, "mm")
		)
}

plot_shap_dependence <- function(df, trait_label, var) {
	
	df <- df %>%
		mutate(leaf_type = recode(leaf_type,
															"broadleaf" = "Broadleaf",
															"coniferous" = "Coniferous"
		))
	
	# --- Human-readable axis labels ---
	var_labels <- c(
		"standage"   = "Stand age",
		"temp_pc"    = "Temperature PC",
		"rain_pc"    = "Precipitation PC",
		"soil_pc"    = "Soil PC",
		"elevation"  = "Elevation",
		"soil_ph"    = "Soil pH"
	)
	
	x_label <- ifelse(var %in% names(var_labels),
										var_labels[[var]],
										paste0(var, " (feature value)"))
	
	ggplot(
		df %>% filter(trait == trait_label, variable == var),
		aes(x = feature_value, y = shap_value,
				color = leaf_type, shape = leaf_type)
	) +
		geom_point(alpha = 0.25, size = 1.4) +
		scale_color_manual(
			name = "Leaf type",
			values = c("Broadleaf" = "#228B22", "Coniferous" = "#d95f02")
		) +
		scale_shape_manual(
			name = "Leaf type",
			values = c("Broadleaf" = 16, "Coniferous" = 17)
		) +
		labs(
			x = x_label,
			y = "SHAP value"
		) +
		theme_bw(base_size = 11) +
		theme(
			legend.position = c(0.82, 0.1),
			legend.title = element_blank(),
			legend.background = element_rect(fill = "white", color = "grey70"),
			panel.grid.minor = element_blank(),
			panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
			axis.text = element_text(size = 9),
			axis.title = element_text(size = 10),
			plot.margin = margin(5, 5, 5, 5),
			plot.title = element_blank()
		)
}

# Identify top 3 environmental drivers per leaf type (excluding standage)
get_top_env_vars <- function(shap_df, trait_label, top_n = 3) {
	shap_df %>%
		filter(trait == trait_label, variable != "standage") %>%
		group_by(leaf_type, variable) %>%
		summarise(mean_abs = mean(abs(shap_value), na.rm = TRUE), .groups = "drop") %>%
		group_by(leaf_type) %>%
		slice_max(order_by = mean_abs, n = top_n, with_ties = FALSE) %>%
		pull(variable) %>%
		unique()
}





























## ----- Partial Dependence; Scenarios; Biomes -----

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


recode_group <- function(g) dplyr::case_when(
	stringr::str_detect(g, "Upper 25%|Upper 10%") ~ "high",
	stringr::str_detect(g, "Lower 25%|Lower 10%") ~ "low",
	TRUE ~ NA_character_
)

# Quantile stratification WITHIN the subset passed (e.g., within a biome)
stratify_within <- function(data, variable, probs = c(0.25, 0.75)) {
	qs <- quantile(data[[variable]], probs, na.rm = TRUE)
	list(
		lower = dplyr::filter(data, .data[[variable]] <= qs[1]),
		upper = dplyr::filter(data, .data[[variable]] >= qs[2])
	)
}

fit_models_on_strata_biome <- function(data, variable, traits, covariates, hyper_grid, num_threads, labels, probs) {
	S <- stratify_within(data, variable, probs)
	lower_models <- purrr::map(traits, ~ fit_rf_model(.x, S$lower, covariates, hyper_grid, num_threads)); names(lower_models) <- traits
	upper_models <- purrr::map(traits, ~ fit_rf_model(.x, S$upper, covariates, hyper_grid, num_threads)); names(upper_models) <- traits
	
	lower_perf <- purrr::map(lower_models, "performance") %>% dplyr::bind_rows() %>% dplyr::mutate(group = labels[1])
	upper_perf <- purrr::map(upper_models, "performance") %>% dplyr::bind_rows() %>% dplyr::mutate(group = labels[2])
	list(models = list(lower = lower_models, upper = upper_models),
			 performance = dplyr::bind_rows(lower_perf, upper_perf),
			 stratified_data = S)
}

# Wrapper: run a full PDP bootstrap inside ONE biome for all env vars
pdp_biome_once <- function(boot_B, traits, covariates, hyper_grid, num_threads, probs) {
	# Adjust labels for 10/90 sensitivity if needed
	qtxt <- if (identical(probs, c(0.10,0.90))) "10%" else "25%"
	qlabels <- list(
		temp_pc = c(paste0("Cold Temperatures (Lower ", qtxt, ")"),
								paste0("Warm Temperatures (Upper ", qtxt, ")")),
		soil_pc = c(paste0("Sandy Soils (Lower ", qtxt, ")"),
								paste0("Water-Retentive Soils (Upper ", qtxt, ")")),
		rain_pc = c(paste0("Low Precipitation (Lower ", qtxt, ")"),
								paste0("High Precipitation (Upper ", qtxt, ")")),
		elevation = c(paste0("Low Elevation (Lower ", qtxt, ")"),
									paste0("High Elevation (Upper ", qtxt, ")")),
		soil_ph = c(paste0("Low Soil pH (Lower ", qtxt, ")"),
								paste0("High Soil pH (Upper ", qtxt, ")"))
	)
	
	temp <- fit_models_on_strata_biome(boot_B, "temp_pc", traits, covariates, hyper_grid, num_threads, qlabels$temp_pc, probs)
	soil <- fit_models_on_strata_biome(boot_B, "soil_pc", traits, covariates, hyper_grid, num_threads, qlabels$soil_pc, probs)
	rain <- fit_models_on_strata_biome(boot_B, "rain_pc", traits, covariates, hyper_grid, num_threads, qlabels$rain_pc, probs)
	elev <- fit_models_on_strata_biome(boot_B, "elevation", traits, covariates, hyper_grid, num_threads, qlabels$elevation, probs)
	ph   <- fit_models_on_strata_biome(boot_B, "soil_ph", traits, covariates, hyper_grid, num_threads, qlabels$soil_ph, probs)
	
	performance <- dplyr::bind_rows(temp$performance, soil$performance, rain$performance, elev$performance, ph$performance) %>%
		dplyr::select(trait, group, rsq, pred_error)
	
	pdp_temp <- calculate_pdp_for_scenario(temp, traits, "standage", qlabels$temp_pc) %>% dplyr::mutate(variable = "Temperature (PC)")
	pdp_soil <- calculate_pdp_for_scenario(soil, traits, "standage", qlabels$soil_pc) %>% dplyr::mutate(variable = "Soil - Water Retention (PC)")
	pdp_rain <- calculate_pdp_for_scenario(rain, traits, "standage", qlabels$rain_pc) %>% dplyr::mutate(variable = "Precipitation (PC)")
	pdp_elev <- calculate_pdp_for_scenario(elev, traits, "standage", qlabels$elevation) %>% dplyr::mutate(variable = "Elevation")
	pdp_ph   <- calculate_pdp_for_scenario(ph, traits, "standage", qlabels$soil_ph) %>% dplyr::mutate(variable = "Soil pH")
	
	pdp_all <- dplyr::bind_rows(pdp_temp, pdp_soil, pdp_rain, pdp_elev, pdp_ph) %>%
		dplyr::left_join(performance, by = c("trait","group"))
	
	pdp_all
}

# Main bootstrap over all biomes
run_pdp_bootstrap_by_biome <- function(data, traits, covariates, hyper_grid,
																			 num_bootstrap = 100, num_threads = 1,
																			 probs = c(0.25,0.75)) {
	biomes <- sort(unique(data$biome_group))
	
	foreach(iter = 1:num_bootstrap, .combine = dplyr::bind_rows,
					.packages = c("dplyr","purrr","tidyr","stringr","ranger","rsample","pdp"),
					.export   = c("pdp_biome_once", "fit_models_on_strata_biome",
												"calculate_pdp_for_scenario", "calculate_partial_dependence",
												"stratify_within", "recode_group", "fit_rf_model")) %dopar% {
													purrr::map_dfr(biomes, function(B) {
														dB <- dplyr::filter(data, biome_group == B)
														if (nrow(dB) < 100) return(NULL)
														boot_B <- dB %>% dplyr::sample_frac(stats::runif(1, 0.6, 0.8), replace = TRUE)
														pdp_all <- pdp_biome_once(boot_B, traits, covariates, hyper_grid, num_threads, probs)
														pdp_all %>% dplyr::mutate(iteration = iter, biome = B)
													})
												}
}

# ----------- Predictability -----------







# XXXXX Switch to VEcv XXXXXXXXXX





