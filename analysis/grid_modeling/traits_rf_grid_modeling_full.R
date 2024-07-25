############################################################################################################################
########################################  MSc Diss. Forest Succession Data Analysis ########################################  
############################################################################################################################

rm(list = ls()); gc() # make sure environment is clean 
set.seed(42)    # set seed for reproducibility

# ----- Session setup -----

# Load necessary libraries
library(caret)
library(ranger)
library(pdp)
library(mgcv)
library(rsample)
library(doParallel)
library(ggbeeswarm)
library(forcats)
library(dggridR)
library(grid)
library(glue)
library(cowplot)
library(psych)
library(factoextra)
library(tidyverse)

# Export plots?
export <- FALSE

# Parallize?
parallel <- TRUE 

# Spatial downsampling? If TRUE set resolution!
downsample <- TRUE 
down_res <- 7

# Check which device is running
node_name <- Sys.info()["nodename"]

# Set file paths conditionally
path_in <- ifelse(node_name == "threadeast", 
									"/home/merlin/RDS_drive/merlin/data/fia_traits", 
									"/Volumes/ritd-ag-project-rd01pr-dmayn10/merlin/data/fia_traits")

path_out <- ifelse(node_name == "threadeast", 
									 "/home/merlin/traits_output", 
									 "/Users/serpent/Documents/MSc/Thesis/Code/analysis")

if (parallel) {
	# For parallel processing: Register cores - 32 if on threadripper; 8 if local
	num_cores <-  ifelse(node_name == "threadeast", 32, 8)
	cl <- makeCluster(num_cores)
	registerDoParallel(cl, cores = num_cores)
	getDoParWorkers()
}

# ---------- Read and prepare input data ----------

data <- read_csv(paste0(path_in, "/plotlevel_data_2024-07-15.csv")) %>%
	# Filter out standage above upper 10% quantiles and managed plots
	filter(standage < quantile(standage, 0.9)) %>%
	filter(managed == 0) 

# Create binary dummy variables based on biome levels with >100 observations; 
# Drop data of remaining biomes.

data <- data %>%
	filter(biome %in% (
		data %>%
			count(biome) %>%
			filter(n > 100) %>%
			pull(biome)))
data <- model.matrix(~ biome - 1, data = data) %>%
	as_tibble() %>%
	rename_with(~ .x %>%
								str_replace_all("biome(?!_)", "biome_") %>%
								str_replace_all(" ", "_") %>%
								str_to_lower()) %>%
	bind_cols(data, .) %>%
	select(-biome)

## ---------- PCA to find environmental axes ----------

## Fit a PCA model to find axes representing underlying climatic and environmental factors. 
## Disregard elevation and soil ph to include in the model and control for their effect
## directly, and since there is no reason to assume underlying factors. 

climate_vars <- data %>% 
	dplyr::select(
		
		annual_mean_temperature,
		max_temperature_of_warmest_month,
		min_temperature_of_coldest_month,
		mean_temperature_of_coldest_quarter,
		mean_temperature_of_warmest_quarter,

		annual_precipitation,
		precipitation_of_driest_quarter,
		precipitation_of_wettest_quarter,
		precipitation_of_coldest_quarter,
		precipitation_of_warmest_quarter,

		sand_content_015cm,
		sand_content_060cm,
		water_capacity_015cm,
		water_capacity_060cm)

KMO(climate_vars) # check KMO
cortest.bartlett(cor(climate_vars), n = ncol(climate_vars)) # check sphericity

# Run PCA
climate_pca <- prcomp(climate_vars, scale. = TRUE, center = TRUE)
summary(climate_pca) 

# Varimax rotation to link variables with components
varimax_rotation <- varimax(climate_pca$rotation[, 1:3])
varimax_rotation <- varimax_rotation$loadings
print(varimax_rotation)

# Filter scores based on varimax loadings
scores <- round(get_pca_var(climate_pca)$coord[,1:3], digits = 3)
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

# Assign scores based on varimax rotated loadings
scores <- filter_scores(scores, varimax_rotation, threshold = 0.2) %>%
	as_tibble(rownames = "variable")

## Plot PCA for appendix 
abbreviations <- c(
	"annual_mean_temperature" = "AMT",
	"max_temperature_of_warmest_month" = "MaxTempWM",
	"min_temperature_of_coldest_month" = "MinTempCM",
	"mean_temperature_of_coldest_quarter" = "MeanTempCQ",
	"mean_temperature_of_warmest_quarter" = "MeanTempWQ",
	"annual_precipitation" = "AP",
	"precipitation_of_driest_quarter" = "PrecDQ",
	"precipitation_of_wettest_quarter" = "PrecWeQ",
	"precipitation_of_coldest_quarter" = "PrecCQ",
	"precipitation_of_warmest_quarter" = "PrecWaQ",
	"sand_content_015cm" = "Sand015",
	"sand_content_060cm" = "Sand060",
	"water_capacity_015cm" = "WaterCap015",
	"water_capacity_060cm" = "WaterCap060")

# Scree plot
scree_plot <- fviz_eig(climate_pca)
scree_plot <- scree_plot$data %>%
	as_tibble() %>%
	ggplot(aes(x = reorder(dim, -eig), y = eig, group = 1)) +  
	geom_bar(stat = "identity", color = "black", fill = "#1B7E74", alpha = 0.8) +
	geom_point() +
	geom_line(stat = "identity") +
	geom_vline(xintercept = 3.5, linetype = "dashed", colour = "red") +
	labs(y = "Percentage of Explained Variance", x = "Principal Components", title = "PC Eigenvalues") +
	theme_bw() +
	theme(text = element_text(family = "sans", size = 8))

# Function to create contribution plots with abbreviated x-axis labels
create_contrib_plot <- function(var_contrib_data, title, fill_color) {
	var_contrib_data %>%
		as_tibble() %>%
		mutate(name = abbreviations[name]) %>%
		ggplot(aes(x = reorder(name, -contrib), y = contrib)) +
		geom_bar(stat = "identity", color = "black", fill = fill_color, alpha = 0.8) +
		geom_hline(yintercept = 7.5, colour = "black", linetype = "dashed") +
		labs(y = "Contributions (%)", x = NULL, title = title) +
		theme_bw() +
		theme(text = element_text(family = "sans", size = 8),
					axis.text.x = element_text(angle = 45, hjust = 1))
}

# Contribution plots for PC1, PC2, and PC3 with abbreviations
var_contrib_pc1 <- fviz_contrib(climate_pca, choice = "var", axes = 1, top = 10)
var_contrib_pc1_plot <- create_contrib_plot(var_contrib_pc1$data, "PC1 - Temperature", "firebrick4")

var_contrib_pc2 <- fviz_contrib(climate_pca, choice = "var", axes = 2, top = 10)
var_contrib_pc2_plot <- create_contrib_plot(var_contrib_pc2$data, "PC2 - Soil / Water Retention", "tan4")

var_contrib_pc3 <- fviz_contrib(climate_pca, choice = "var", axes = 3, top = 10)
var_contrib_pc3_plot <- create_contrib_plot(var_contrib_pc3$data, "PC3 - Precipitation", "dodgerblue3")

pca_plot <- plot_grid(scree_plot, var_contrib_pc1_plot, var_contrib_pc2_plot, var_contrib_pc3_plot,
											ncol = 2, nrow = 2, labels = "auto", label_fontfamily = "sans", label_size = 8)

if (export) {
	ggsave(paste0(path_out, "/plots/appendix/pca_plot.png"),
				 plot = pca_plot,
				 bg = "white",
				 width = 200,  
				 height = 130, 
				 units = "mm",
				 dpi = 600)
	
	write_csv(scores, file = paste0(path_out, "/tables/pca_scores.csv"))
}

# Write components
climate_pca <- climate_pca$x %>% 
	as_tibble() %>% 
	select(PC1, PC2, PC3) %>%
	rename(temp_pc = PC1,
				 soil_pc = PC2,
				 rain_pc = PC3) %>%
	mutate(PID_rep = data$PID_rep)

# Build database for further analysis 
data <- data %>% 
	left_join(climate_pca, by = "PID_rep") %>%
	select(starts_with("wmean_"), standage, temp_pc, soil_pc, rain_pc,
				 elevation, soil_ph_015cm, 
				 biome_boreal_forests_or_taiga, biome_flooded_grasslands, biome_mediterranean_woodlands, biome_temperate_broadleaf_forests,
				 biome_temperate_conifer_forests, biome_temperate_grasslands,biome_tundra, biome_xeric_shrublands, LAT, LON) %>%
	rename_with(~ gsub("_015cm", "", .), ends_with("_015cm")) %>%
	rename_with(~ gsub("wmean_", "", .), starts_with("wmean_")) %>%
	filter(complete.cases(.))

rm(climate_pca, climate_vars, pca_plot, scores, scree_plot, var_contrib_pc1, 
	 var_contrib_pc1_plot, var_contrib_pc2, var_contrib_pc2_plot, var_contrib_pc3, 
	 var_contrib_pc3_plot, abbreviations, varimax_rotation); gc()

# ----- Spatial downsampling (if TRUE) -----

# Use spatial downsampling to speed up runtime and to deal with spatial autocorrelation? 
if (downsample) {
	grid_res <- down_res
	dggs <- dgconstruct(res = grid_res, metric = TRUE, resround = 'nearest')
	
	data <- data %>%
		mutate(cell = dgGEO_to_SEQNUM(dggs, in_lon_deg = LON, in_lat_deg = LAT)$seqnum) %>%
		dplyr::select(-LAT, -LON) %>%
		group_by(cell) %>%
		dplyr::summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = 'drop')
}
	
## ---------- Model training and validation ----------

# Split into training and test sets
split <- initial_split(data, prop = 0.8)

# Extract training and testing datasets
train_data <- training(split)
test_data  <- testing(split)

# Define traits
traits <- c("wood_density", "bark_thickness", "conduit_diam", "leaf_n", 
						"specific_leaf_area", "seed_dry_mass", "shade_tolerance", "height")

# Define covariates 
covariates <- c("standage", "temp_pc", "soil_pc", "rain_pc", "elevation", "soil_ph", "biome_boreal_forests_or_taiga", 
								"biome_flooded_grasslands", "biome_mediterranean_woodlands", "biome_temperate_broadleaf_forests",
								"biome_temperate_conifer_forests", "biome_temperate_grasslands","biome_tundra", "biome_xeric_shrublands")

# Define a grid of hyperparameters including num.trees, mtry, and min.node.size
hyper_grid <- expand.grid(
	num.trees = c(500, 1000, 1500),
	mtry = 2:4,
	min.node.size = c(1, 10, 20))

# Define function to perform cross-validation and find the best hyperparameters
tune_rf_model <- function(trait, data, covariates, hyper_grid, labels) {
	
	# Define model formula
	formula <- as.formula(paste(trait, "~", paste(covariates, collapse = " + ")))
	
	# First: Find best set of hyperparameters in parallel processing
	prediction_error <- foreach(
		
		num.trees = hyper_grid$num.trees,
		mtry = hyper_grid$mtry,
		min.node.size = hyper_grid$min.node.size,
		.combine = 'c', 
		.packages = "ranger") %dopar% {
			
			# Fit model
			mod <- ranger::ranger(
				formula = formula,
				data = data,
				dependent.variable.name = trait,
				num.trees = num.trees,
				mtry = mtry,
				min.node.size = min.node.size,
				num.threads = 1)
			
			# Returning prediction
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
		dplyr::arrange(prediction_error) %>% 
		dplyr::slice(1)
	
	return(list(error_plot = error_plot, best_hyperparameters = best_hyperparameters))
}

# Define function to fit a random forest with the best hyperparameters per trait
fit_rf_model <- function(traits, data, covariates, hyper_parameters) {
	
	# Set threads to 4 if on threadeast (so 32 cores are used with 8 traits) or 1 if locally 
	num_threads <- ifelse(node_name == "threadeast", 4, 1)
	
	# Fit trait models 
	trait_mod <- foreach(
		
		trait = traits,
		.combine = 'c', 
		.packages = c("ranger", "dplyr")
		
	) %dopar% {
		
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
		
		# Returning model and performance metrics 
		return(list(trait_mod = trait_mod, performance = performance))
	}
	
	return(trait_mod)
}

# Clean memory before parallel processing
rm(dggs, split); gc()

# Get tuning results 
tuning_result <- map(traits, ~ tune_rf_model(.x, train_data, covariates, hyper_grid))
names(tuning_result) <- traits # assign trait names 
hyper_grid <- map(tuning_result, "best_hyperparameters") %>% bind_rows()
tuning_error_plot <- plot_grid(plotlist = map(tuning_result, "error_plot"), ncol = 2, nrow = 4)

if (export) {
ggsave(paste0(path_out, "/plots/appendix/tuning_plot.png"),
			 plot = tuning_error_plot,
			 bg = "white",
			 width = 200,  
			 height = 130, 
			 units = "mm",
			 dpi = 600)
}

# Clean memory before parallel processing
rm(tuning_error_plot, tuning_result); gc()

# Now fit trait models
best_models <- map(traits, ~ fit_rf_model(.x, train_data, covariates, hyper_grid))

# Examine performance of hyper parameters
performance_metrics <- map(best_models, "performance") %>% bind_rows()
best_models <- map(best_models, "trait_mod")

# Clean memory before parallel processing
rm(hyper_grid); gc()

## ---------- Calculate Shapley values for each trait model ----------

# Prediction function for fastshap
predict_fn <- function(object, newdata) {
	predict(object, data = newdata)$predictions
}

# Loop through final trait models to calculate shapley values 
shap_values_list <- list()

# Loop over each model and trait with inner parallelization
for (i in seq_along(best_models)) {
	model <- best_models[[i]]
	trait <- traits[i]
	
	shap_values <- fastshap::explain(model, 
																	 X = test_data %>% select(all_of(covariates)) %>% as.matrix(),
																	 pred_wrapper = predict_fn, 
																	 nsim = 10, 
																	 parallel = TRUE)  
	
	shap_values <- as.data.frame(shap_values)
	shap_values$trait <- trait
	
	# Append to the list
	shap_values_list[[i]] <- shap_values
}

# Combine all results into a single data frame
shap_values <- do.call(rbind, shap_values_list) %>% as_tibble()

## Visualize shapley values; first get a tibble for plotting
shap_long <- shap_values %>%
	pivot_longer(cols = -c(trait), names_to = "feature", values_to = "shap") %>%
	as_tibble() %>%
	mutate(feature = gsub("_", " ", feature),
				 trait_name = case_when(
				 	trait == "wood_density" ~ "Wood Density",
				 	trait == "bark_thickness" ~ "Bark Thickness",
				 	trait == "conduit_diam" ~ "Conduit Diameter",
				 	trait == "leaf_n" ~ "Leaf Nitrogen",
				 	trait == "specific_leaf_area" ~ "Specific Leaf Area",
				 	trait == "seed_dry_mass" ~ "Seed Dry Mass",
				 	trait == "shade_tolerance" ~ "Shade Tolernce",
				 	trait == "height" ~ "Tree Height",
				 	TRUE ~ NA_character_)) 

# Sum absolute Shapley values to determine overall importance
feature_importance <- shap_long %>%
	group_by(trait, feature) %>%
	summarize(importance = sum(abs(shap)), .groups = "drop") %>%
	arrange(trait, importance) %>%
	right_join(shap_long, by = c("trait", "feature")) %>%
	left_join(
		performance_metrics %>%
			mutate(trait_name = case_when(
				trait == "wood_density" ~ "Wood Density",
				trait == "bark_thickness" ~ "Bark Thickness",
				trait == "conduit_diam" ~ "Conduit Diameter",
				trait == "leaf_n" ~ "Leaf Nitrogen",
				trait == "specific_leaf_area" ~ "Specific Leaf Area",
				trait == "seed_dry_mass" ~ "Seed Dry Mass",
				trait == "shade_tolerance" ~ "Shade Tolernce",
				trait == "height" ~ "Tree Height",
				TRUE ~ NA_character_)) %>%
			dplyr::select(trait, rsq) %>%
			mutate(rsq = round(rsq, digits = 3)) %>%
			distinct(),
		by = "trait") 

# Shorten feature labels
feature_labels <- c(
	"standage" = "STAGE",
	"temp pc" = "TEMP (PC)",
	"soil pc" = "SOIL (PC)",
	"rain pc"= "PRCP (PC)",
	"elevation" = "ELEV",
	"soil ph" = "pH",
	"biome temperate conifer forests" = "B-TC",
	"biome temperate broadleaf forests" = "B-TB",
	"biome temperate grasslands" = "B-TG",
	"biome xeric shrublands" = "B-XS",
	"biome mediterranean woodlands" = "B-MW",
	"biome boreal forests or taiga" = "B-BF",
	"biome tundra" = "B-TU",
	"biome flooded grasslands" = "B-FG")

# Function to create shapley plots for a single trait
plot_shapley_for_trait <- function(trait_id, hide_x_axis_labels = FALSE) {
	
	feature_importance %>%
		filter(trait == trait_id) %>%
		{
			shap_min <- min(.$shap)
			shap_max <- max(.$shap)
			shap_range <- shap_max - shap_min
			y_lim_min <- shap_min - 0.1 * shap_range
			y_lim_max <- shap_max + 0.1 * shap_range
			
			# Plot shapley bee swarms
			p1 <- ggplot(., aes(x = reorder(feature, importance), y = shap, color = shap)) +
				geom_quasirandom(alpha = 0.5) +
				scale_color_viridis_c(option = "viridis",
															name = "Feature\nValue",
															breaks = c(shap_min + 0.4 * shap_range,
																				 shap_max - 0.4 * shap_range),
															labels = c("low", "high"),
															guide = "none") +
				scale_x_discrete(labels = feature_labels) +
				scale_y_continuous(n.breaks = 5, limits = c(y_lim_min, y_lim_max)) +
				coord_flip() +
				labs(x = NULL, 
						 y = "Shapley Value", 
						 title = glue("{distinct(., trait_name)$trait_name}")) +
				theme_bw(base_line_size = .3,
								 base_rect_size = .5) +
				theme(text = element_text(family = "sans", size = 6),
							plot.title = element_text(face = "bold", size = 6),
							strip.text = element_text(face = "bold"),
							plot.margin = margin(t=7.5, b=5, r=1, l=5)) +
				theme(axis.title.x = if(hide_x_axis_labels) element_blank() else element_text())
							
			# Plot absolute shapley importance 
			p2 <- ggplot(distinct(., feature, importance, rsq), aes(x = reorder(feature, importance), y = importance)) +
				geom_bar(stat = "identity", fill = "#1B7E74", color = "black", alpha = 0.7, linewidth = 0.2) +
				coord_flip() +
				scale_y_continuous(n.breaks = 5) +
				labs(x = NULL,
						 y = "Feature Importance",
						 title = glue("R² = {distinct(., rsq)$rsq}")) +
				theme_bw(base_line_size = .3,
								 base_rect_size = .5) +
				theme(text = element_text(family = "sans", size = 6),
							plot.title = element_text(size = 6),
							axis.text.y = element_blank(),
							axis.ticks.y = element_blank(),
							plot.margin = margin(t = 15, b = 2.5, r = 5, l = 1)) +
				theme(axis.title.x = if(hide_x_axis_labels) element_blank() else element_text())
			
			# Set plot structure 
			plot_grid(p1, p2, 
								ncol = 2, 
								nrow = 1, 
								rel_widths = c(.6, .4),
								align = "h", 
								axis = c("l"), 
								greedy = TRUE)
		}
}

# Add legend and build compound figure 
legend <- get_legend(
	feature_importance %>%
		filter(trait_name == "Bark Thickness") %>%
		ggplot(aes(x = reorder(feature, importance), y = shap, color = shap)) +
		geom_quasirandom(alpha = 0.5) +
		scale_color_viridis_c(option = "viridis",
													name = "Feature\nValue",
													breaks = c(
														min(feature_importance$shap) + 0.05 * (max(feature_importance$shap) - min(feature_importance$shap)), 
														max(feature_importance$shap)- 0.05 * (max(feature_importance$shap) - min(feature_importance$shap))),
													labels = c("low", "high"),
													guide = guide_colorbar(
														label.position = "right",
														barwidth = unit(3, "mm"),
														barheight = unit(100, "mm"))) +
		theme(legend.position = "right", 
					text = element_text(family = "sans", size = 6),
					legend.title = element_text(size = 6),
					legend.text = element_text(size = 6),
					legend.margin = margin(t=0, b=-5, r=0, l=0)))

# Create a list of plots with a flag to hide x-axis labels for the upper row
plots <- map2(traits, rep(c(TRUE, FALSE), each = length(traits) / 2), plot_shapley_for_trait)

shap_plot <- plot_grid(
	plot_grid(
		plotlist = plots, 
		ncol = 4, nrow = 2, align = "v"),
	legend, ncol = 2, rel_widths = c(1, .07))

if (export) {
	ggsave(paste0(path_out, "/plots/shap_plot.png"),
				 plot = shap_plot,
				 bg = "white",
				 width = 200,  
				 height = 130, 
				 units = "mm",
				 dpi = 600)
	
	write.csv(shap_values, file = paste0(path_out, "/tables/shap_values.csv"))
}

# Clean memory before parallel processing
rm(feature_importance, legend, model, shap_long, plots,
	 shap_plot, shap_values, shap_values_list, feature_labels, i, trait); gc()

## -------- Calculate partial dependence curves ----------

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
		quantile_names = names
	)
}

## >>> Look at environmental quantiles <<<

# Set Quantile values
quantiles_25 <- c(0.25, 0.75)

# Clean labels for traits
trait_labels <- c(
	wood_density = "Wood\nDens.",
	bark_thickness = "Bark\nThick.",
	conduit_diam = "Conduit\nDiam.",
	leaf_n = "Leaf\nNitrogen",
	specific_leaf_area = "SLA",
	seed_dry_mass = "Seed\nDry Mass",
	shade_tolerance = "Shade\nTol.",
	height = "Tree\nHeight")

# Define covariates and quantile names
variables <- list(
	temp_pc = c("Cold Temperatures (lower 25%)", "Warm Temperatures (upper 25%)"),
	soil_pc = c("Sandy Soils (lower 25%)", "Water-Retentive Soils (upper 25%)"),
	rain_pc = c("Low Precipitation (lower 25%)", "High Precipitation (upper 25%)"),
	elevation = c("Low Elevation (lower 25%)", "High Elevation (upper 25%)"),
	soil_ph = c("Low Soil pH (lower 25%)", "High soil pH (upper 25%)"))

# Get partial dependence for each covariate of
pdp_results <- lapply(names(variables), function(var) {
	get_pdp(var, variables[[var]])
})

# Assign results to individual variables
pdp_temp_25 <- pdp_results[[1]] 
pdp_soil_25 <- pdp_results[[2]]
pdp_prcp_25 <- pdp_results[[3]]
pdp_elev_25 <- pdp_results[[4]]
pdp_ph_25 <- pdp_results[[5]]
rm(pdp_results); gc() # clean-up 

## >>> Plot pdp curves for appendix <<<

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

# Define plot parameters and generate plots
plot_params <- list(
	list(data = pdp_temp_25, levels = variables$temp_pc, colors = c("darkslateblue", "darkred"), labels = variables$temp_pc, x_lab = NULL, y_lab = "Predicted Trait Value (log-scaled)"),
	list(data = pdp_prcp_25, levels = variables$rain_pc, colors = c("tan4", "cyan4"), labels = variables$rain_pc, x_lab = NULL, y_lab = NULL),
	list(data = pdp_soil_25, levels = variables$soil_pc, colors = c("khaki4", "springgreen3"), labels = variables$soil_pc, x_lab = "Standage", y_lab = NULL),
	list(data = pdp_elev_25, levels = variables$elevation, colors = c("bisque4", "royalblue4"), labels = variables$elevation, x_lab = "Standage", y_lab = "Predicted Trait Value (log-scaled)"),
	list(data = pdp_ph_25, levels = variables$soil_ph, colors = c("orangered1", "royalblue4"), labels = variables$soil_ph, x_lab = "Standage", y_lab = NULL))

pdp_plots <- lapply(plot_params, function(params) {
	create_pdp_plot(params$data, params$levels, params$colors, params$labels, params$x_lab, params$y_lab)
})

# Build compound figure 
pdp_25 <- plot_grid(
	plotlist = pdp_plots,
	ncol = 3,
	nrow = 2,
	rel_widths = c(1, 1, 1.05),
	align = "h",
	axis = c("l"),
	greedy = TRUE,
	labels = "auto",
	label_fontfamily = "sans",
	label_size = 8)

if (export) {
	ggsave(filename = paste0(path_out, "/plots/appendix/pdp_quantiles.png"),
				 plot = pdp_25, 
				 bg = "white",
				 width = 265, 
				 height = 205, 
				 units = "mm", 
				 dpi = 600)
	
	write.csv(pdp_temp_25, file = paste0(path_out, "/tables/pdp_temp_25.csv"))
	write.csv(pdp_soil_25, file = paste0(path_out, "/tables/pdp_soil_25.csv"))
	write.csv(pdp_prcp_25, file = paste0(path_out, "/tables/pdp_prcp_25.csv"))
	write.csv(pdp_elev_25, file = paste0(path_out, "/tables/pdp_elev_25.csv"))
	write.csv(pdp_ph_25,   file = paste0(path_out, "/tables/pdp_ph_25.csv"))
}

## Done with parallel processing - stop the cluster
stopCluster(cl)
rm(plot_params, pdp_plots, variables); gc()

## ---------- Conclusive analysis ----------

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

# Calculate slopes and deltas
slopes <- bind_rows(
	process_pdp_data(pdp_temp_25) %>% mutate(variable = "Temperature"),
	process_pdp_data(pdp_soil_25) %>% mutate(variable = "Soil"),
	process_pdp_data(pdp_prcp_25) %>% mutate(variable = "Precipitation"),
	process_pdp_data(pdp_elev_25) %>% mutate(variable = "Elevation"),
	process_pdp_data(pdp_ph_25) %>% mutate(variable = "pH"))


# Summarise for traits 
final <- slopes %>%
	group_by(trait) %>%
	summarize(
		mean_slope = mean(slope),
		sd_slope = sd(slope),
		mean_delta = mean(delta),
		sd_delta = sd(delta)) %>%
	ungroup() %>%
	left_join(performance_metrics %>% select(trait, rsq), by = "trait") %>%
	mutate(trait = case_when(
				 	trait == "wood_density" ~ "Wood Density",
				 	trait == "bark_thickness" ~ "Bark Thickness",
				 	trait == "conduit_diam" ~ "Conduit Diameter",
				 	trait == "leaf_n" ~ "Leaf Nitrogen",
				 	trait == "specific_leaf_area" ~ "Specific Leaf Area",
				 	trait == "seed_dry_mass" ~ "Seed Dry Mass",
				 	trait == "shade_tolerance" ~ "Shade Tolernce",
				 	trait == "height" ~ "Tree Height",
				 	TRUE ~ NA_character_)) 

# ------ Final Figure -----
shapes <- c(15, 16, 17, 18, 19, 15, 16, 17) 
final <- final %>%
	ggplot(aes(x = mean_delta, y = rsq, color = trait, shape = trait, label = trait)) +
	geom_point(size = 3, fill = "white") +
	geom_errorbar(aes(xmin = mean_delta - sd_delta, xmax = mean_delta + sd_delta), width = 0.005) +
	scale_colour_viridis_d() +
	scale_shape_manual(values = shapes) +
	labs(x = "Mean Difference Across Succession (log-scaled)",
			 y = "R²",
			 color = "Trait",
			 shape = "Trait") +
	geom_vline(xintercept = median(final$mean_delta), linetype = "dashed", color = "black") +
	geom_hline(yintercept = median(final$rsq), linetype = "dashed", color = "black") +
	theme_bw() +
	theme(legend.position = "right",
				text = element_text(family = "sans", size = 8))

# ------ For appendix -----
rest <- slopes %>%
	left_join(performance_metrics %>% select(trait, rsq), by = "trait") %>%
	mutate(group = factor(group, levels = c("Cold Temperatures (lower 25%)", "Warm Temperatures (upper 25%)",
																					"Low Precipitation (lower 25%)", "High Precipitation (upper 25%)",
																					"Sandy Soils (lower 25%)", "Water-Retentive Soils (upper 25%)",
																					"Low Elevation (lower 25%)", "High Elevation (upper 25%)",
																					"Low Soil pH (lower 25%)", "High soil pH (upper 25%)"))) %>%
	ggplot(aes(x = delta, y = rsq, color = trait, shape = trait, label = trait)) +
	geom_point(size = 3, fill = "white") +
	scale_colour_viridis_d() +
	scale_shape_manual(values = shapes) +
	labs(x = "Mean Difference Across Succession (log-scaled)",
			 y = "R²",
			 color = "Trait",
			 shape = "Trait") +
	geom_vline(xintercept = median(final$mean_delta), linetype = "dashed", color = "black") +
	geom_hline(yintercept = median(final$rsq), linetype = "dashed", color = "black") +
	theme_bw() +
	theme(legend.position = c(.75,.15),
				legend.direction = "horizontal",
				legend.title = element_blank(),
				legend.background = element_rect(fill = "transparent", colour = "transparent"),
				text = element_text(family = "sans", size = 8),
				strip.background = element_rect(fill = "white", color = "black", linewidth = .75)) +
	facet_wrap(~group, ncol = 4, nrow = 3)

if (export) {
	
	ggsave(filename = paste0(path_out, "/plots/succession_plot.png"),
				 plot = final, 
				 bg = "white",
				 width = 200, 
				 height = 170, 
				 units = "mm", 
				 dpi = 600)
	
	ggsave(filename = paste0(path_out, "/plots/appendix/quantiles_succession.png"),
				 plot = rest, 
				 bg = "white",
				 width = 200, 
				 height = 170, 
				 units = "mm", 
				 dpi = 600)
	
	write.csv(slopes, file = paste0(path_out, "/tables/succession.csv") )
}

rm(final, pdp_25, rest, slopes, shapes, trait_labels, pdp_ph_25, 
	 pdp_elev_25, pdp_prcp_25, pdp_soil_25, pdp_temp_25); gc()
