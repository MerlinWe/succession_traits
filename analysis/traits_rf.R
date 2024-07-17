############################################################################################################################
########################################  MSc Diss. Forest Succession Data Analysis ########################################  
############################################################################################################################

rm(list = ls()) # make sure environment is clean 
set.seed(42)    # set seed for reproducibility

# Load necessary libraries
library(caret)
library(ranger)
library(pdp)
library(mgcv)
library(rsample)
library(doParallel)
library(ggbeeswarm)
library(patchwork)
library(dggridR)
library(grid)
library(glue)
library(cowplot)
library(tidyverse)

# Spatial downsampling? 
downsample <- TRUE 

# Check which device is running
node_name <- Sys.info()["nodename"]

# For parallel processing: Register cores - 32 if on threadripper; 8 if local
num_cores <-  ifelse(node_name == "threadeast", 32, 8)
cl <- makeCluster(num_cores)
registerDoParallel(cl, cores = num_cores)
getDoParWorkers()

# Set paths conditionally
path_in <- ifelse(node_name == "threadeast", 
									"/home/merlin/RDS_drive/merlin/data/fia_traits", 
									"/Volumes/ritd-ag-project-rd01pr-dmayn10/merlin/data/fia_traits")

path_out <- ifelse(node_name == "threadeast", 
									 "/home/merlin/traits_output", 
									 "/Users/serpent/Documents/MSc/Thesis/Code/analysis")

# ---------- Read and prepare input data ----------

# Read data and make some minor adjustments
data <- read_csv(paste0(path_in, "/plotlevel_data_2024-07-15.csv")) %>%
	
	# Keep only columns relevant for the traits analysis 
	select(starts_with("wmean_"), standage, biome, managed,
				 annual_mean_temperature, annual_precipitation, temperature_seasonality, mean_diurnal_range,
				 min_temperature_of_coldest_month, max_temperature_of_warmest_month,
				 elevation, pop_density, sand_content_015cm, soil_ph_015cm, water_capacity_015cm, LAT, LON) %>%
	filter(complete.cases(.)) %>%
	
	# Adjust some column names 
	rename_with(~ gsub("wmean_", "", .), starts_with("wmean_")) %>%
	rename_with(~ gsub("_015cm", "", .), ends_with("_015cm")) %>%
	rename(min_temperature = min_temperature_of_coldest_month,
				 max_temperature = max_temperature_of_warmest_month) %>%
	
	# Filter standage by the upper 10% quantiles
	filter(standage < quantile(standage, 0.9)) %>%
	# Disregard plots that have been managed 
	filter(managed == 0) %>%
	select(-managed)

# Create binary dummy variables based on biome levels with >100 observations; drop data of remaining biomes 

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

# Use spatial downsampling to speed up runtime and to deal with spatial autocorrelation? 
if (downsample) {
	grid_res <- 8
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
covariates <- c("standage", "annual_mean_temperature", "annual_precipitation", 
								"temperature_seasonality", "mean_diurnal_range", "min_temperature", 
								"max_temperature", "elevation", "pop_density", "sand_content", "soil_ph", 
								"water_capacity", "biome_boreal_forests_or_taiga", "biome_flooded_grasslands",        
								"biome_mediterranean_woodlands", "biome_temperate_broadleaf_forests", "biome_tundra",
								"biome_temperate_conifer_forests", "biome_temperate_grasslands", "biome_xeric_shrublands")


# Define a grid of hyperparameters including num.trees, mtry, and min.node.size
hyper_grid <- expand.grid(
	num.trees = c(500, 1000, 1500),
	mtry = 2:4,
	min.node.size = c(1, 10, 20))

# Define function to perform cross-validation and find the best hyperparameters
tune_rf_model <- function(trait, data, covariates, hyper_grid) {
	
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
			return(mod$prediction.error * 100)
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
		theme(text = element_text(family = "Arial", size = 6),
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

# Get tuning results 
tuning_result <- map(traits, ~ tune_rf_model(.x, train_data, covariates, hyper_grid))
names(tuning_result) <- traits # assign trait names 
hyper_grid <- map(tuning_result, "best_hyperparameters") %>% bind_rows()
tuning_error_plot <- plot_grid(plotlist = map(tuning_result, "error_plot"), ncol = 2, nrow = 4)
ggsave(paste0(path_out, "/plots/tuning_plot.png"),
			 plot = tuning_error_plot,
			 bg = "white",
			 width = 200,  
			 height = 130, 
			 units = "mm",
			 dpi = 600)

# Now fit trait models
trait_mods <- map(traits, ~ fit_rf_model(.x, train_data, covariates, hyper_grid))

# Examine performance of hyper parameters
performance_metrics <- map(trait_mods, "performance") %>% bind_rows()
best_models <- map(trait_mods, "trait_mod")

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
	"annual precipitation" = "AP",
	"annual mean temperature" = "AMT",
	"elevation" = "EL",
	"standage" = "SA",
	"temperature seasonality" = "TS",
	"min temperature" = "MinT",
	"biome temperate conifer forests" = "BTC",
	"max temperature" = "MaxT",
	"soil ph" = "SpH",
	"sand content" = "SC",
	"mean diurnal range" = "MDR",
	"water capacity" = "WC",
	"pop density" = "PD",
	"biome temperate broadleaf forests" = "BTB",
	"biome temperate grasslands" = "BTG",
	"biome xeric shrublands" = "BXS",
	"biome mediterranean woodlands" = "BMW",
	"biome boreal forests or taiga" = "BBF",
	"biome tundra" = "BTU",
	"biome flooded grasslands" = "BFG")

# Function to create shapley plots for a single trait
plot_shapley_for_trait <- function(trait_id) {
	
	feature_importance %>%
		filter(trait == trait_id) %>%
		{
			# Plot shapley bee swarms
			p1 <- ggplot(., aes(x = reorder(feature, importance), y = shap, color = shap)) +
				geom_quasirandom(alpha = 0.5) +
				scale_color_viridis_c(option = "viridis",
															name = "Feature\nValue",
															breaks = c(min(feature_importance$shap) + 0.4 * (max(feature_importance$shap) - min(feature_importance$shap)),
																				 max(feature_importance$shap) - 0.4 * (max(feature_importance$shap) - min(feature_importance$shap))),
															labels = c("low", "high"),
															guide = "none") +
				scale_x_discrete(labels = feature_labels) +
				scale_y_continuous(n.breaks = 5) +
				coord_flip() +
				labs(x = NULL, 
						 y = "Shapley Value", 
						 title = glue("{distinct(., trait_name)$trait_name}")) +
				theme_bw(base_line_size = .3,
								 base_rect_size = .5) +
				theme(text = element_text(family = "Arial", size = 6),
							plot.title = element_text(face = "bold", size = 6),
							strip.text = element_text(face = "bold"),
							plot.margin = margin(t=7.5, b=5, r=1, l=5))
			
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
				theme(text = element_text(family = "Arial", size = 6),
							plot.title = element_text(size = 6),
							axis.text.y = element_blank(),
							axis.ticks.y = element_blank(),
							plot.margin = margin(t = 15, b = 2.5, r = 5, l = 1))
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
														min(feature_importance$shap) + 0.3 * (max(feature_importance$shap) - min(feature_importance$shap)), 
														max(feature_importance$shap)- 0.4 * (max(feature_importance$shap) - min(feature_importance$shap))),
													labels = c("low", "high"),
													guide = guide_colorbar(
														label.position = "right",
														barwidth = unit(3, "mm"),
														barheight = unit(100, "mm"))) +
		theme(legend.position = "right", 
					text = element_text(family = "Arial", size = 6),
					legend.title = element_text(size = 6),
					legend.text = element_text(size = 6),
					legend.margin = margin(t=0, b=-5, r=0, l=0)))

shap_plot <- plot_grid(
	plot_grid(
		plotlist = map(traits, plot_shapley_for_trait), 
		ncol = 4, nrow = 2, align = "v"),
	legend, ncol = 2, rel_widths = c(1, .07))

ggsave(paste0(path_out, "/plots/shap_plot.png"),
			 plot = shap_plot,
			 bg = "white",
			 width = 200,  
			 height = 130, 
			 units = "mm",
			 dpi = 600)

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

##### Climate quantiles #####

# Clean labels for traits
trait_labels <- c(
	wood_density = "Wood Density",
	bark_thickness = "Bark Thickness",
	conduit_diam = "Conduit Diameter",
	leaf_n = "Leaf Nitrogen",
	specific_leaf_area = "Specific Leaf Area",
	seed_dry_mass = "Seed Dry Mass",
	shade_tolerance = "Shade Tolerance",
	height = "Tree Height"
)

# >>> 10% Quantiles for temperature and precipitation  <<<

# Set Quantile values
quantiles_10 <- c(0.1, 0.9)

# Annual mean temperature 
temp_names_10 <- c("Coldest 10%", "Hottest 10%")
pdp_temp_10 <- get_partial_dependence(
	best_models = best_models,
	train_data = train_data,
	traits = traits,
	pred.vars = c("standage", "annual_mean_temperature"),
	quantiles = quantiles_10,
	quantile_names = temp_names_10)

# Annual Precipitation 
prcp_names_10 <- c("Driest 10%", "Wettest 10%")
pdp_prcp_10 <- get_partial_dependence(
	best_models = best_models,
	train_data = train_data,
	traits = traits,
	pred.vars = c("standage", "annual_precipitation"),
	quantiles = quantiles_10,
	quantile_names = prcp_names_10)

# Plot curves 
climate_10 <- plot_grid(
	
	# Temperature
	pdp_temp_10 %>%
	ggplot(aes(x = standage, y = yhat, color = group, shape = group)) +
		geom_ribbon(aes(ymin = lower, ymax = upper, fill = group), alpha = 0.2, linewidth = .3) +
		geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, linewidth = .7) +
		scale_color_manual(values = setNames(c("darkslateblue", "darkred"), temp_names_10)) +
		scale_fill_manual(values = setNames(c("darkslateblue", "darkred"), temp_names_10)) +
		labs(x = "Standage", y = "Predicted Trait Value (log-scaled)", color = "Quantile", shape = "Quantile", fill = "Quantile") + 
		theme_bw(base_line_size = .3,
						 base_rect_size = .5) +
		theme(legend.position = "top",
					legend.justification = "center",
					text = element_text(family = "Arial", size = 8),
					legend.title = element_blank(),
					strip.background.y = element_blank(),
					strip.text.y = element_blank(),
					strip.background = element_rect(fill = "white", color = "black", linewidth = .75),
					plot.margin = margin(t=0, b=5, l=1, r=5),
					legend.margin = margin(t=5, b=-3, l=0, r=0),
					legend.key.size = unit(0.5, "cm"),
					legend.text = element_text(size = 7, margin = margin(l = -25, r = 10)),
					legend.spacing.x = unit(1, "cm"), 
					legend.spacing.y = unit(0.5, "cm")) + 
		facet_grid(trait ~ group, scales = "free", labeller = labeller(trait = trait_labels)),
	
	# Precipitation
	pdp_prcp_10 %>%
		ggplot(aes(x = standage, y = yhat, color = group, shape = group)) +
		geom_ribbon(aes(ymin = lower, ymax = upper, fill = group), alpha = 0.2, linewidth = .3) +
		geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, linewidth = .7) +
		scale_color_manual(values = setNames(c("tan4", "cyan4"), prcp_names_10)) +
		scale_fill_manual(values = setNames(c("tan4", "cyan4"), prcp_names_10)) +
		labs(x = "Standage", y = NULL, color = "Quantile", shape = "Quantile", fill = "Quantile") + 
		theme_bw(base_line_size = .3,
						 base_rect_size = .5) +
		theme(legend.position = "top",
					legend.justification = "center",
					text = element_text(family = "Arial", size = 8),
					legend.title = element_blank(),
					strip.text.y = element_text(size = 8),
					strip.background = element_rect(fill = "white", color = "black", linewidth = .75),
					plot.margin = margin(t=0, b=5, l=1, r=5),
					legend.margin = margin(t=5, b=-3, l=0, r=0),
					legend.key.size = unit(0.5, "cm"),
					legend.text = element_text(size = 7, margin = margin(l = -25, r = 10)),
					legend.spacing.x = unit(1, "cm"), 
					legend.spacing.y = unit(0.5, "cm")) + 
		facet_grid(trait ~ group, scales = "free", labeller = labeller(trait = trait_labels)),
	
	ncol = 2, 
	nrow = 1, 
	rel_widths = c(1, 1.05),
	align = "h", 
	axis = c("l"), 
	greedy = TRUE)

ggsave(filename = paste0(path_out, "/plots/prcp_10.png"),
			 plot = climate_10, 
			 bg = "white",
			 width = 200, 
			 height = 250, 
			 units = "mm", 
			 dpi = 600)

# >>> 25% Quantiles for Temperature <<<

# Set Quantile values
quantiles_25 <- c(0.25, 0.75)

# Annual mean temperature 
temp_names_25 <- c("Coldest 25%", "Hottest 25%")
pdp_temp_25 <- get_partial_dependence(
	best_models = best_models,
	train_data = train_data,
	traits = traits,
	pred.vars = c("standage", "annual_mean_temperature"),
	quantiles = quantiles_25,
	quantile_names = temp_names_25)

# Annual Precipitation 
prcp_names_25 <- c("Driest 25%", "Wettest 25%")
pdp_prcp_25 <- get_partial_dependence(
	best_models = best_models,
	train_data = train_data,
	traits = traits,
	pred.vars = c("standage", "annual_precipitation"),
	quantiles = quantiles_25,
	quantile_names = prcp_names_25)

# Plot curves 
climate_25 <- plot_grid(
	
	# Temperature
	pdp_temp_25 %>%
		ggplot(aes(x = standage, y = yhat, color = group, shape = group)) +
		geom_ribbon(aes(ymin = lower, ymax = upper, fill = group), alpha = 0.2, linewidth = .3) +
		geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, linewidth = .7) +
		scale_color_manual(values = setNames(c("darkslateblue", "darkred"), temp_names_25)) +
		scale_fill_manual(values = setNames(c("darkslateblue", "darkred"), temp_names_25)) +
		labs(x = "Standage", y = "Predicted Trait Value (log-scaled)", color = "Quantile", shape = "Quantile", fill = "Quantile") + 
		theme_bw(base_line_size = .3,
						 base_rect_size = .5) +
		theme(legend.position = "top",
					legend.justification = "center",
					text = element_text(family = "Arial", size = 8),
					legend.title = element_blank(),
					strip.background.y = element_blank(),
					strip.text.y = element_blank(),
					strip.background = element_rect(fill = "white", color = "black", linewidth = .75),
					plot.margin = margin(t=0, b=5, l=1, r=5),
					legend.margin = margin(t=5, b=-3, l=0, r=0),
					legend.key.size = unit(0.5, "cm"),
					legend.text = element_text(size = 7, margin = margin(l = -25, r = 10)),
					legend.spacing.x = unit(1, "cm"), 
					legend.spacing.y = unit(0.5, "cm")) + 
		facet_grid(trait ~ group, scales = "free", labeller = labeller(trait = trait_labels)),
	
	# Precipitation
	pdp_prcp_25 %>%
		ggplot(aes(x = standage, y = yhat, color = group, shape = group)) +
		geom_ribbon(aes(ymin = lower, ymax = upper, fill = group), alpha = 0.2, linewidth = .3) +
		geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, linewidth = .7) +
		scale_color_manual(values = setNames(c("tan4", "cyan4"), prcp_names_25)) +
		scale_fill_manual(values = setNames(c("tan4", "cyan4"), prcp_names_25)) +
		labs(x = "Standage", y = NULL, color = "Quantile", shape = "Quantile", fill = "Quantile") + 
		theme_bw(base_line_size = .3,
						 base_rect_size = .5) +
		theme(legend.position = "top",
					legend.justification = "center",
					text = element_text(family = "Arial", size = 8),
					legend.title = element_blank(),
					strip.text.y = element_text(size = 8),
					strip.background = element_rect(fill = "white", color = "black", linewidth = .75),
					plot.margin = margin(t=0, b=5, l=1, r=5),
					legend.margin = margin(t=5, b=-3, l=0, r=0),
					legend.key.size = unit(0.5, "cm"),
					legend.text = element_text(size = 7, margin = margin(l = -25, r = 10)),
					legend.spacing.x = unit(1, "cm"), 
					legend.spacing.y = unit(0.5, "cm")) + 
		facet_grid(trait ~ group, scales = "free", labeller = labeller(trait = trait_labels)),
	
	ncol = 2, 
	nrow = 1, 
	rel_widths = c(1, 1.05),
	align = "h", 
	axis = c("l"), 
	greedy = TRUE)

ggsave(filename = paste0(path_out, "/plots/prcp_25.png"),
			 plot = climate_25, 
			 bg = "white",
			 width = 200, 
			 height = 250, 
			 units = "mm", 
			 dpi = 600)

## >>>>> Now look at different biomes <<<<<

# General function to create partial dependence data for each trait based on biomes
get_partial_dependence_biomes <- function(best_models, train_data, traits, pred.var, biome_vars, biome_names = NULL) {
	
	if (is.null(biome_names)) {
		biome_names <- biome_vars
	}
	
	# Print the biome names
	for (biome in biome_vars) {
		cat("Processing:", biome, "\n")
	}
	
	# Separate data for each biome
	data_list <- lapply(biome_vars, function(biome) {
		train_data %>% filter(!!sym(biome) == 1)
	})
	
	# Perform partial dependence calculation in parallel
	pdp_data_combined <- foreach(trait = traits, .combine = 'rbind', .packages = c('pdp', 'dplyr')) %dopar% {
		model <- best_models[[which(traits == trait)]]
		
		# Combine partial dependence data for each biome
		pdp_data_list <- lapply(seq_along(biome_vars), function(i) {
			data_subset <- data_list[[i]]
			
			# Build grid for partial dependence
			pred_grid <- expand.grid(
				standage = unique(data_subset$standage))
			
			# Calculate partial dependence
			pdp_data <- pdp::partial(
				object = model,
				pred.var = pred.var,
				pred.grid = pred_grid,
				train = data_subset,
				plot = FALSE)
			
			# Label the data with the biome name
			pdp_data %>% mutate(group = biome_names[i])
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

# Set biome names and variable names; we also consider if a system is dominated by gymnosperms or angiosperms
biomes <- tibble(
	biome_vars = c("biome_boreal_forests_or_taiga", "biome_flooded_grasslands", "biome_mediterranean_woodlands", 
								 "biome_temperate_broadleaf_forests", "biome_temperate_conifer_forests", "biome_temperate_grasslands", 
								 "biome_tundra", "biome_xeric_shrublands"),
	biome_names = c("Boreal Forests/Taiga", "Flooded Grasslands", "Mediterranean Woodlands", 
									"Temperate Broadleaf", "Temperate Conifer", "Temperate Grasslands", 
									"Tundra", "Xeric Shrublands")) %>% 
	mutate(
	dominance = c("Gymnosperm", "Angiosperm", "Angiosperm", 
								"Angiosperm", "Gymnosperm", "Angiosperm", 
								"Gymnosperm", "Gymnosperm"))

# Calculate partial dependence towards standage for all traits for the biomes
pdp_biomes <- get_partial_dependence_biomes(
	best_models = best_models,
	train_data = train_data,
	traits = traits,
	pred.var = "standage",
	biome_vars = biomes$biome_vars,
	biome_names = biomes$biome_names) %>%
	# Add biome dominance
	left_join(
		biomes %>% 
			rename(group = biome_names) %>% 
			dplyr::select(-biome_vars),
		by = "group")

# Plot partial dependence across biomes 
dominance_colors <- c("Angiosperm" = "forestgreen", "Gymnosperm" = "saddlebrown")
biome_plot <- ggplot(pdp_biomes, aes(x = standage, y = yhat, color = dominance, fill = dominance, group = group)) +
	geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .1, linewidth = .1) +
	geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, linewidth = .7) +
	scale_color_manual(values = dominance_colors) +
	scale_fill_manual(values = dominance_colors) +
	labs(x = "Standage", y = "Predicted Trait Value (log-scaled)", color = "Dominance", fill = "Dominance", title = NULL) + 
	theme_bw() +
	theme(
		legend.position = "top",
		legend.justification = "left",
		text = element_text(family = "Arial", size = 8),
		strip.text.y = element_text(size = 6.5),
		strip.text.x = element_text(size = 7),
		strip.background = element_rect(fill = "white", color = "black", linewidth = .75),
		plot.margin = margin(t=5, b=5, l=5, r=5),
		legend.margin = margin(t=5, b=-3, l=0, r=0)) +
	facet_grid(trait ~ group, scales = "free_y", labeller = labeller(trait = trait_labels))

ggsave(filename = paste0(path_out, "/plots/biomes_plot.png"),
			 plot = biome_plot, 
			 bg = "white",
			 width = 200, 
			 height = 200, 
			 units = "mm", 
			 dpi = 300)

## Done - stop the cluster
stopCluster(cl)
