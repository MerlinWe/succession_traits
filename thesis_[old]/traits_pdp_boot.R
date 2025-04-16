############################################################################################################################
########################################  MSc Diss. Forest Succession Data Analysis ########################################  
############################################################################################################################

rm(list = ls())   # make sure environment is clean 
set.seed(42)      # set seed for reproducibility

# ----- Session set-up -----

# Load necessary libraries
library(caret)
library(ranger)
library(coin)
library(rsample)
library(doParallel)
library(forcats)
library(glue)
library(gtable)
library(furrr)
library(foreach)
library(tidyverse)

# Export plots and data?
export = TRUE

# Parallize? 32 cores if on threadripper - 8 if local
parallel = TRUE 

# Check which device is running
node_name <- Sys.info()["nodename"]

# Set file paths conditionally
path_out <-  "~/Git/merlin/traits_output"

## ---------- Fit Stratified RF Models and Calculate PDP ----------

# Read input data
data <- read_csv("~/Git/merlin/traits_merlin.csv")   

# Split into training and test sets
split <- initial_split(data, prop = 0.8)

# Extract training and testing datasets
train_data <- training(split)
test_data  <- testing(split)

# Define traits
traits <- c("wood_density", "bark_thickness", "conduit_diam", "leaf_n", 
            "specific_leaf_area", "seed_dry_mass", "shade_tolerance", "height")

# Define covariates 
covariates <- c("standage", "temp_pc", "soil_pc", "rain_pc", "elevation", "soil_ph", 
                "biome_boreal_forests_or_taiga", "biome_flooded_grasslands", 
                "biome_mediterranean_woodlands", "biome_temperate_broadleaf_forests",
                "biome_temperate_conifer_forests", "biome_temperate_grasslands",
                "biome_tundra", "biome_xeric_shrublands")

# Hyperparameter grid with reduced num.trees to speed things up 
hyper_grid <- tibble(
    trait = traits,
    num.trees = sample(c(100, 150, 200), 1),  # Randomly choose trees between 100-200
    mtry = rep(4, length(traits)),
    min.node.size = rep(1, length(traits)))

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

# Set Quantile values; covariates, num_threads; and quantile labels
quantiles_25 <- c(0.25, 0.75)
num_threads <- ifelse(node_name == "threadeast", 4, 1)

quantile_labels <- list(
    temp_pc = c("Cold Temperatures (Lower 25%)", "Warm Temperatures (Upper 25%)"),
    soil_pc = c("Sandy Soils (Lower 25%)", "Water-Retentive Soils (Upper 25%)"),
    rain_pc = c("Low Precipitation (Lower 25%)", "High Precipitation (Upper 25%)"),
    elevation = c("Low Elevation (Lower 25%)", "High Elevation (Upper 25%)"),
    soil_ph = c("Low Soil pH (Lower 25%)", "High Soil pH (Upper 25%)"))

# Set up parallel backend
num_bootstrap <- 100  # Number of bootstrap iterations
num_threads <- ifelse(node_name == "threadeast", 10, 1)

# set cluster 
if (parallel) { 
    num_cores <-  ifelse(node_name == "threadeast", 4, 8)
    cl <- makeCluster(num_cores)
    registerDoParallel(cl, cores = num_cores)
    getDoParWorkers()
}

# Ensure all functions and variables are available in each worker
clusterExport(cl, c("fit_rf_model", "fit_models_on_strata", "calculate_pdp_for_scenario", 
                    "calculate_partial_dependence", "stratify", "traits", "covariates", 
                    "hyper_grid", "quantile_labels"))

# Bootstrap loop
bootstrap_results <- foreach(iter = 1:num_bootstrap, .combine = bind_rows, .packages = c("ranger", "foreach", "doParallel", "rsample", "tidyverse")) %dopar% {
    
    # Draw bootstrap sample (60-80% of the data)
    boot_data <- data %>% sample_frac(runif(1, 0.6, 0.8), replace = TRUE)
    
    # Fit models on stratified bootstrap data
    temp <- fit_models_on_strata(boot_data, "temp_pc", traits, covariates, hyper_grid, num_threads, quantile_labels$temp_pc)
    soil <- fit_models_on_strata(boot_data, "soil_pc", traits, covariates, hyper_grid, num_threads, quantile_labels$soil_pc)
    prcp <- fit_models_on_strata(boot_data, "rain_pc", traits, covariates, hyper_grid, num_threads, quantile_labels$rain_pc)
    elev <- fit_models_on_strata(boot_data, "elevation", traits, covariates, hyper_grid, num_threads, quantile_labels$elevation)
    ph   <- fit_models_on_strata(boot_data, "soil_ph", traits, covariates, hyper_grid, num_threads, quantile_labels$soil_ph)
    
    # Extract R² values for each model
    performance <- bind_rows(temp$performance, soil$performance, prcp$performance,
                             elev$performance, ph$performance) %>%
        dplyr::select(trait, group, rsq, pred_error) %>%
        mutate(iteration = iter)  # Add iteration ID
    
    # Compute Partial Dependence for each trait and scenario
    pdp_temp <- calculate_pdp_for_scenario(temp, traits, "standage", quantile_labels$temp_pc) %>% mutate(variable = "Temperature (PC)")
    pdp_soil <- calculate_pdp_for_scenario(soil, traits, "standage", quantile_labels$soil_pc) %>% mutate(variable = "Soil - Water Retention (PC)")
    pdp_prcp <- calculate_pdp_for_scenario(prcp, traits, "standage", quantile_labels$rain_pc) %>% mutate(variable = "Precipitation (PC)")
    pdp_elev <- calculate_pdp_for_scenario(elev, traits, "standage", quantile_labels$elevation) %>% mutate(variable = "Elevation")
    pdp_ph   <- calculate_pdp_for_scenario(ph, traits, "standage", quantile_labels$soil_ph) %>% mutate(variable = "Soil pH")
    
    # Combine all PDP results
    pdp_data <- bind_rows(pdp_temp, pdp_soil, pdp_prcp, pdp_elev, pdp_ph) %>%
        left_join(performance, by = c("trait", "group")) %>%
        mutate(iteration = iter)  # Add iteration ID
    
    return(pdp_data)
}

stopImplicitCluster()

if (export) {
    write_csv(bootstrap_results, file = paste0(path_out, "/bootstrap_pdp.csv"))
}

### END 

