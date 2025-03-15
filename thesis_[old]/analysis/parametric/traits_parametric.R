############################################################################################################################
########################################  MSc Diss. Forest Succession Data Analysis ########################################  
############################################################################################################################

rm(list = ls())   # make sure environment is clean 
set.seed(42)      # set seed for reproducibility

# ----- Session set-up -----

# Load necessary libraries
library(nlme)
library(car)
library(MuMIn)
library(emmeans)
library(mgcv)
library(gstat)
library(ggeffects)
library(cowplot)
library(gridExtra)
library(psych)
library(spdep)
library(gstat)
library(sp)
library(corrplot)
library(factoextra)
library(doParallel)
library(foreach)
library(tidyverse)

export = FALSE # export? 

# Check which device is running
node_name <- Sys.info()["nodename"]

# Set file paths conditionally
path <- ifelse(node_name == "threadeast", "/home/merlin/traits_output/parametric", 
									 "/Users/serpent/Documents/MSc/Thesis/Code/analysis/parametric")

parallel = TRUE # parallel?

if (parallel) { # set cluster 
	num_cores <-  ifelse(node_name == "threadeast", 32, 8)
	cl <- makeCluster(num_cores)
	registerDoParallel(cl, cores = num_cores)
	getDoParWorkers()
}

dat <- read_csv(paste0(path, "/traits_parametric.csv")) %>%
	
	# Make sure biome is a factor 
	mutate(biome = as.factor(biome)) %>%
	
	# Add a small jitter to coords
	mutate(
		LAT_jitter = LAT + runif(n(), -0.0001, 0.0001),
		LON_jitter = LON + runif(n(), -0.0001, 0.0001))


# *******************************************************************************************************************
# *************************** Full Factorial Approach using Interactions in Global Models ***************************
# *******************************************************************************************************************

# Define the analysis with the interaction term 
analyze_trait <- function(dat, trait) {
	
	cat("\nAnalyzing trait:", trait, "\n")
	
	# --- Fit the LME Model ---
	model <- lme(
		as.formula(paste(trait, "~ standage * (temp_pc + rain_pc + soil_pc + elevation + soil_ph) + biome")),
		random = ~ 1 | PID,
		correlation = corExp(form = ~ LAT_jitter + LON_jitter, nugget = TRUE),
		control = lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100),
		data = dat
	)
	
	# Print model summary
	cat("\nModel Summary:\n")
	print(summary(model))
	
	# --- Extract Residuals and Fitted Values ---
	
	residuals <- residuals(model, type = "normalized")
	fitted_values <- fitted(model)
	
	# --- Diagnostic Plots ---
	
	# Q-Q Plot
	p1 <- ggplot(data.frame(residuals), aes(sample = residuals)) +
		stat_qq() +
		stat_qq_line(color = "red") +
		ggtitle(paste("Q-Q Plot for", trait)) +
		theme_minimal()
	
	# Residuals vs Fitted Plot
	p2 <- ggplot(data.frame(fitted = fitted_values, residuals = residuals),
							 aes(x = fitted, y = residuals)) +
		geom_point(alpha = 0.5) +
		geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
		ggtitle(paste("Residuals vs Fitted for", trait)) +
		xlab("Fitted Values") +
		ylab("Residuals") +
		theme_minimal()
	
	# Combine diagnostic plots
	grid.arrange(p1, p2, ncol = 2)
	
	# --- Moran's I Test for Spatial Autocorrelation ---
	
	cat("\nPerforming Moran's I test...\n")
	coordinates <- cbind(dat$LAT_jitter, dat$LON_jitter)
	
	# Handle disconnected points by creating a connected neighbours list
	suppressWarnings({
		nb <- knn2nb(knearneigh(coordinates, k = 5))
		if (any(card(nb) == 0)) {
			cat("Warning: Isolated points detected. Connecting all points using a minimum distance threshold.\n")
			dnb <- dnearneigh(coordinates, 0, max(dist(coordinates)))
			nb <- dnb
		}
	})
	
	listw <- nb2listw(nb, style = "W", zero.policy = TRUE)  # Handle zero neighbours
	moran_test <- moran.test(residuals, listw, zero.policy = TRUE)
	print(moran_test)
	
	# --- Semivariogram ---
	
	cat("\nGenerating Semivariogram...\n")
	dat_copy <- dat %>%
		mutate(residuals = residuals) 
	coordinates(dat_copy) <- ~ LAT_jitter + LON_jitter
	
	variogram_res <- variogram(residuals ~ 1, dat_copy)
	plot(variogram_res, main = paste("Semivariogram of Residuals for", trait))
	
	# --- Extract Coefficients for Interactions ---
	
	cat("\nExtracting coefficients for interactions...\n")
	coef_summary <- summary(model)$tTable
	
	# Create a data frame of all coefficients
	coefficients_df <- data.frame(
		Term = rownames(coef_summary),
		Estimate = coef_summary[, "Value"],
		StdError = coef_summary[, "Std.Error"],
		pValue = coef_summary[, "p-value"]
	)
	
	# Add a column to distinguish interactions and main effects
	coefficients_df$EffectType <- ifelse(grepl(":", coefficients_df$Term), "Interaction", "Main Effect")
	
	# Print the coefficients
	print(coefficients_df)
	
	# Return all results as a list
	return(list(
		model = model,
		residuals = residuals,
		moran_test = moran_test,
		variogram = variogram_res,
		coefficients = coefficients_df
	))
}

# Define traits
traits <- c("wood_density", "bark_thickness", "conduit_diam", "leaf_n", "specific_leaf_area", "seed_dry_mass", "shade_tolerance", "height")

# Initialize a list to store results
global_models <- list()

# Loop through traits in parallel 
global_models <- foreach(trait = traits, .packages = c("nlme", "car", "MuMIn", "emmeans", "mgcv", "gstat",
																											 "cowplot", "gridExtra", "psych", "spdep", "sp", "tidyverse")) %dopar% {
	
	tryCatch({
		cat("\nFitting global model for:", trait, "\n")
		analyze_trait(dat, trait)
	}, error = function(e) {
		cat("\nError for trait:", trait, "\n", conditionMessage(e), "\n")
		return(NULL)  
	})
}

# Stop the cluster after computation
stopCluster(cl)

# Check Moran's I significance for spatial autocorrelation
check_autocorrelation <- function(model) {
	if (!is.null(model) && !is.null(model$moran_test)) {
		p_value <- model$moran_test$p.value
		if (p_value < 0.05) {
			return("Yes")
		} else {
			return("No")
		}
	} else {
		return("No Model/No Moran Test")
	}
}

# Name the list of results and check spatial autocorrelation
names(global_models) <- traits
sapply(global_models, check_autocorrelation)

# Extract coefficients and add the trait name
coefficients_combined <- do.call(rbind, lapply(names(global_models), function(trait) {
	model <- global_models[[trait]]
	if (!is.null(model) && !is.null(model$coefficients)) {
		coeffs <- model$coefficients
		coeffs$Trait <- trait
		return(coeffs)
	} else {
		return(NULL)
	}
}))

# Convert coefficients to a large tibble
coefficients_combined <- as_tibble(coefficients_combined) %>%
	mutate(
		EffectType = ifelse(grepl(":", Term), "Interaction", "Main Effect"))

# Plot main effects
global_main_effects <- coefficients_combined %>% 
	filter(EffectType == "Main Effect") %>%
	filter(!startsWith(Term, "biome")) %>% # exclude biomes? 
	filter(!startsWith(Term, "(Intercept)")) %>% # exclude intercept?
	ggplot(aes(x = Trait, y = Estimate, fill = Trait)) +
	geom_bar(stat = "identity", colour = "black", alpha = .8) +
	geom_errorbar(aes(ymin = Estimate - 1.96 * StdError, ymax = Estimate + 1.96 * StdError), 
								width = 0.2, color = "black") +
	facet_wrap(~ Term, scales = "free", ncol = 5, nrow = 3) +
	scale_fill_viridis_d() +
	theme_bw() +
	labs(title = "Main Effects", x = NULL, y = "Estimate") +
	theme(axis.text.x = element_text(angle = 45, hjust = 1),
				legend.position = "none",
				text = element_text(size = 8))

# Plot interactions
global_interactions <- coefficients_combined %>% 
	filter(EffectType == "Interaction") %>%
	ggplot(aes(x = Trait, y = Estimate, fill = Term)) +
	geom_bar(stat = "identity", colour = "black", alpha = .8) +
	geom_errorbar(aes(ymin = Estimate - 1.96 * StdError, ymax = Estimate + 1.96 * StdError), 
								width = 0.2, color = "black") +
	facet_wrap(~ Term, scales = "free", ncol = 3, nrow = 2) +
	scale_fill_viridis_d() +
	theme_bw() +
	labs(title = "Interaction", x = NULL, y = "Estimate") +
	theme(axis.text.x = element_text(angle = 45, hjust = 1),
				legend.position = "none",
				text = element_text(size = 8))

if (export) {
	
	ggsave(paste0(path, "/global_main_effects.png"),
				 plot = global_main_effects,
				 bg = "white",
				 width = 250,  
				 height = 200, 
				 units = "mm",
				 dpi = 600)
	
	ggsave(paste0(path, "/global_interactions.png"),
				 plot = global_interactions,
				 bg = "white",
				 width = 250,  
				 height = 200, 
				 units = "mm",
				 dpi = 600)
}

# **********************************************************************************************************
# *************************** Additive Approach using Stratified Scenario Models ***************************
# **********************************************************************************************************

# Define analysis without interaction term 
analyze_trait <- function(dat, trait) {
	cat("\nAnalyzing trait:", trait, "\n")
	
	# --- Fit the LME Model (Main Effects Only) ---
	model <- lme(
		as.formula(paste(trait, "~ standage + temp_pc + rain_pc + soil_pc + elevation + soil_ph + biome")),
		random = ~ 1 | PID,
		correlation = corExp(form = ~ LAT_jitter + LON_jitter, nugget = TRUE),
		control = lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100),
		data = dat
	)
	
	# Print model summary
	cat("\nModel Summary:\n")
	print(summary(model))
	
	# --- Extract Residuals and Fitted Values ---
	residuals <- residuals(model, type = "normalized")
	fitted_values <- fitted(model)
	
	# --- Diagnostic Plots ---
	p1 <- ggplot(data.frame(residuals), aes(sample = residuals)) +
		stat_qq() +
		stat_qq_line(color = "red") +
		ggtitle(paste("Q-Q Plot for", trait)) +
		theme_minimal()
	
	p2 <- ggplot(data.frame(fitted = fitted_values, residuals = residuals),
							 aes(x = fitted, y = residuals)) +
		geom_point(alpha = 0.5) +
		geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
		ggtitle(paste("Residuals vs Fitted for", trait)) +
		xlab("Fitted Values") +
		ylab("Residuals") +
		theme_minimal()
	
	grid.arrange(p1, p2, ncol = 2)
	
	# --- Return Model ---
	return(list(model = model, residuals = residuals))
}

# Define stratification 
stratify <- function(data, variable, quantiles) {
	
	lower_quantile <- quantile(data[[variable]], quantiles[1])
	upper_quantile <- quantile(data[[variable]], quantiles[2])
	
	lower_data <- data %>% filter(.[[variable]] <= lower_quantile)
	upper_data <- data %>% filter(.[[variable]] >= upper_quantile)
	
	return(list(lower = lower_data, upper = upper_data))
}

# Helper function to extract slope and intercept 
extract_intercept_slope <- function(model, trait, group, variable) {
	if (is.null(model) || is.null(model$model)) {
		return(tibble(
			Trait = trait,
			Group = group,
			Variable = variable,
			Intercept = NA_real_,
			Slope = NA_real_
		))
	}
	
	# Extract coefficients safely
	coef_summary <- tryCatch(
		summary(model$model)$tTable,
		error = function(e) return(NULL)
	)
	
	# Check if coef_summary is valid
	if (is.null(coef_summary)) {
		return(tibble(
			Trait = trait,
			Group = group,
			Variable = variable,
			Intercept = NA_real_,
			Slope = NA_real_
		))
	}
	
	# Safely retrieve intercept and slope
	intercept <- if ("(Intercept)" %in% rownames(coef_summary)) coef_summary["(Intercept)", "Value"] else NA_real_
	slope <- if ("standage" %in% rownames(coef_summary)) coef_summary["standage", "Value"] else NA_real_
	
	# Return tibble
	return(tibble(
		Trait = trait,
		Group = group,
		Variable = variable,
		Intercept = intercept,
		Slope = slope
	))
}

# Define Workflow to fit additive scenario models 
fit_mixed_models_on_strata <- function(data, variable, traits, quantiles, labels) {
	# Step 1: Stratify the data by the feature scenario
	stratified_data <- stratify(data, variable, quantiles)
	
	# Step 2: Fit models for the lower quantile group
	cat("\nFitting models for the lower quantile of", variable, "...\n")
	lower_models <- map(traits, ~ tryCatch(
		analyze_trait(stratified_data$lower, .x),
		error = function(e) {
			cat("\nError for trait:", .x, "\n", conditionMessage(e), "\n")
			return(NULL)
		}
	))
	names(lower_models) <- traits
	
	# Step 3: Fit models for the upper quantile group
	cat("\nFitting models for the upper quantile of", variable, "...\n")
	upper_models <- map(traits, ~ tryCatch(
		analyze_trait(stratified_data$upper, .x),
		error = function(e) {
			cat("\nError for trait:", .x, "\n", conditionMessage(e), "\n")
			return(NULL)
		}
	))
	names(upper_models) <- traits
	
	# Step 4: Extract intercepts and slopes for both groups
	lower_coefficients <- map_dfr(lower_models, ~ {
		if (is.null(.x)) {
			return(tibble(
				Trait = "unknown",
				Group = labels[1],
				Variable = variable,
				Intercept = NA_real_,
				Slope = NA_real_
			))
		}
		extract_intercept_slope(.x, .y, labels[1], variable)
	}, .id = "Trait")
	
	upper_coefficients <- map_dfr(upper_models, ~ {
		if (is.null(.x)) {
			return(tibble(
				Trait = "unknown",
				Group = labels[2],
				Variable = variable,
				Intercept = NA_real_,
				Slope = NA_real_
			))
		}
		extract_intercept_slope(.x, .y, labels[2], variable)
	}, .id = "Trait")
	
	# Combine coefficients for both strata
	combined_coefficients <- bind_rows(lower_coefficients, upper_coefficients)
	
	# Return results
	return(list(
		models = list(lower = lower_models, upper = upper_models),
		coefficients = combined_coefficients,
		stratified_data = stratified_data
	))
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


# Fit models stratified by temperature
temp_results <- fit_mixed_models_on_strata(dat, "temp_pc", traits, quantiles_25, quantile_labels$temp_pc)

# Fit models stratified by precipitation
precip_results <- fit_mixed_models_on_strata(dat, "rain_pc", traits, quantiles_25, quantile_labels$rain_pc)

# Fit models stratified by soil characteristics
soil_results <- fit_mixed_models_on_strata(dat, "soil_pc", traits, quantiles_25, quantile_labels$soil_pc)

# Fit models stratified by elevation
elevation_results <- fit_mixed_models_on_strata(dat, "elevation", traits, quantiles_25, quantile_labels$elevation)

# Fit models stratified by soil pH
ph_results <- fit_mixed_models_on_strata(dat, "soil_ph", traits, quantiles_25, quantile_labels$soil_ph)

all_coefficients <- bind_rows(
	temp_results$coefficients %>% mutate(Variable = "Temperature"),
	precip_results$coefficients %>% mutate(Variable = "Precipitation"),
	soil_results$coefficients %>% mutate(Variable = "Soil Characteristics"),
	elevation_results$coefficients %>% mutate(Variable = "Elevation"),
	ph_results$coefficients %>% mutate(Variable = "Soil pH"))

# View combined coefficients
head(all_coefficients)