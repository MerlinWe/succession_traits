
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
	select(-Name) 

## ----- Preliminary analysis of shade tolerance ----------

# Get input data for fitting a random forest 
shade <- fia_clean %>% 
	filter(standage < 500) %>%
	
	select(wmean_shade_tolerance, standage, ecoregion, managed, 
				 annual_mean_temperature, annual_precipitation, temperature_seasonality, mean_diurnal_range, 
				 min_temperature_of_coldest_month, max_temperature_of_warmest_month,
				 elevation, pop_density, sand_content_015cm, soil_ph_015cm, water_capacity_015cm) %>%
	
	filter(complete.cases(.)) %>%
	as_tibble() 

ggplot(shade, aes(x = wmean_shade_tolerance)) +
	geom_histogram(binwidth = 0.05, fill = 'blue', color = 'black', alpha = .8) +
	theme_bw() +
	labs(x = 'Shade tolerance', y = 'Frequency')

ggplot(shade, aes(x = standage, y = wmean_shade_tolerance)) +
	geom_hex(alpha = .8) +
	geom_smooth(method = "loess", se = FALSE) +
	theme_bw() +
	labs(x = 'Standage', y = 'Shade Tolerance')

## Random forest regression workflow 
set.seed(42)

# withhold 10% for final testing; get the test ids at random
test_ids <- sample(1:nrow(shade), round(.10*nrow(shade)))
test_data <- shade %>% slice(test_ids)

# the rest of the data go into testing/training
train_data <- shade %>% slice(-test_ids)

## ----- 1. Initial Model Fit (as a Baseline) -----

# Initial model fit with default parameters
baseline_rf <- ranger(
	# RF formula 
	wmean_shade_tolerance ~ standage + 
		managed + annual_mean_temperature + annual_precipitation + temperature_seasonality +
		mean_diurnal_range + min_temperature_of_coldest_month + max_temperature_of_warmest_month +
		elevation + pop_density + sand_content_015cm + soil_ph_015cm + water_capacity_015cm,
	# Data and variable importance
	data = train_data,
	importance = "impurity")

## SHAPLEY stuff 

library(colorspace)
library(doParallel)
library(foreach)
library(viridis)
library(latex2exp)
library(cowplot)
library(fastshap)
library(ggbeeswarm)

# specify the number of out-of-fit test points, and number of simulations for fastshap. the more sims the better
ntest <- 5000
nsim <- 100

# register parallel cores
registerDoParallel(min(detectCores()-2, 72))

# prediction function for fastshap
pfun <- function(object, newdata) {
	predict(object, data = newdata)$predictions
}

set.seed(15)

# create the test/train split
test_id <- treg_use %>% slice(sample(1:nrow(.), ntest, replace = FALSE)) %>% select(id) %>% unlist() %>% as.numeric()
train_data <- treg_use %>% filter(!id%in%test_id)
test_data <- treg_use %>% filter(id%in%test_id) 

# Fit the random forest models to the training data
r1 <- ranger(PC1~., data = train_data %>% select(PC1, contains("Env_")), case.weights = train_data$wt)
r2 <- ranger(PC2~., data = train_data %>% select(PC2, contains("Env_")), case.weights = train_data$wt)

#  Estimate the shapley values on the test data
fs1 <- fastshap::explain(r1, X = test_data %>% select(contains("Env_")) %>% as.matrix(), pred_wrapper = pfun, nsim = nsim, .parallel = TRUE, adjust = TRUE)
fs2 <- fastshap::explain(r2, X = test_data %>% select(contains("Env_")) %>% as.matrix(), pred_wrapper = pfun, nsim = nsim, .parallel = TRUE, adjust = TRUE)

# combine the models and merge in the environmental covariates
shap_df <- bind_rows(fs1 %>% as_tibble() %>% mutate(id = test_data$id) %>% gather(var, shap, -id) %>% mutate(type = "shap", axis = 1),
										 fs2 %>% as_tibble() %>% mutate(id = test_data$id) %>% gather(var, shap, -id) %>% mutate(type = "shap", axis = 2)) %>%
	left_join(test_data %>% select(LAT, LON, accepted_bin, angio, id, contains("Env_")) %>% gather(var, value, -id, -angio, -accepted_bin, -LAT, -LON)) 


