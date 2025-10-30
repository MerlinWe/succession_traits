## OLD 

#
Define function to stratify data into high/low environments and calculate MSE
compute_mse_by_bin <- function(trait, model, data, env_vars) {
	
	data$prediction <- predict(model, data)$predictions
	data$residual_sq <- (data[[trait]] - data$prediction)^2
	
	# Create standage bins internally
	data <- data %>%
		mutate(standage_bin = cut(standage, breaks = seq(0, 150, by = 10),
															include.lowest = TRUE, right = FALSE))
	
	# Loop over environmental variables
	map_dfr(env_vars, function(env) {
		
		# Stratify into quantiles
		env_qs <- quantile(data[[env]], probs = c(0.25, 0.75), na.rm = TRUE)
		data_strat <- data %>%
			mutate(env_group = case_when(
				!!sym(env) <= env_qs[1] ~ "low",
				!!sym(env) >= env_qs[2] ~ "high",
				TRUE ~ NA_character_
			)) %>%
			filter(!is.na(env_group))
		
		# Summarize MSE
		data_strat %>%
			group_by(trait = trait, variable = env, env_group, standage_bin) %>%
			summarise(mse = mean(residual_sq, na.rm = TRUE), .groups = "drop")
	})
}

# NEW: 



env_vars <- c("elevation", "rain_pc", "soil_pc", "soil_ph", "temp_pc")



ve_summaries <- map_dfr(traits, ~ oof_preds_by_bins(.x, data, env_vars,
                                                    standage_breaks = seq(0,150,10),
                                                    v = 10, repeats = 100))

# Aggregate over repeats, keep uncertainty
ve_final <- ve_summaries %>%
  group_by(trait, variable, env_group, standage_bin) %>%
  summarise(
    n_med = median(n),
    VEcv_med = median(VEcv, na.rm = TRUE),
    VEcv_lwr = quantile(VEcv, 0.025, na.rm = TRUE),
    VEcv_upr = quantile(VEcv, 0.975, na.rm = TRUE),
    E1_med   = median(E1, na.rm = TRUE),
    E1_lwr   = quantile(E1, 0.025, na.rm = TRUE),
    E1_upr   = quantile(E1, 0.975, na.rm = TRUE),
    .groups = "drop"
  )



### 


library(dplyr)
library(purrr)
library(rsample)
library(ranger)
library(stringr)

# VEcv: variance-explained by cross-validated predictions
VEcv <- function(obs, pred) {
  1 - mean((obs - pred)^2, na.rm = TRUE) / var(obs, na.rm = TRUE)
}

# Legates-McCabe E1: linear-error skill
E1 <- function(obs, pred) {
  1 - (sum(abs(obs - pred), na.rm = TRUE) / sum(abs(obs - mean(obs, na.rm = TRUE)), na.rm = TRUE))
}

# Make buffered folds 
make_folds <- function(data, v = 10) {
  vfold_cv(data, v = v)  # rsample
}

oof_skill_by_bins <- function(trait, data, env_vars,
                              standage_breaks = seq(0, 150, 10),
                              v = 10, repeats = 50,
                              rf_args = list(num.trees = 1000, min.node.size = 5, num.threads = parallel::detectCores())) {
  
  # Prepare bins
  data <- data %>%
    mutate(standage_bin = cut(standage, breaks = standage_breaks, include.lowest = TRUE, right = FALSE))
  
  reps <- seq_len(repeats)
  
  map_dfr(reps, function(rp) {
    folds <- make_folds(data, v = v)
    preds_all <- vector("list", length(folds$splits))
    
    for (i in seq_along(folds$splits)) {
      train <- analysis(folds$splits[[i]])
      test  <- assessment(folds$splits[[i]])
      
      predictors <- setdiff(names(data), traits)  # or explicitly your predictor set
      form <- reformulate(termlabels = predictors, response = trait)
      
      rf_fit <- do.call(ranger::ranger, c(list(
        formula = form,
        data    = train %>% select(all_of(c(trait, predictors))) %>% drop_na()
      ), rf_args))
      
      test$prediction <- predict(rf_fit, data = test)$predictions
      preds_all[[i]] <- test
    }
    
    oof <- bind_rows(preds_all)
    
    # Stratify per env variable
    map_dfr(env_vars, function(env) {
      qs <- quantile(oof[[env]], c(0.25, 0.75), na.rm = TRUE)
      oof_strat <- oof %>%
        mutate(env_group = case_when(
          .data[[env]] <= qs[1] ~ "low",
          .data[[env]] >= qs[2] ~ "high",
          TRUE ~ NA_character_
        ),
        variable = env) %>%
        filter(!is.na(env_group))
      
      # Summarise skill per bin
      oof_strat %>%
        group_by(trait = trait, variable, env_group, standage_bin) %>%
        summarise(
          VEcv = VEcv(.data[[trait]], prediction),
          E1   = E1(.data[[trait]], prediction),
          n    = n(),
          .groups = "drop"
        ) %>%
        mutate(repeat_id = rp)
    })
  })
}

env_vars <- c("elevation","rain_pc","soil_pc","soil_ph","temp_pc")

ve_results <- map_dfr(traits, ~ oof_skill_by_bins(.x, data, env_vars, repeats = 100))
