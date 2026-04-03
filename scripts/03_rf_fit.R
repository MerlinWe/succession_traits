################################################################################
## succession_traits: 03 — Random forest model fitting

## Fit one random forest model per trait × leaf type, with optional
## hyperparameter tuning via OOB grid search. Evaluates performance on a
## held-out test set and saves model objects for downstream SHAP and PDP

################################################################################

rm(list = ls())
set.seed(42)

# ── Libraries ─────────────────────────────────────────────────────────────────
library(ranger)
library(doParallel)
library(foreach)
library(tidyverse)

# ── Configuration ─────────────────────────────────────────────────────────────
TUNING      <- TRUE   # TRUE: run grid search; FALSE: load saved hyperparameters
PARALLEL    <- TRUE
N_CORES     <- parallel::detectCores(logical = FALSE) - 1L  # leave one core free
N_TEST_MAX  <- 5000L   # cap on test set size per leaf type

PATH_IN         <- "data_processed/fia_traits_clean.rds"
PATH_MODELS     <- "models"
PATH_TABLES     <- "tables"
PATH_PLOTS      <- "plots/tuning"

# ── Functions ─────────────────────────────────────────────────────────────────
source("scripts/functions.R")

# ── Parallel backend ──────────────────────────────────────────────────────────
if (PARALLEL) {
	cl <- makeCluster(N_CORES)
	registerDoParallel(cl)
	message(sprintf("Parallel backend: %d cores", N_CORES))
} else {
	registerDoSEQ()
	message("Running sequentially")
}

# ══════════════════════════════════════════════════════════════════════════════
# 1. Load data and define modelling vocabulary
# ══════════════════════════════════════════════════════════════════════════════

data <- read_rds(PATH_IN)

# Derive leaf type from biome dummies
# Only broadleaf and coniferous forest plots enter the analysis

data <- data %>%
	mutate(
		leaf_type = case_when(
			biome_boreal_forests_or_taiga   == 1 |
				biome_temperate_conifer_forests == 1 ~ "coniferous",
			biome_temperate_broadleaf_forests == 1 |
				biome_mediterranean_woodlands   == 1 ~ "broadleaf",
			TRUE ~ NA_character_
		)
	) %>%
	filter(!is.na(leaf_type))

message(sprintf(
	"Modelling dataset: %d plots  (%d broadleaf / %d coniferous)",
	nrow(data),
	sum(data$leaf_type == "broadleaf"),
	sum(data$leaf_type == "coniferous")
))

# Traits to model (response variables)
TRAIT_LABELS <- c(
	"bark_thickness"     = "Bark Thickness",
	"conduit_diam"       = "Conduit Diameter",
	"height"             = "Tree Height",
	"leaf_density"       = "Leaf Density",
	"leaf_k"             = "Leaf Potassium",
	"root_depth"         = "Root Depth",
	"seed_dry_mass"      = "Seed Dry Mass",
	"shade_tolerance"    = "Shade Tolerance",
	"specific_leaf_area" = "Specific Leaf Area"
)
TRAITS <- names(TRAIT_LABELS)

# Predictors (covariates)
COVARIATES <- c("standage", "temp_pc", "soil_pc", "rain_pc", "elevation", "soil_ph")

# Leaf types to loop over
LEAF_TYPES <- c("broadleaf", "coniferous")


# ══════════════════════════════════════════════════════════════════════════════
# 2. Hyperparameter grid
# ══════════════════════════════════════════════════════════════════════════════

if (TUNING) {
	# Grid to search during tuning phase
	BASE_GRID <- expand.grid(
		num.trees     = c(500L, 1000L, 1500L),
		mtry          = 2:4,
		min.node.size = c(1L, 10L, 20L),
		stringsAsFactors = FALSE
	)
	message(sprintf("Tuning grid: %d combinations per trait", nrow(BASE_GRID)))
} else {
	# Load previously tuned hyperparameters
	# Expected columns: trait, num_trees, mtry, min_node_size
	tuned_params <- list(
		broadleaf  = read_csv(file.path(PATH_TABLES, "perf_broadleaf.csv"),
													show_col_types = FALSE),
		coniferous = read_csv(file.path(PATH_TABLES, "perf_coniferous.csv"),
													show_col_types = FALSE)
	)
	message("Using pre-tuned hyperparameters")
}


# ══════════════════════════════════════════════════════════════════════════════
# 3. Fit models per leaf type
# ══════════════════════════════════════════════════════════════════════════════

perf_all <- list()

for (lt in LEAF_TYPES) {
	
	message(sprintf("\n══ Leaf type: %s ══════════════════════════════════════════", lt))
	
	# ── 3.1 Subset and split ────────────────────────────────────────────────────
	df_lt <- data %>%
		filter(leaf_type == lt) %>%
		mutate(.id = row_number())
	
	if (nrow(df_lt) < 1000L)
		warning(sprintf("Only %d rows for '%s' — results may be unstable.", nrow(df_lt), lt))
	
	n_test  <- min(N_TEST_MAX, floor(0.2 * nrow(df_lt)))
	splits  <- split_train_test(df_lt, ntest = n_test, seed = 42L)
	
	train_data <- splits$train
	test_data  <- splits$test
	
	message(sprintf(
		"  Train: %d plots  |  Test: %d plots",
		nrow(train_data), nrow(test_data)
	))
	
	# Save splits for downstream scripts (SHAP runs on test set)
	saveRDS(
		list(train = train_data, test = test_data),
		file = file.path(PATH_MODELS, sprintf("splits_%s.rds", lt))
	)
	
	# ── 3.2 Hyperparameter tuning (optional) ────────────────────────────────────
	if (TUNING) {
		message("  Tuning hyperparameters (OOB grid search)...")
		tuned_list <- vector("list", length(TRAITS))
		names(tuned_list) <- TRAITS
		
		for (tr in TRAITS) {
			message(sprintf("    · %s", tr))
			tune_res <- tune_rf_model(
				trait      = tr,
				data       = train_data,
				covariates = COVARIATES,
				hyper_grid = BASE_GRID
			)
			tuned_list[[tr]] <- tune_res$best_hyperparameters
			
			ggsave(
				filename = file.path(PATH_PLOTS, sprintf("tuning_%s_%s.png", lt, tr)),
				plot     = tune_res$error_plot,
				width = 4, height = 3, dpi = 300
			)
		}
		
		hyper_grid <- bind_rows(tuned_list)
		
		hyper_grid <- bind_rows(tuned_list) %>%
			rename(num_trees     = num.trees,
						 min_node_size = min.node.size)
		
	} else {
		hyper_grid <- tuned_params[[lt]] %>%
			dplyr::select(trait, num_trees, mtry, min_node_size)
	}
	
	# ── 3.3 Fit models ───────────────────────────────────────────────────────────
	message("  Fitting models...")
	
	# Parallelise across traits; use num_threads = 1 inside ranger to avoid
	# oversubscription (outer foreach already uses all available cores)
	fitted <- foreach(
		tr        = TRAITS,
		.packages = c("ranger", "dplyr", "tibble"),
		.export   = c("fit_rf_model")
	) %dopar% {
		fit_rf_model(
			trait            = tr,
			df_train         = train_data,
			covariates       = COVARIATES,
			hyper_parameters = hyper_grid,
			num_threads      = 1L          # <-- 1 thread per worker; workers use cores
		)
	}
	names(fitted) <- TRAITS
	
	# ── 3.4 Evaluate on held-out test set ────────────────────────────────────────
	# OOB R² (from ranger) is an in-sample estimate and can be optimistic.
	# We report test-set R² and RMSE as the primary performance metrics.
	
	message("  Evaluating on test set...")
	
	perf_test <- map_dfr(TRAITS, function(tr) {
		mod   <- fitted[[tr]]$trait_mod
		preds <- predict(mod, data = test_data)$predictions
		obs   <- test_data[[tr]]
		
		# Remove NA pairs (shade_tolerance may have missing values)
		complete <- !is.na(obs) & !is.na(preds)
		obs   <- obs[complete]
		preds <- preds[complete]
		
		ss_res <- sum((obs - preds)^2)
		ss_tot <- sum((obs - mean(obs))^2)
		
		tibble(
			trait         = tr,
			leaf_type     = lt,
			n_train       = nrow(train_data),
			n_test        = sum(complete),
			# Test-set metrics (primary)
			r2_test       = 1 - ss_res / ss_tot,
			rmse_test     = sqrt(mean((obs - preds)^2)),
			# OOB metrics (in-sample reference)
			r2_oob        = fitted[[tr]]$performance$rsq,
			rmse_oob      = sqrt(fitted[[tr]]$performance$pred_error),
			# Hyperparameters used
			num_trees     = fitted[[tr]]$performance$num_trees,
			mtry          = fitted[[tr]]$performance$mtry,
			min_node_size = fitted[[tr]]$performance$min_node_size
		)
	})
	
	perf_all[[lt]] <- perf_test
	
	# ── 3.5 Save models and performance ─────────────────────────────────────────
	walk(TRAITS, function(tr) {
		saveRDS(
			fitted[[tr]]$trait_mod,
			file = file.path(PATH_MODELS, sprintf("rf_%s_%s.rds", tr, lt))
		)
	})
	
	write_csv(
		perf_test,
		file.path(PATH_TABLES, sprintf("perf_%s.csv", lt))
	)
	
	message(sprintf("  ✓ %s complete", lt))
}

# Stop cluster
if (PARALLEL) stopCluster(cl)


# ══════════════════════════════════════════════════════════════════════════════
# 4. Performance summary
# ══════════════════════════════════════════════════════════════════════════════

perf_combined <- bind_rows(perf_all)

message("\n── Model performance (test-set R²) ────────────────────────────────")
perf_combined %>%
	dplyr::select(trait, leaf_type, r2_test, rmse_test, r2_oob) %>%
	mutate(across(where(is.numeric), ~ round(., 3))) %>%
	arrange(leaf_type, desc(r2_test)) %>%
	print(n = Inf)

# Flag any poorly performing models
poor <- perf_combined %>% filter(r2_test < 0.4)
if (nrow(poor) > 0) {
	warning(sprintf(
		"%d model(s) have test R² < 0.4 — review before proceeding:",
		nrow(poor)
	))
	poor %>%
		dplyr::select(trait, leaf_type, r2_test) %>%
		pwalk(~ message(sprintf("  ⚠ %s (%s): R² = %.3f", ..1, ..2, ..3)))
}

# Save combined performance table
write_csv(
	perf_combined,
	file.path(PATH_TABLES, "perf_combined.csv")
)

message(sprintf(
	"\nModels saved to: %s/rf_<trait>_<leaftype>.rds",
	PATH_MODELS
))
message(sprintf(
	"Performance saved to: %s/perf_combined.csv",
	PATH_TABLES
))


