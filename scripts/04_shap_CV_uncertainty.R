Sys.setenv(SHAP_SMOKE_TEST = "false")

################################################################################
## succession_traits: 04 — SHAP analysis
## Computes Shapley values for each trait × leaf type model on the held-out
## test set, quantifies the relative importance of environmental vs
## successional filtering, and exports results for figure building.
##
## SHAP is computed once on the held-out test set (same split used in
## 03_rf_fit.R). No bootstrapping is applied — with n ≈ 5000 test plots and
## stable RF fits, SHAP sampling variance is negligible. Robustness of the
## environmental vs successional ranking is demonstrated by replication across
## the two independent leaf-type models (broadleaf / coniferous).
##
## Input:  data_processed/fia_traits_clean.rds
##         models/rf_<trait>_<leaftype>.rds     (from 03_rf_fit.R)
##         models/splits_<leaftype>.rds         (from 03_rf_fit.R)
##
## Output: tables/shap_values.rds              (full SHAP long table)
##         tables/shap_importance.rds           (env vs succ ratio summary)
##         tables/shap_per_var.rds              (per-variable importance)
##         tables/shap_importance_ci.rds        (env/succ ratio 95% CI via repeated CV; optional)
##         tables/shap_cv_checkpoints/           (resumable per-job CV checkpoints)
##
## Smoke test (does not overwrite production outputs):
##   Rscript scripts/04_shap_with_CV_uncertainty.R --smoke-test
##
## Author: M. Weiss @ Maynard Lab UCL / ETH Zürich
################################################################################

rm(list = ls())

BASE_SEED <- 42L
RNGkind("L'Ecuyer-CMRG")
set.seed(BASE_SEED)

# ── Libraries ─────────────────────────────────────────────────────────────────
REQUIRED_PACKAGES <- c(
	"fastshap", "ranger", "tidyverse", "doParallel", "foreach"
)
missing_packages <- REQUIRED_PACKAGES[
	!vapply(REQUIRED_PACKAGES, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_packages) > 0L) {
	stop(
		"Missing required R package(s): ", paste(missing_packages, collapse = ", "),
		". Install them before running this script.",
		call. = FALSE
	)
}

library(fastshap)
library(ranger)
library(doParallel)
library(foreach)
library(tidyverse)

# ── Configuration ─────────────────────────────────────────────────────────────
cli_args <- commandArgs(trailingOnly = TRUE)
env_flag <- function(name, default = FALSE) {
	x <- tolower(trimws(Sys.getenv(name, unset = as.character(default))))
	x %in% c("1", "true", "t", "yes", "y", "on")
}

# Activate without editing this file, using either --smoke-test or the
# SHAP_SMOKE_TEST=true environment variable. Smoke outputs have a _smoke suffix
# and use a separate checkpoint directory, so they cannot overwrite a full run.
SMOKE_TEST <- "--smoke-test" %in% cli_args || env_flag("SHAP_SMOKE_TEST")
RUN_SHAP_CI <- !("--skip-ci" %in% cli_args || env_flag("SHAP_SKIP_CI"))

N_SIM           <- if (SMOKE_TEST) 5L   else 100L
MAIN_EXPLAIN_MAX <- if (SMOKE_TEST) 200L else NULL
PARALLEL         <- !SMOKE_TEST

PATH_DATA   <- "data_processed/fia_traits_clean.rds"
PATH_MODELS <- "models"
PATH_TABLES <- "tables"
PATH_PLOTS  <- "figures/supplementary"

OUTPUT_SUFFIX <- if (SMOKE_TEST) "_smoke" else ""
output_rds <- function(stem) {
	file.path(PATH_TABLES, paste0(stem, OUTPUT_SUFFIX, ".rds"))
}

dir.create(PATH_TABLES, recursive = TRUE, showWarnings = FALSE)

FUNCTIONS_FILE <- "scripts/functions.R"
if (!file.exists(FUNCTIONS_FILE)) {
	stop("Required helper file not found: ", FUNCTIONS_FILE, call. = FALSE)
}
source(FUNCTIONS_FILE)
if (!exists("fit_rf_model", mode = "function") ||
		!exists("predict_fn", mode = "function")) {
	stop(
		"scripts/functions.R must define fit_rf_model() and predict_fn().",
		call. = FALSE
	)
}

# ── Vocabulary (must match 03_rf_fit.R) ───────────────────────────────────────

TRAITS <- c("bark_thickness", "conduit_diam", "height", "leaf_density",
						"leaf_k", "root_depth", "seed_dry_mass", "shade_tolerance",
						"specific_leaf_area")

COVARIATES <- c("standage", "temp_pc", "soil_pc", "rain_pc", "elevation", "soil_ph")
LEAF_TYPES <- c("broadleaf", "coniferous")
ENV_VARS  <- c("temp_pc", "soil_pc", "rain_pc", "elevation", "soil_ph")
SUCC_VARS <- "standage"

# Define labels here rather than relying on objects left in the global
# environment or on a plotting script. These mappings are also explicitly
# exported to CV workers below.
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

VAR_LABELS <- c(
	"standage" = "Stand age",
	"temp_pc" = "Temperature PC",
	"soil_pc" = "Soil water retention PC",
	"rain_pc" = "Precipitation PC",
	"elevation" = "Elevation",
	"soil_ph" = "Soil pH"
)

# ── Preflight validation ───────────────────────────────────────────────────────
assert_files_exist <- function(paths, description) {
	missing <- paths[!file.exists(paths)]
	if (length(missing) > 0L) {
		stop(
			description, " missing:\n  ", paste(missing, collapse = "\n  "),
			call. = FALSE
		)
	}
}

assert_columns <- function(x, required, description) {
	missing <- setdiff(required, names(x))
	if (length(missing) > 0L) {
		stop(
			description, " is missing required column(s): ",
			paste(missing, collapse = ", "),
			call. = FALSE
		)
	}
}

model_paths <- unlist(lapply(LEAF_TYPES, function(lt) {
	file.path(PATH_MODELS, sprintf("rf_%s_%s.rds", TRAITS, lt))
}), use.names = FALSE)
split_paths <- file.path(PATH_MODELS, sprintf("splits_%s.rds", LEAF_TYPES))
perf_paths  <- file.path(PATH_TABLES, sprintf("perf_%s.csv", LEAF_TYPES))

assert_files_exist(model_paths, "Random-forest model file(s)")
assert_files_exist(split_paths, "Train/test split file(s)")
if (RUN_SHAP_CI) {
	assert_files_exist(PATH_DATA, "Modelling dataset")
	assert_files_exist(perf_paths, "Tuned-hyperparameter file(s)")
}

if (!setequal(names(TRAIT_LABELS), TRAITS)) {
	stop("TRAIT_LABELS must contain exactly one label for every trait.", call. = FALSE)
}
if (!all(c(COVARIATES, TRAITS) %in% c(names(VAR_LABELS), names(TRAIT_LABELS)))) {
	stop("Label mappings are incomplete.", call. = FALSE)
}

message("\nMode: ", if (SMOKE_TEST) "SMOKE TEST (reduced workload)" else "NORMAL")
message("Preflight: required packages, helper functions, and input files found.")

biome_columns <- c(
	"biome_boreal_forests_or_taiga",
	"biome_temperate_conifer_forests",
	"biome_temperate_broadleaf_forests",
	"biome_mediterranean_woodlands"
)
data_ci_source <- NULL
hyper_ci_source <- NULL
if (RUN_SHAP_CI) {
	data_ci_source <- read_rds(PATH_DATA)
	assert_columns(
		data_ci_source,
		c("PID", TRAITS, COVARIATES, biome_columns),
		"Modelling dataset"
	)
	numeric_columns <- c(TRAITS, COVARIATES)
	non_numeric <- numeric_columns[
		!vapply(data_ci_source[numeric_columns], is.numeric, FUN.VALUE = logical(1))
	]
	if (length(non_numeric) > 0L) {
		stop(
			"These modelling columns must be numeric: ",
			paste(non_numeric, collapse = ", "),
			call. = FALSE
		)
	}
	if (anyNA(data_ci_source$PID)) {
		stop("Modelling dataset column PID contains missing values.", call. = FALSE)
	}

	hyper_ci_source <- map(LEAF_TYPES, function(lt) {
		path <- file.path(PATH_TABLES, sprintf("perf_%s.csv", lt))
		hg <- read_csv(path, show_col_types = FALSE)
		assert_columns(
			hg, c("trait", "num_trees", "mtry", "min_node_size"),
			paste0("Hyperparameter file ", path)
		)
		hg <- hg %>%
			dplyr::select(trait, num_trees, mtry, min_node_size) %>%
			dplyr::filter(trait %in% TRAITS)
		counts <- table(factor(hg$trait, levels = TRAITS))
		if (any(counts != 1L)) {
			stop(
				"Hyperparameter file must contain exactly one row per trait: ", path,
				call. = FALSE
			)
		}
		hg
	}) %>% set_names(LEAF_TYPES)

	message(sprintf(
		"Preflight: CV dataset (%d rows) and hyperparameters validated.",
		nrow(data_ci_source)
	))
}


# ══════════════════════════════════════════════════════════════════════════════
# 1. Load models and test splits
# ══════════════════════════════════════════════════════════════════════════════

models <- map(LEAF_TYPES, function(lt) {
	map(TRAITS, function(tr) {
		path <- file.path(PATH_MODELS, sprintf("rf_%s_%s.rds", tr, lt))
		mod <- read_rds(path)
		if (!inherits(mod, "ranger")) {
			stop("Expected a ranger model in: ", path, call. = FALSE)
		}
		mod
	}) %>% set_names(TRAITS)
}) %>% set_names(LEAF_TYPES)

splits <- map(LEAF_TYPES, function(lt) {
	path <- file.path(PATH_MODELS, sprintf("splits_%s.rds", lt))
	x <- read_rds(path)
	if (!is.list(x) || !all(c("train", "test") %in% names(x))) {
		stop("Split file must contain named train and test data frames: ", path,
				 call. = FALSE)
	}
	assert_columns(x$test, c("PID_rep", COVARIATES), paste0("Test split ", path))
	if (nrow(x$test) == 0L) stop("Test split is empty: ", path, call. = FALSE)
	x
}) %>% set_names(LEAF_TYPES)

message("Models loaded: ", paste(LEAF_TYPES, collapse = " / "))
message(sprintf("fastshap version: %s", packageVersion("fastshap")))


# ══════════════════════════════════════════════════════════════════════════════
# 2. Compute SHAP values on held-out test set
# ══════════════════════════════════════════════════════════════════════════════
# One explain() call per trait × leaf type.
# X is the test set feature matrix (covariates only).
# SHAP values are computed out-of-fit: the model never saw these observations
# during training, matching the approach used for performance evaluation in
# 03_rf_fit.R and following Maynard et al. (2022).

message(sprintf("\nComputing SHAP values (nsim = %d per model)...", N_SIM))

shap_values <- map_dfr(LEAF_TYPES, function(lt) {
	
	test_df <- splits[[lt]]$test
	
	X_test <- test_df %>%
		dplyr::select(all_of(COVARIATES)) %>%
		as.matrix()
	if (!is.null(MAIN_EXPLAIN_MAX) && nrow(X_test) > MAIN_EXPLAIN_MAX) {
		set.seed(BASE_SEED + match(lt, LEAF_TYPES))
		keep <- sample.int(nrow(X_test), MAIN_EXPLAIN_MAX, replace = FALSE)
		test_df <- test_df[keep, , drop = FALSE]
		X_test <- X_test[keep, , drop = FALSE]
	}
	
	message(sprintf("  %s: n_test = %d", lt, nrow(test_df)))
	
	map_dfr(TRAITS, function(tr) {
		
		message(sprintf("    · %s", tr))
		
		mod <- models[[lt]][[tr]]
		set.seed(
			BASE_SEED + 100L * match(lt, LEAF_TYPES) + match(tr, TRAITS)
		)
		
		shap_mat <- fastshap::explain(
			object       = mod,
			X            = X_test,
			pred_wrapper = predict_fn,
			nsim         = N_SIM,
			parallel     = PARALLEL,
			adjust       = TRUE
		)
		
		# SHAP values to long format
		shap_long <- as.data.frame(shap_mat) %>%
			rowid_to_column(".row") %>%
			pivot_longer(-.row, names_to = "variable", values_to = "shap_value")
		
		# Feature values to long format (for dependence plots)
		feat_long <- as.data.frame(X_test) %>%
			rowid_to_column(".row") %>%
			mutate(PID_rep = test_df$PID_rep) %>%
			pivot_longer(-c(.row, PID_rep),
									 names_to  = "variable",
									 values_to = "feature_value")
		
		left_join(shap_long, feat_long, by = c(".row", "variable")) %>%
			mutate(trait = tr, leaf_type = lt) %>%
			dplyr::select(trait, leaf_type, PID_rep, variable, shap_value, feature_value)
	})
})

write_rds(shap_values, output_rds("shap_values"))
message(sprintf("\nSHAP complete: %d rows saved to %s",
								nrow(shap_values), output_rds("shap_values")))


# ══════════════════════════════════════════════════════════════════════════════
# 3. Variable importance: sum of absolute SHAP per predictor
# ══════════════════════════════════════════════════════════════════════════════
# Importance = sum(|SHAP|) across test plots, following Maynard et al. (2022).

shap_per_var <- shap_values %>%
	group_by(trait, leaf_type, variable) %>%
	summarise(
		sum_abs_shap  = sum(abs(shap_value),  na.rm = TRUE),
		mean_abs_shap = mean(abs(shap_value), na.rm = TRUE),
		mean_shap     = mean(shap_value,      na.rm = TRUE),
		.groups = "drop"
	) %>%
	mutate(
		variable_label = recode(variable, !!!VAR_LABELS),
		trait_label    = recode(trait,    !!!TRAIT_LABELS)
	)

write_rds(shap_per_var, output_rds("shap_per_var"))


# ══════════════════════════════════════════════════════════════════════════════
# 4. RQ1: Environmental vs successional filtering ratio
# ══════════════════════════════════════════════════════════════════════════════
# Sum absolute SHAP within two ecological categories, then compute their ratio.
#
# Ratio = total_environmental_shap / total_successional_shap
#
# Ratio > 1: environment explains more variation than stand age.
# Ratio = 1: equal importance.
# Ratio < 1: succession dominates (expected for height, shade tolerance).
#
# Using a ratio rather than proportions avoids the zero-sum framing of
# percentages and produces a single interpretable number per trait:
# "environmental filtering explains X times more variation than succession."
#
# The category-level summation (5 env predictors summed before dividing by
# 1 successional predictor) means the ratio reflects ecological signal, not
# predictor count.

shap_importance <- shap_per_var %>%
	mutate(category = case_when(
		variable %in% ENV_VARS  ~ "environmental",
		variable %in% SUCC_VARS ~ "successional",
		TRUE ~ NA_character_
	)) %>%
	filter(!is.na(category)) %>%
	group_by(trait, leaf_type, category) %>%
	summarise(category_shap = sum(sum_abs_shap), .groups = "drop") %>%
	pivot_wider(names_from = category, values_from = category_shap) %>%
	mutate(
		env_succ_ratio = environmental / successional,
		total_shap     = environmental + successional,
		prop_env       = environmental / total_shap,
		prop_succ      = successional  / total_shap,
		trait_label    = recode(trait, !!!TRAIT_LABELS)
	) %>%
	arrange(leaf_type, desc(env_succ_ratio))

write_rds(shap_importance, output_rds("shap_importance"))

message("\n── RQ1: Environmental vs successional filtering ratio ──────────────")
shap_importance %>%
	dplyr::select(trait_label, leaf_type, env_succ_ratio, prop_env) %>%
	mutate(across(where(is.numeric), ~ round(., 2))) %>%
	print(n = Inf)


# ══════════════════════════════════════════════════════════════════════════════
# 4b. Uncertainty on the env/succ ratio via repeated cross-validation (optional)
# ══════════════════════════════════════════════════════════════════════════════
# The ratio above is a single point estimate from one held-out test set. To show
# how robust it is to sampling and model-fit variability, we repeat k-fold CV a
# number of times, refit the RF on each training partition, recompute SHAP on the
# held-out fold, and form the env/succ ratio from that fold's summed |SHAP|.
# Folds are assigned at the PLOT level (PID), so all repeated inventories of a
# plot stay together and never leak between train and test. The 2.5/97.5
# percentiles across the repeat × fold ratios give a 95% interval per trait ×
# leaf type. The point estimates reported in the text are unchanged; this only
# adds error bars to Fig. 2a.
#
# Cost: refits / re-explains 18 models × N_REPEATS_CI × N_FOLDS_CI times. nsim is
# lower than the main run (aggregation across folds compensates) and the explained
# set can be capped (EXPLAIN_MAX). On a 32-core node the defaults run in a few
# hours; lower N_REPEATS_CI / EXPLAIN_MAX to trade precision for speed, or set
# RUN_SHAP_CI <- FALSE to skip entirely.

N_FOLDS_CI   <- if (SMOKE_TEST) 2L    else 10L
N_REPEATS_CI <- if (SMOKE_TEST) 1L    else 5L
N_SIM_CI     <- if (SMOKE_TEST) 3L    else 50L
EXPLAIN_MAX  <- if (SMOKE_TEST) 100L  else 1500L
SMOKE_MAX_PLOTS_PER_LEAF <- 1000L
MIN_TRAIN_CI <- if (SMOKE_TEST) 20L else 200L
MIN_TEST_CI  <- if (SMOKE_TEST) 10L else 50L
CI_LEVEL     <- 0.95

requested_cores <- suppressWarnings(as.integer(Sys.getenv("SHAP_CV_CORES", "32")))
if (length(requested_cores) != 1L || is.na(requested_cores) || requested_cores < 1L) {
	stop("SHAP_CV_CORES must be a positive integer.", call. = FALSE)
}
MAX_CV_CORES <- if (SMOKE_TEST) min(2L, requested_cores) else requested_cores

explain_tag <- if (is.null(EXPLAIN_MAX)) "all" else as.character(EXPLAIN_MAX)
CV_RUN_TAG <- sprintf(
	"%s_f%02d_r%02d_s%03d_n%s_seed%d",
	if (SMOKE_TEST) "smoke" else "normal",
	N_FOLDS_CI, N_REPEATS_CI, N_SIM_CI, explain_tag, BASE_SEED
)
CHECKPOINT_ROOT <- file.path(PATH_TABLES, "shap_cv_checkpoints")
CHECKPOINT_DIR  <- file.path(CHECKPOINT_ROOT, CV_RUN_TAG)

job_id <- function(trait, leaf_type, rep) {
	sprintf("%s__%s__rep%03d", leaf_type, trait, as.integer(rep))
}

job_checkpoint_path <- function(trait, leaf_type, rep) {
	file.path(CHECKPOINT_DIR, paste0(job_id(trait, leaf_type, rep), ".rds"))
}

error_checkpoint_path <- function(trait, leaf_type, rep) {
	stamp <- format(Sys.time(), "%Y%m%dT%H%M%S")
	file.path(
		CHECKPOINT_DIR,
		sprintf(
			"%s__ERROR__%s__pid%d.rds",
			job_id(trait, leaf_type, rep), stamp, Sys.getpid()
		)
	)
}

# Write to a temporary file in the destination directory, then rename. A valid
# checkpoint is never overwritten; invalid checkpoints are quarantined before
# a job is resubmitted.
atomic_write_rds <- function(x, path) {
	dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
	tmp <- tempfile(
		pattern = paste0(".", basename(path), "."),
		tmpdir = dirname(path), fileext = ".tmp"
	)
	on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
	readr::write_rds(x, tmp)
	if (!file.rename(tmp, path)) {
		stop("Could not atomically move temporary file into place: ", path)
	}
	invisible(path)
}

read_valid_checkpoint <- function(path, expected_job = NULL) {
	if (!file.exists(path)) return(NULL)
	payload <- tryCatch(readRDS(path), error = function(e) NULL)
	if (!is.list(payload) ||
			!identical(payload$schema_version, 1L) ||
			!identical(payload$config_tag, CV_RUN_TAG) ||
			!is.data.frame(payload$result) ||
			nrow(payload$result) != N_FOLDS_CI) return(NULL)

	required_result_columns <- c(
		"trait", "leaf_type", "rep", "fold", "n_train", "n_test",
		"n_explained", "env_sum", "succ_sum", "env_succ_ratio", "prop_env"
	)
	if (!all(required_result_columns %in% names(payload$result))) return(NULL)
	if (!setequal(payload$result$fold, seq_len(N_FOLDS_CI))) return(NULL)
	if (any(!is.finite(payload$result$env_succ_ratio)) ||
			any(!is.finite(payload$result$prop_env))) return(NULL)

	if (!is.null(expected_job)) {
		if (!all(payload$result$trait == expected_job$trait) ||
				!all(payload$result$leaf_type == expected_job$leaf_type) ||
				!all(payload$result$rep == expected_job$rep)) return(NULL)
	}
	payload
}

run_cv_job <- function(job) {
	tr <- as.character(job$trait[[1L]])
	lt <- as.character(job$leaf_type[[1L]])
	rp <- as.integer(job$rep[[1L]])
	id <- job_id(tr, lt, rp)
	out_path <- job_checkpoint_path(tr, lt, rp)

	message(sprintf("[%s] starting", id))
	tryCatch({
		df_lt <- dplyr::filter(data_ci, leaf_type == lt)
		hg <- hyper_ci[[lt]]
		if (nrow(df_lt) == 0L) stop("No rows available for leaf type '", lt, "'.")
		if (is.null(hg)) stop("No hyperparameter table available for leaf type '", lt, "'.")

		# Stable plot-level folds: repeated inventories of a PID cannot leak across
		# train and test. Seeds vary by trait, leaf type, repeat, and fold.
		fold_seed <- BASE_SEED +
			100000L * match(lt, LEAF_TYPES) +
			1000L * match(tr, TRAITS) + rp
		set.seed(fold_seed)
		pids <- sort(unique(df_lt$PID))
		if (length(pids) < N_FOLDS_CI) {
			stop("Only ", length(pids), " unique PIDs for ", N_FOLDS_CI, " folds.")
		}
		fold_of_pid <- stats::setNames(
			sample(rep(seq_len(N_FOLDS_CI), length.out = length(pids))),
			as.character(pids)
		)
		df_lt$fold_id <- unname(fold_of_pid[as.character(df_lt$PID)])

		fold_results <- lapply(seq_len(N_FOLDS_CI), function(k) {
			train_df <- df_lt %>%
				dplyr::filter(fold_id != k) %>%
				dplyr::select(dplyr::all_of(c(tr, COVARIATES))) %>%
				tidyr::drop_na()
			test_df <- df_lt %>%
				dplyr::filter(fold_id == k) %>%
				dplyr::select(dplyr::all_of(c(tr, COVARIATES))) %>%
				tidyr::drop_na()

			if (nrow(train_df) < MIN_TRAIN_CI || nrow(test_df) < MIN_TEST_CI) {
				stop(sprintf(
					"fold %d has insufficient complete cases (train=%d, test=%d; minima=%d/%d)",
					k, nrow(train_df), nrow(test_df), MIN_TRAIN_CI, MIN_TEST_CI
				))
			}

			mod <- fit_rf_model(
				tr, train_df, COVARIATES, hg, num_threads = 1L
			)$trait_mod

			X_test <- dplyr::select(test_df, dplyr::all_of(COVARIATES))
			shap_seed <- fold_seed + 100L * k
			set.seed(shap_seed)
			if (!is.null(EXPLAIN_MAX) && nrow(X_test) > EXPLAIN_MAX) {
				keep <- sample.int(nrow(X_test), EXPLAIN_MAX, replace = FALSE)
				X_test <- X_test[keep, , drop = FALSE]
			}

			sm <- fastshap::explain(
				mod,
				X = as.matrix(X_test),
				pred_wrapper = predict_fn,
				nsim = N_SIM_CI,
				parallel = FALSE,
				adjust = TRUE
			)
			abs_sum <- colSums(abs(as.data.frame(sm)), na.rm = TRUE)
			missing_shap_columns <- setdiff(COVARIATES, names(abs_sum))
			if (length(missing_shap_columns) > 0L) {
				stop(
					"fold ", k, " SHAP output is missing: ",
					paste(missing_shap_columns, collapse = ", ")
				)
			}
			env_sum  <- sum(abs_sum[ENV_VARS], na.rm = TRUE)
			succ_sum <- sum(abs_sum[SUCC_VARS], na.rm = TRUE)
			if (!is.finite(env_sum) || !is.finite(succ_sum) || succ_sum <= 0) {
				stop("fold ", k, " produced a non-finite or zero SHAP denominator.")
			}

			tibble::tibble(
				trait = tr,
				leaf_type = lt,
				rep = rp,
				fold = k,
				n_train = nrow(train_df),
				n_test = nrow(test_df),
				n_explained = nrow(X_test),
				env_sum = env_sum,
				succ_sum = succ_sum,
				env_succ_ratio = env_sum / succ_sum,
				prop_env = env_sum / (env_sum + succ_sum)
			)
		})

		result <- dplyr::bind_rows(fold_results)
		payload <- list(
			schema_version = 1L,
			config_tag = CV_RUN_TAG,
			job_id = id,
			completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
			result = result
		)
		atomic_write_rds(payload, out_path)
		message(sprintf("[%s] complete; checkpoint saved", id))
		tibble::tibble(job_id = id, status = "complete", error = NA_character_)
	}, error = function(e) {
		err_path <- error_checkpoint_path(tr, lt, rp)
		error_payload <- list(
			schema_version = 1L,
			config_tag = CV_RUN_TAG,
			job_id = id,
			failed_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
			message = conditionMessage(e),
			call = paste(deparse(conditionCall(e)), collapse = " ")
		)
		try(atomic_write_rds(error_payload, err_path), silent = TRUE)
		message(sprintf("[%s] FAILED: %s", id, conditionMessage(e)))
		tibble::tibble(job_id = id, status = "failed", error = conditionMessage(e))
	})
}

run_shap_ci <- function() {
	dir.create(CHECKPOINT_DIR, recursive = TRUE, showWarnings = FALSE)

	# Full modelling dataset with leaf type derived exactly as in scripts 03/05/06.
	data_ci <- data_ci_source %>%
		mutate(leaf_type = case_when(
			biome_boreal_forests_or_taiga == 1 |
				biome_temperate_conifer_forests == 1 ~ "coniferous",
			biome_temperate_broadleaf_forests == 1 |
				biome_mediterranean_woodlands == 1 ~ "broadleaf",
			TRUE ~ NA_character_
		)) %>%
		filter(!is.na(leaf_type))

	if (SMOKE_TEST) {
		set.seed(BASE_SEED + 900000L)
		keep_pid <- purrr::map_dfr(LEAF_TYPES, function(lt) {
			pids <- sort(unique(data_ci$PID[data_ci$leaf_type == lt]))
			n_keep <- min(length(pids), SMOKE_MAX_PLOTS_PER_LEAF)
			tibble::tibble(
				leaf_type = lt,
				PID = sample(pids, n_keep, replace = FALSE)
			)
		})
		data_ci <- dplyr::semi_join(data_ci, keep_pid, by = c("leaf_type", "PID"))
		message(sprintf(
			"4b smoke-test data cap: %d rows across at most %d PIDs per leaf type.",
			nrow(data_ci), SMOKE_MAX_PLOTS_PER_LEAF
		))
	}

	leaf_counts <- table(factor(data_ci$leaf_type, levels = LEAF_TYPES))
	if (any(leaf_counts == 0L)) {
		stop(
			"No modelling rows after leaf-type derivation for: ",
			paste(names(leaf_counts)[leaf_counts == 0L], collapse = ", "),
			call. = FALSE
		)
	}

	# Already validated during preflight, before the ordinary SHAP work began.
	hyper_ci <- hyper_ci_source

	# Validate complete-case availability before any workers or expensive fits.
	for (lt in LEAF_TYPES) {
		for (tr in TRAITS) {
			x <- data_ci %>%
				dplyr::filter(leaf_type == lt) %>%
				dplyr::select(PID, dplyr::all_of(c(tr, COVARIATES))) %>%
				tidyr::drop_na()
			if (nrow(x) < MIN_TRAIN_CI + MIN_TEST_CI ||
					dplyr::n_distinct(x$PID) < N_FOLDS_CI) {
				stop(sprintf(
					"Insufficient complete data for %s/%s: %d rows, %d PIDs.",
					lt, tr, nrow(x), dplyr::n_distinct(x$PID)
				), call. = FALSE)
			}
		}
	}

	jobs <- expand.grid(
		trait = TRAITS,
		leaf_type = LEAF_TYPES,
		rep = seq_len(N_REPEATS_CI),
		stringsAsFactors = FALSE
	) %>%
		dplyr::arrange(leaf_type, trait, rep)

	# Detect valid checkpoints. Corrupt/incompatible files are preserved with an
	# .invalid timestamp and the corresponding job is safely recomputed.
	is_complete <- logical(nrow(jobs))
	for (i in seq_len(nrow(jobs))) {
		job <- jobs[i, , drop = FALSE]
		path <- job_checkpoint_path(job$trait, job$leaf_type, job$rep)
		is_complete[i] <- !is.null(read_valid_checkpoint(path, job))
		if (file.exists(path) && !is_complete[i]) {
			quarantine <- paste0(
				path, ".invalid_", format(Sys.time(), "%Y%m%dT%H%M%S")
			)
			if (!file.rename(path, quarantine)) {
				stop("Invalid checkpoint could not be quarantined: ", path, call. = FALSE)
			}
			message("Quarantined invalid checkpoint: ", basename(quarantine))
		}
	}

	pending_jobs <- jobs[!is_complete, , drop = FALSE]
	message(sprintf(
		"\n4b. SHAP CV uncertainty [%s]: %d/%d jobs already complete; %d pending.",
		CV_RUN_TAG, sum(is_complete), nrow(jobs), nrow(pending_jobs)
	))
	message("Checkpoint directory: ", CHECKPOINT_DIR)

	status <- tibble::tibble(
		job_id = character(), status = character(), error = character()
	)
	if (nrow(pending_jobs) > 0L) {
		detected_cores <- suppressWarnings(parallel::detectCores(logical = FALSE))
		if (length(detected_cores) != 1L || is.na(detected_cores)) detected_cores <- 1L
		n_cores <- max(1L, min(MAX_CV_CORES, detected_cores, nrow(pending_jobs)))

		cl <- NULL
		on.exit({
			if (!is.null(cl)) try(parallel::stopCluster(cl), silent = TRUE)
			foreach::registerDoSEQ()
		}, add = TRUE)

		cl <- parallel::makePSOCKcluster(n_cores, outfile = "")
		doParallel::registerDoParallel(cl)
		message(sprintf("Dispatching %d pending jobs on %d workers.",
								nrow(pending_jobs), n_cores))

		worker_exports <- c(
			"run_cv_job", "atomic_write_rds", "job_id", "job_checkpoint_path",
			"error_checkpoint_path", "fit_rf_model", "predict_fn",
			"data_ci", "hyper_ci", "COVARIATES", "ENV_VARS", "SUCC_VARS",
			"TRAITS", "LEAF_TYPES", "TRAIT_LABELS", "VAR_LABELS",
			"BASE_SEED", "N_FOLDS_CI", "N_SIM_CI", "EXPLAIN_MAX",
			"MIN_TRAIN_CI", "MIN_TEST_CI", "CV_RUN_TAG", "CHECKPOINT_DIR"
		)
		parallel::clusterExport(cl, worker_exports, envir = environment())
		parallel::clusterEvalQ(cl, {
			library(ranger)
			library(fastshap)
			library(dplyr)
			library(tidyr)
			library(tibble)
			library(readr)
			NULL
		})

		job_list <- split(pending_jobs, seq_len(nrow(pending_jobs)))
		status <- foreach::foreach(
			job = job_list,
			.combine = dplyr::bind_rows,
			.noexport = worker_exports,
			.errorhandling = "pass",
			.options.snow = list(preschedule = FALSE)
		) %dopar% {
			run_cv_job(job)
		}

		parallel::stopCluster(cl)
		cl <- NULL
		foreach::registerDoSEQ()
	}

	# Re-read every job from disk; worker return values alone never establish
	# completion. This also makes a rerun and a fresh run follow the same path.
	payloads <- vector("list", nrow(jobs))
	valid <- logical(nrow(jobs))
	for (i in seq_len(nrow(jobs))) {
		job <- jobs[i, , drop = FALSE]
		path <- job_checkpoint_path(job$trait, job$leaf_type, job$rep)
		payloads[[i]] <- read_valid_checkpoint(path, job)
		valid[i] <- !is.null(payloads[[i]])
	}
	ratio_folds <- dplyr::bind_rows(lapply(payloads[valid], `[[`, "result"))

	if (!all(valid)) {
		partial_path <- output_rds("shap_importance_cv_raw_partial")
		if (file.exists(partial_path)) {
			partial_backup <- paste0(
				partial_path, ".previous_",
				format(Sys.time(), "%Y%m%dT%H%M%S"), "_pid", Sys.getpid()
			)
			if (!file.rename(partial_path, partial_backup)) {
				stop("Could not preserve the previous partial CV output: ",
						 partial_path, call. = FALSE)
			}
		}
		atomic_write_rds(ratio_folds, partial_path)
		failed_jobs <- jobs[!valid, , drop = FALSE] %>%
			mutate(job_id = job_id(trait, leaf_type, rep))
		failure_report <- list(
			config_tag = CV_RUN_TAG,
			failed_jobs = failed_jobs,
			worker_status = status,
			partial_raw_path = partial_path,
			checkpoint_dir = CHECKPOINT_DIR
		)
		failure_path <- file.path(CHECKPOINT_DIR, "latest_failure_summary.rds")
		if (file.exists(failure_path)) {
			failure_path <- file.path(
				CHECKPOINT_DIR,
				paste0("failure_summary_", format(Sys.time(), "%Y%m%dT%H%M%S"), ".rds")
			)
		}
		atomic_write_rds(failure_report, failure_path)
		stop(
			sprintf(
				paste0(
					"SHAP CV stopped with %d incomplete job(s). Completed checkpoints are safe. ",
					"Partial raw results: %s. Failure details: %s. Rerun the same command to retry only incomplete jobs."
				),
				sum(!valid), partial_path, failure_path
			),
			call. = FALSE
		)
	}

	expected_rows <- nrow(jobs) * N_FOLDS_CI
	if (nrow(ratio_folds) != expected_rows) {
		stop(sprintf(
			"Internal completeness check failed: expected %d CV rows, found %d.",
			expected_rows, nrow(ratio_folds)
		), call. = FALSE)
	}

	# Save the combined raw CV results before calculating any interval summary.
	raw_path <- output_rds("shap_importance_cv_raw")
	if (file.exists(raw_path)) {
		backup_path <- paste0(
			raw_path, ".previous_", format(Sys.time(), "%Y%m%dT%H%M%S")
		)
		if (!file.rename(raw_path, backup_path)) {
			stop("Could not preserve the previous raw CV output: ", raw_path,
					 call. = FALSE)
		}
	}
	atomic_write_rds(ratio_folds, raw_path)
	message(sprintf("Combined raw CV results saved first: %s (%d rows).",
							raw_path, nrow(ratio_folds)))

	alpha <- (1 - CI_LEVEL) / 2
	shap_importance_ci <- ratio_folds %>%
		group_by(trait, leaf_type) %>%
		summarise(
			n_cv = n(),
			n_repeats = n_distinct(rep),
			n_folds = n_distinct(fold),
			ratio_cv_med = median(env_succ_ratio),
			ratio_lwr = quantile(env_succ_ratio, alpha, names = FALSE, type = 8),
			ratio_upr = quantile(env_succ_ratio, 1 - alpha, names = FALSE, type = 8),
			prop_env_med = median(prop_env),
			prop_env_lwr = quantile(prop_env, alpha, names = FALSE, type = 8),
			prop_env_upr = quantile(prop_env, 1 - alpha, names = FALSE, type = 8),
			.groups = "drop"
		) %>%
		left_join(
			shap_importance %>%
				dplyr::select(
					trait, leaf_type,
					ratio_point = env_succ_ratio,
					prop_env_point = prop_env
				),
			by = c("trait", "leaf_type")
		) %>%
		mutate(
			trait_label = recode(trait, !!!TRAIT_LABELS),
			ci_level = CI_LEVEL,
			cv_config = CV_RUN_TAG
		) %>%
		dplyr::relocate(
			trait, trait_label, leaf_type, ratio_point, ratio_lwr, ratio_upr,
			prop_env_point, prop_env_lwr, prop_env_upr
		)

	ci_path <- output_rds("shap_importance_ci")
	if (file.exists(ci_path)) {
		backup_path <- paste0(
			ci_path, ".previous_", format(Sys.time(), "%Y%m%dT%H%M%S")
		)
		if (!file.rename(ci_path, backup_path)) {
			stop("Could not preserve the previous CI output: ", ci_path, call. = FALSE)
		}
	}
	atomic_write_rds(shap_importance_ci, ci_path)

	message("\n── RQ1 env/succ point estimates + repeated-CV 95% intervals ───────")
	shap_importance_ci %>%
		dplyr::select(
			trait_label, leaf_type, ratio_point,
			ratio_cv_med, ratio_lwr, ratio_upr, n_cv
		) %>%
		mutate(across(where(is.numeric), ~ round(., 2))) %>%
		arrange(leaf_type, desc(ratio_point)) %>%
		print(n = Inf)

	message(sprintf(
		"4b complete: all %d jobs and %d fold estimates validated.",
		nrow(jobs), nrow(ratio_folds)
	))
}

if (RUN_SHAP_CI) {
	run_shap_ci()
} else {
	message("\n4b. Skipping SHAP CV uncertainty (--skip-ci / SHAP_SKIP_CI=true).")
}


# ══════════════════════════════════════════════════════════════════════════════
# 5. Console diagnostics
# ══════════════════════════════════════════════════════════════════════════════

message("\n── Top 2 predictors per trait (mean |SHAP|) ────────────────────────")
shap_per_var %>%
	group_by(trait_label, leaf_type) %>%
	slice_max(mean_abs_shap, n = 2) %>%
	dplyr::select(trait_label, leaf_type, variable_label, mean_abs_shap) %>%
	mutate(mean_abs_shap = round(mean_abs_shap, 4)) %>%
	arrange(leaf_type, trait_label) %>%
	print(n = Inf)

message("\n── Stand age SHAP summary by leaf type ─────────────────────────────")
shap_values %>%
	filter(variable == "standage") %>%
	group_by(leaf_type) %>%
	summarise(
		mean_abs = mean(abs(shap_value), na.rm = TRUE),
		sd_shap  = sd(shap_value,        na.rm = TRUE),
		.groups  = "drop"
	) %>%
	mutate(across(where(is.numeric), ~ round(., 5))) %>%
	print()

message("\n── 04_shap.R complete ──────────────────────────────────────────────")
message(sprintf("Tables: %s/", PATH_TABLES))
message(sprintf("Plots:  %s/", PATH_PLOTS))
