################################################################################
## succession_traits: 05 — grouped partial-dependence bootstrap
##
## Tests whether environmental context modifies trait trajectories through
## succession. The bootstrap sampling unit is the permanent FIA plot (PID), so
## all repeat inventories from a plot are sampled together. Each
## trait × forest type × bootstrap iteration is independently checkpointed.
##
## Normal run:
##   Rscript scripts/05_pdp.R
## Smoke test:
##   Rscript scripts/05_pdp.R --smoke-test
## Parallel smoke test (also exercises the server worker path):
##   Rscript scripts/05_pdp.R --parallel-smoke
## Optional worker override:
##   PDP_N_CORES=16 Rscript scripts/05_pdp.R
##
## Output root:
##   tables/pdp_grouped_pid/<configuration>/
################################################################################

rm(list = ls())

required_packages <- c("ranger", "doParallel", "foreach", "tidyverse")
missing_packages <- required_packages[
	!vapply(required_packages, requireNamespace, quietly = TRUE,
				FUN.VALUE = logical(1))
]
if (length(missing_packages) > 0L) {
	stop("Missing required package(s): ", paste(missing_packages, collapse = ", "),
			 call. = FALSE)
}

library(ranger)
library(doParallel)
library(foreach)
library(tidyverse)

source("scripts/functions.R")
source("scripts/plot_theme.R")

# ── Configuration ─────────────────────────────────────────────────────────────

env_flag <- function(name) {
	tolower(Sys.getenv(name, unset = "false")) %in% c("1", "true", "yes", "y")
}

env_int <- function(name, default) {
	x <- suppressWarnings(as.integer(Sys.getenv(name, unset = "")))
	if (is.na(x) || x < 1L) as.integer(default) else x
}

cli_args <- commandArgs(trailingOnly = TRUE)
PARALLEL_SMOKE <- "--parallel-smoke" %in% cli_args
SMOKE_TEST <- "--smoke-test" %in% cli_args || PARALLEL_SMOKE ||
	env_flag("PDP_SMOKE_TEST")
ANALYSIS_VERSION <- "grouped-pid-v1"
BASE_SEED <- 42L
GROUP_COL <- "PID"

TRAITS_ALL <- c(
	"bark_thickness", "conduit_diam", "height", "leaf_density", "leaf_k",
	"root_depth", "seed_dry_mass", "shade_tolerance", "specific_leaf_area"
)
LEAF_TYPES <- c("broadleaf", "coniferous")
COVARIATES <- c("standage", "temp_pc", "soil_pc", "rain_pc", "elevation", "soil_ph")
ENV_VARS_ALL <- c("temp_pc", "soil_pc", "rain_pc", "elevation", "soil_ph")

TRAITS <- if (SMOKE_TEST) "height" else TRAITS_ALL
ENV_VARS <- if (SMOKE_TEST) "temp_pc" else ENV_VARS_ALL
N_BOOT <- if (SMOKE_TEST) 1L else 100L
PROBS <- c(0.25, 0.75)
EARLY_AGE <- 10L
LATE_AGE <- 100L
PDP_GRID <- seq(EARLY_AGE, LATE_AGE, by = 5L)
MIN_STRATUM_ROWS <- if (SMOKE_TEST) 15L else 50L
MIN_STRATUM_GROUPS <- if (SMOKE_TEST) 8L else 25L
SMOKE_MAX_GROUPS_PER_LEAF <- 1000L

default_cores <- max(1L, min(16L, parallel::detectCores(logical = FALSE) - 1L))
N_CORES <- if (SMOKE_TEST) {
	if (PARALLEL_SMOKE) 2L else 1L
} else {
	env_int("PDP_N_CORES", default_cores)
}
PARALLEL <- N_CORES > 1L

PATH_DATA <- "data_processed/fia_traits_clean.rds"
PATH_TABLES <- "tables"
HYPER_FILES <- file.path(PATH_TABLES, sprintf("perf_%s.csv", LEAF_TYPES))

CONFIG_ID <- sprintf(
	"%s_pid_b%03d_q%02d_seed%d",
	if (PARALLEL_SMOKE) "smoke_parallel" else if (SMOKE_TEST) "smoke" else "normal",
	N_BOOT, as.integer(100 * PROBS[1]), BASE_SEED
)
OUTPUT_DIR <- file.path(PATH_TABLES, "pdp_grouped_pid", CONFIG_ID)
CHECKPOINT_DIR <- file.path(OUTPUT_DIR, "checkpoints")
ERROR_DIR <- file.path(OUTPUT_DIR, "errors")
dir.create(CHECKPOINT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(ERROR_DIR, recursive = TRUE, showWarnings = FALSE)

message("\n── Grouped PDP bootstrap ───────────────────────────────────────────")
message("Mode: ", if (SMOKE_TEST) "SMOKE TEST" else "NORMAL")
message("Configuration: ", CONFIG_ID)
message("Bootstrap unit: ", GROUP_COL)
message("Comparison ages: ", EARLY_AGE, " and ", LATE_AGE, " years")
message("Output: ", OUTPUT_DIR)

# ── Validation and input loading ──────────────────────────────────────────────

required_files <- c(
	PATH_DATA, HYPER_FILES, "scripts/functions.R", "scripts/plot_theme.R"
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
	stop("Required file(s) missing:\n  ", paste(missing_files, collapse = "\n  "),
			 call. = FALSE)
}

data <- readr::read_rds(PATH_DATA)
required_columns <- unique(c(
	GROUP_COL, "PID_rep", TRAITS_ALL, COVARIATES,
	"biome_boreal_forests_or_taiga", "biome_temperate_conifer_forests",
	"biome_temperate_broadleaf_forests", "biome_mediterranean_woodlands"
))
missing_columns <- setdiff(required_columns, names(data))
if (length(missing_columns) > 0L) {
	stop("Input data lack required column(s): ",
			 paste(missing_columns, collapse = ", "), call. = FALSE)
}
if (anyNA(data[[GROUP_COL]])) {
	stop("Grouping column ", GROUP_COL, " contains missing values.", call. = FALSE)
}

data <- data %>%
	mutate(
		leaf_type = case_when(
			biome_boreal_forests_or_taiga == 1 |
				biome_temperate_conifer_forests == 1 ~ "coniferous",
			biome_temperate_broadleaf_forests == 1 |
				biome_mediterranean_woodlands == 1 ~ "broadleaf",
			TRUE ~ NA_character_
		)
	) %>%
	filter(leaf_type %in% LEAF_TYPES)

if (anyDuplicated(data$PID_rep)) {
	stop("PID_rep must uniquely identify inventory rows.", call. = FALSE)
}

if (SMOKE_TEST) {
	data <- map_dfr(seq_along(LEAF_TYPES), function(i) {
		lt <- LEAF_TYPES[i]
		d <- filter(data, leaf_type == lt)
		groups <- unique(d[[GROUP_COL]])
		set.seed(BASE_SEED + i)
		keep <- sample(groups, min(length(groups), SMOKE_MAX_GROUPS_PER_LEAF))
		filter(d, .data[[GROUP_COL]] %in% keep)
	})
	message("Smoke-test subset: at most ", SMOKE_MAX_GROUPS_PER_LEAF,
				" plots per forest type; one trait, one environmental axis and one bootstrap.")
}

group_counts <- data %>%
	group_by(leaf_type) %>%
	summarise(n_rows = n(), n_groups = n_distinct(.data[[GROUP_COL]]), .groups = "drop")
if (any(group_counts$n_groups < 2L)) {
	stop("Each forest type must contain at least two permanent plots.", call. = FALSE)
}
print(group_counts, n = Inf)

hyper_params <- map(HYPER_FILES, function(path) {
	x <- readr::read_csv(path, show_col_types = FALSE)
	required <- c("trait", "num_trees", "mtry", "min_node_size")
	missing <- setdiff(required, names(x))
	if (length(missing) > 0L) {
		stop(basename(path), " lacks: ", paste(missing, collapse = ", "),
				 call. = FALSE)
	}
	x %>% select(all_of(required)) %>% distinct()
}) %>% set_names(LEAF_TYPES)

for (lt in LEAF_TYPES) {
	if (!setequal(hyper_params[[lt]]$trait, TRAITS_ALL)) {
		stop("Hyperparameter table for ", lt,
				 " must contain exactly one row for every trait.", call. = FALSE)
	}
}
if (SMOKE_TEST) {
	hyper_params <- map(hyper_params, ~ mutate(.x, num_trees = pmin(num_trees, 50L)))
}

INPUT_SIGNATURE <- paste(
	unname(tools::md5sum(c(PATH_DATA, HYPER_FILES))), collapse = "__"
)
SETTINGS_SIGNATURE <- paste(
	ANALYSIS_VERSION, paste(TRAITS, collapse = ","),
	paste(ENV_VARS, collapse = ","), N_BOOT, paste(PROBS, collapse = ","),
	paste(PDP_GRID, collapse = ","), MIN_STRATUM_ROWS, MIN_STRATUM_GROUPS,
	SMOKE_MAX_GROUPS_PER_LEAF, BASE_SEED,
	PARALLEL_SMOKE,
	sep = "__"
)

# ── Checkpoint and seed helpers ───────────────────────────────────────────────

job_id <- function(trait, leaf_type, iteration) {
	sprintf("%s__%s__boot%03d", leaf_type, trait, as.integer(iteration))
}

checkpoint_path <- function(id) file.path(CHECKPOINT_DIR, paste0(id, ".rds"))
error_path <- function(id) file.path(ERROR_DIR, paste0(id, "__ERROR.rds"))

atomic_write_rds <- function(x, path) {
	dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
	tmp <- tempfile(pattern = paste0(basename(path), "__"), tmpdir = dirname(path))
	on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
	readr::write_rds(x, tmp)
	if (!file.rename(tmp, path)) stop("Atomic rename failed for ", path)
	invisible(path)
}

read_valid_checkpoint <- function(path, expected_id = NULL) {
	if (!file.exists(path)) return(NULL)
	x <- tryCatch(readr::read_rds(path), error = function(e) NULL)
	if (is.null(x) || !is.list(x) || !identical(x$status, "complete") ||
			!identical(x$analysis_version, ANALYSIS_VERSION) ||
			!identical(x$config_id, CONFIG_ID) ||
			!identical(x$input_signature, INPUT_SIGNATURE) ||
			!identical(x$settings_signature, SETTINGS_SIGNATURE) ||
			!is.data.frame(x$result) || nrow(x$result) == 0L) return(NULL)
	if (!is.null(expected_id) && !identical(x$job_id, expected_id)) return(NULL)
	x
}

bootstrap_seed <- function(leaf_type, iteration) {
	as.integer(BASE_SEED + 100000L * match(leaf_type, LEAF_TYPES) + iteration)
}

model_seed <- function(trait, leaf_type, variable, iteration, group) {
	as.integer(
		BASE_SEED + 10000000L * match(leaf_type, LEAF_TYPES) +
			1000000L * match(trait, TRAITS_ALL) +
			10000L * match(variable, ENV_VARS_ALL) +
			100L * as.integer(iteration) + match(group, c("low", "high"))
	)
}

environment_group_labels <- function(variable) {
	switch(
		variable,
		temp_pc = c(low = "Lower temperature", high = "Higher temperature"),
		soil_pc = c(low = "Lower water retention", high = "Higher water retention"),
		rain_pc = c(low = "Lower precipitation", high = "Higher precipitation"),
		elevation = c(low = "Lower elevation", high = "Higher elevation"),
		soil_ph = c(low = "Lower soil pH", high = "Higher soil pH"),
		c(low = "Lower quantile", high = "Upper quantile")
	)
}

# ── One independently resumable bootstrap job ─────────────────────────────────

run_pdp_job <- function(job) {
	tr <- as.character(job$trait)
	lt <- as.character(job$leaf_type)
	iter <- as.integer(job$iteration)
	id <- job_id(tr, lt, iter)
	out_path <- checkpoint_path(id)

	message("[", id, "] starting")
	tryCatch({
		df_job <- data %>%
			dplyr::filter(leaf_type == lt) %>%
			dplyr::select(dplyr::all_of(unique(c(GROUP_COL, tr, COVARIATES)))) %>%
			tidyr::drop_na()
		if (dplyr::n_distinct(df_job[[GROUP_COL]]) < 2L) {
			stop("Fewer than two complete plot groups remain for this trait.")
		}

		boot_data <- cluster_bootstrap(
			df_job, GROUP_COL, seed = bootstrap_seed(lt, iter)
		)

		result <- purrr::map_dfr(ENV_VARS, function(env) {
			strata <- stratify_within(boot_data, env, PROBS)
			names(strata) <- c("low", "high")

			support <- purrr::imap_dfr(strata, function(d, group_name) {
				tibble::tibble(
					group = group_name,
					n_rows = nrow(d),
					n_plot_groups = dplyr::n_distinct(d[[GROUP_COL]]),
					n_bootstrap_clusters = dplyr::n_distinct(d$.bootstrap_cluster),
					min_age = min(d$standage, na.rm = TRUE),
					max_age = max(d$standage, na.rm = TRUE)
				)
			})

			if (any(support$n_rows < MIN_STRATUM_ROWS) ||
					any(support$n_bootstrap_clusters < MIN_STRATUM_GROUPS)) {
				stop(sprintf(
					"%s strata too small: %s",
					env,
					paste0(
						support$group, "=", support$n_rows, " rows/",
						support$n_bootstrap_clusters, " sampled clusters",
						collapse = "; "
					)
				))
			}
			if (any(support$min_age > EARLY_AGE) || any(support$max_age < LATE_AGE)) {
				stop(sprintf(
					"%s lacks observed age support for the %d–%d year PDP grid.",
					env, EARLY_AGE, LATE_AGE
				))
			}

			labels <- environment_group_labels(env)
			purrr::imap_dfr(strata, function(d, group_name) {
				fit <- fit_rf_model(
					trait = tr,
					df_train = d,
					covariates = COVARIATES,
					hyper_parameters = hyper_params[[lt]],
					num_threads = 1L,
					seed = model_seed(tr, lt, env, iter, group_name)
				)$trait_mod

				compute_pdp(
					fit, d, "standage", grid_values = PDP_GRID
				) %>%
					mutate(
						group = group_name,
						group_label = unname(labels[group_name])
					) %>%
					left_join(support, by = "group")
			}) %>%
				mutate(variable = env)
		}) %>%
			mutate(
				iteration = iter,
				leaf_type = lt,
				trait = tr,
				trait_label = recode(tr, !!!TRAIT_LABELS),
				variable_label = recode(variable, !!!ENV_LABELS),
				job_id = id
			) %>%
			select(
				iteration, leaf_type, trait, trait_label, variable, variable_label,
				group, group_label, standage, yhat, n_rows, n_plot_groups,
				n_bootstrap_clusters, min_age, max_age, job_id
			)

		payload <- list(
			status = "complete",
			analysis_version = ANALYSIS_VERSION,
			config_id = CONFIG_ID,
			input_signature = INPUT_SIGNATURE,
			settings_signature = SETTINGS_SIGNATURE,
			job_id = id,
			job = list(trait = tr, leaf_type = lt, iteration = iter),
			bootstrap_seed = bootstrap_seed(lt, iter),
			completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
			result = result
		)
		atomic_write_rds(payload, out_path)
		message("[", id, "] complete; checkpoint saved")
		list(job_id = id, status = "complete", error = NA_character_)
	}, error = function(e) {
		error_payload <- list(
			status = "failed",
			analysis_version = ANALYSIS_VERSION,
			config_id = CONFIG_ID,
			input_signature = INPUT_SIGNATURE,
			settings_signature = SETTINGS_SIGNATURE,
			job_id = id,
			job = list(trait = tr, leaf_type = lt, iteration = iter),
			bootstrap_seed = bootstrap_seed(lt, iter),
			failed_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
			error = conditionMessage(e),
			calls = paste(utils::capture.output(sys.calls()), collapse = "\n")
		)
		try(atomic_write_rds(error_payload, error_path(id)), silent = TRUE)
		message("[", id, "] FAILED: ", conditionMessage(e))
		list(job_id = id, status = "failed", error = conditionMessage(e))
	})
}

# ── Execute or resume jobs ────────────────────────────────────────────────────

jobs <- tidyr::expand_grid(
	trait = TRAITS,
	leaf_type = LEAF_TYPES,
	iteration = seq_len(N_BOOT)
) %>%
	mutate(job_id = purrr::pmap_chr(list(trait, leaf_type, iteration), job_id))

complete_before <- vapply(
	jobs$job_id,
	function(id) !is.null(read_valid_checkpoint(checkpoint_path(id), id)),
	FUN.VALUE = logical(1)
)
remaining <- jobs[!complete_before, , drop = FALSE]
message(sprintf(
	"Jobs: %d total; %d already complete; %d remaining.",
	nrow(jobs), sum(complete_before), nrow(remaining)
))

run_remaining_jobs <- function() {
	if (nrow(remaining) == 0L) return(list())
	if (!PARALLEL) {
		foreach::registerDoSEQ()
		return(lapply(seq_len(nrow(remaining)), function(i) run_pdp_job(remaining[i, ])))
	}

	workers <- min(N_CORES, nrow(remaining))
	cl <- parallel::makeCluster(workers, outfile = "")
	on.exit({
		try(parallel::stopCluster(cl), silent = TRUE)
		foreach::registerDoSEQ()
	}, add = TRUE)
	doParallel::registerDoParallel(cl)
	message("Parallel backend: ", workers, " workers")

	export_vars <- c(
		"run_pdp_job", "job_id", "checkpoint_path", "error_path",
		"atomic_write_rds", "bootstrap_seed", "model_seed",
		"environment_group_labels", "cluster_bootstrap", "stratify_within",
		"fit_rf_model", "compute_pdp", "data", "hyper_params", "GROUP_COL",
		"COVARIATES", "ENV_VARS", "ENV_VARS_ALL", "TRAITS_ALL", "LEAF_TYPES",
		"PROBS", "PDP_GRID", "EARLY_AGE", "LATE_AGE", "MIN_STRATUM_ROWS",
		"MIN_STRATUM_GROUPS", "BASE_SEED", "TRAIT_LABELS", "ENV_LABELS",
		"ANALYSIS_VERSION", "CONFIG_ID", "INPUT_SIGNATURE",
		"SETTINGS_SIGNATURE",
		"CHECKPOINT_DIR", "ERROR_DIR"
	)
	parallel::clusterExport(cl, export_vars, envir = environment())
	parallel::clusterEvalQ(cl, {
		library(ranger)
		library(dplyr)
		library(tidyr)
		library(purrr)
		library(readr)
		library(stringr)
		library(tibble)
		NULL
	})

	job_list <- split(remaining, seq_len(nrow(remaining)))
	foreach::foreach(
		job = job_list,
		.packages = c(
			"ranger", "dplyr", "tidyr", "purrr", "readr", "stringr", "tibble"
		),
		.noexport = export_vars
	) %dopar% run_pdp_job(job)
}

job_status_new <- run_remaining_jobs()
job_status <- if (length(job_status_new) == 0L) {
	tibble(job_id = character(), status = character(), error = character())
} else {
	bind_rows(job_status_new)
}

payloads <- lapply(jobs$job_id, function(id) {
	read_valid_checkpoint(checkpoint_path(id), id)
})
complete_after <- !vapply(payloads, is.null, FUN.VALUE = logical(1))

status_table <- jobs %>%
	transmute(
		job_id, trait, leaf_type, iteration,
		status = if_else(complete_after, "complete", "incomplete")
	) %>%
	left_join(select(job_status, job_id, error), by = "job_id")
readr::write_csv(status_table, file.path(OUTPUT_DIR, "job_status.csv"))

completed_results <- payloads[complete_after] %>% map_dfr("result")
atomic_write_rds(completed_results, file.path(OUTPUT_DIR, "pdp_raw_partial.rds"))

if (!all(complete_after)) {
	failure_report <- status_table %>% filter(status != "complete")
	atomic_write_rds(failure_report, file.path(OUTPUT_DIR, "failure_report.rds"))
	stop(
		sprintf(
			"Grouped PDP bootstrap stopped with %d incomplete job(s). Completed checkpoints and the partial raw table are safe. Rerun the same command to resume.",
			sum(!complete_after)
		),
		call. = FALSE
	)
}

# Save completed raw bootstrap curves before any downstream summarisation.
pdp_raw <- completed_results
atomic_write_rds(pdp_raw, file.path(OUTPUT_DIR, "pdp_raw.rds"))

pdp_stats <- pdp_raw %>%
	group_by(iteration, leaf_type, trait, trait_label, variable, variable_label, group) %>%
	summarise(
		slope = unname(coef(lm(yhat ~ standage))[2]),
		intercept = unname(coef(lm(yhat ~ standage))[1]),
		yhat_early = yhat[standage == EARLY_AGE][1],
		yhat_late = yhat[standage == LATE_AGE][1],
		.groups = "drop"
	) %>%
	pivot_wider(
		names_from = group,
		values_from = c(slope, intercept, yhat_early, yhat_late)
	) %>%
	mutate(
		slope_diff = slope_high - slope_low,
		abs_rate_diff = abs(slope_high) - abs(slope_low),
		intercept_diff = intercept_high - intercept_low,
		gap_early = abs(yhat_early_high - yhat_early_low),
		gap_late = abs(yhat_late_high - yhat_late_low),
		gap_change = gap_late - gap_early
	)
atomic_write_rds(pdp_stats, file.path(OUTPUT_DIR, "pdp_stats.rds"))

pdp_summary <- pdp_stats %>%
	group_by(leaf_type, trait, trait_label, variable, variable_label) %>%
	summarise(
		slope_high_median = median(slope_high, na.rm = TRUE),
		slope_low_median = median(slope_low, na.rm = TRUE),
		slope_median = median(slope_diff, na.rm = TRUE),
		slope_lwr = quantile(slope_diff, 0.025, na.rm = TRUE),
		slope_upr = quantile(slope_diff, 0.975, na.rm = TRUE),
		abs_rate_diff_median = median(abs_rate_diff, na.rm = TRUE),
		intercept_median = median(intercept_diff, na.rm = TRUE),
		intercept_lwr = quantile(intercept_diff, 0.025, na.rm = TRUE),
		intercept_upr = quantile(intercept_diff, 0.975, na.rm = TRUE),
		gap_early_med = median(gap_early, na.rm = TRUE),
		gap_early_lwr = quantile(gap_early, 0.025, na.rm = TRUE),
		gap_early_upr = quantile(gap_early, 0.975, na.rm = TRUE),
		gap_late_med = median(gap_late, na.rm = TRUE),
		gap_late_lwr = quantile(gap_late, 0.025, na.rm = TRUE),
		gap_late_upr = quantile(gap_late, 0.975, na.rm = TRUE),
		gap_change_med = median(gap_change, na.rm = TRUE),
		gap_change_lwr = quantile(gap_change, 0.025, na.rm = TRUE),
		gap_change_upr = quantile(gap_change, 0.975, na.rm = TRUE),
		n_boot = n_distinct(iteration),
		.groups = "drop"
	) %>%
	mutate(
		slope_direction_stable = slope_lwr > 0 | slope_upr < 0,
		intercept_direction_stable = intercept_lwr > 0 | intercept_upr < 0,
		# Backward-compatible aliases used by the current plotting scripts.
		slope_robust = slope_direction_stable,
		intercept_robust = intercept_direction_stable,
		gap_change_class = case_when(
			gap_change_lwr > 0 ~ "widening",
			gap_change_upr < 0 ~ "narrowing",
			TRUE ~ "uncertain"
		),
		modification_magnitude_100yr = 100 * abs(slope_median),
		intervals_provisional = n_boot < 20L
	)
atomic_write_rds(pdp_summary, file.path(OUTPUT_DIR, "pdp_summary.rds"))

run_metadata <- list(
	analysis_version = ANALYSIS_VERSION,
	config_id = CONFIG_ID,
	smoke_test = SMOKE_TEST,
	bootstrap_unit = GROUP_COL,
	n_boot = N_BOOT,
	traits = TRAITS,
	environmental_variables = ENV_VARS,
	quantile_probabilities = PROBS,
	pdp_grid = PDP_GRID,
	early_age = EARLY_AGE,
	late_age = LATE_AGE,
	minimum_stratum_rows = MIN_STRATUM_ROWS,
	minimum_stratum_groups = MIN_STRATUM_GROUPS,
	base_seed = BASE_SEED,
	input_signature = INPUT_SIGNATURE,
	settings_signature = SETTINGS_SIGNATURE,
	script_md5 = unname(tools::md5sum("scripts/05_pdp.R")),
	functions_md5 = unname(tools::md5sum("scripts/functions.R")),
	plot_theme_md5 = unname(tools::md5sum("scripts/plot_theme.R")),
	completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
	session_info = utils::capture.output(sessionInfo())
)
atomic_write_rds(run_metadata, file.path(OUTPUT_DIR, "run_metadata.rds"))

message("\nGrouped PDP bootstrap complete.")
message("  Jobs completed: ", sum(complete_after), "/", nrow(jobs))
message("  Raw curve rows: ", nrow(pdp_raw))
message("  Summary combinations: ", nrow(pdp_summary))
message("  Output directory: ", OUTPUT_DIR)
