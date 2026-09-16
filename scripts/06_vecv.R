################################################################################
## succession_traits: 06 — grouped predictability analysis (VEcv)
##
## Repeated inventories from the same FIA plot are kept in the same fold.
## Each trait × forest type × repeat is independently checkpointed and can be
## resumed. Production outputs are versioned and never overwrite the legacy
## row-wise VEcv files in tables/.
##
## Normal run:
##   Rscript scripts/06_vecv.R
## Smoke test:
##   Rscript scripts/06_vecv.R --smoke-test
## Parallel smoke test (also exercises the server worker path):
##   Rscript scripts/06_vecv.R --parallel-smoke
## Optional worker override:
##   VECV_N_CORES=16 Rscript scripts/06_vecv.R
##
## Output root:
##   tables/vecv_grouped_pid/<configuration>/
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
	env_flag("VECV_SMOKE_TEST")
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
N_REPEATS <- if (SMOKE_TEST) 1L else 30L
N_FOLDS <- if (SMOKE_TEST) 2L else 10L
MIN_BIN_N <- if (SMOKE_TEST) 3L else 30L
PROBS <- c(0.25, 0.75)
STANDAGE_BREAKS <- seq(0, 150, by = 10)
SMOKE_MAX_GROUPS_PER_LEAF <- 250L

default_cores <- max(1L, min(16L, parallel::detectCores(logical = FALSE) - 1L))
N_CORES <- if (SMOKE_TEST) {
	if (PARALLEL_SMOKE) 2L else 1L
} else {
	env_int("VECV_N_CORES", default_cores)
}
PARALLEL <- N_CORES > 1L

PATH_DATA <- "data_processed/fia_traits_clean.rds"
PATH_TABLES <- "tables"
HYPER_FILES <- file.path(PATH_TABLES, sprintf("perf_%s.csv", LEAF_TYPES))

CONFIG_ID <- sprintf(
	"%s_pid_f%02d_r%02d_seed%d",
	if (PARALLEL_SMOKE) "smoke_parallel" else if (SMOKE_TEST) "smoke" else "normal",
	N_FOLDS, N_REPEATS, BASE_SEED
)
OUTPUT_DIR <- file.path(PATH_TABLES, "vecv_grouped_pid", CONFIG_ID)
CHECKPOINT_DIR <- file.path(OUTPUT_DIR, "checkpoints")
ERROR_DIR <- file.path(OUTPUT_DIR, "errors")
dir.create(CHECKPOINT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(ERROR_DIR, recursive = TRUE, showWarnings = FALSE)

message("\n── Grouped VEcv analysis ───────────────────────────────────────────")
message("Mode: ", if (SMOKE_TEST) "SMOKE TEST" else "NORMAL")
message("Configuration: ", CONFIG_ID)
message("Resampling unit: ", GROUP_COL)
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
				" plots per forest type; one trait and one environmental axis.")
}

group_counts <- data %>%
	group_by(leaf_type) %>%
	summarise(n_rows = n(), n_groups = n_distinct(.data[[GROUP_COL]]), .groups = "drop")
if (any(group_counts$n_groups < N_FOLDS)) {
	stop("At least one forest type contains fewer plot groups than CV folds.",
			 call. = FALSE)
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
	x %>%
		select(all_of(required)) %>%
		distinct()
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
	paste(ENV_VARS, collapse = ","), N_FOLDS, N_REPEATS,
	paste(PROBS, collapse = ","), paste(STANDAGE_BREAKS, collapse = ","),
	MIN_BIN_N, SMOKE_MAX_GROUPS_PER_LEAF, BASE_SEED,
	PARALLEL_SMOKE,
	sep = "__"
)

# ── Checkpoint helpers ────────────────────────────────────────────────────────

job_id <- function(trait, leaf_type, repeat_id) {
	sprintf("%s__%s__rep%03d", leaf_type, trait, as.integer(repeat_id))
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

job_seed <- function(trait, leaf_type, repeat_id) {
	as.integer(
		BASE_SEED + 1000000L * match(trait, TRAITS_ALL) +
			10000L * match(leaf_type, LEAF_TYPES) + 100L * as.integer(repeat_id)
	)
}

run_vecv_job <- function(job) {
	tr <- as.character(job$trait)
	lt <- as.character(job$leaf_type)
	rp <- as.integer(job$repeat_id)
	id <- job_id(tr, lt, rp)
	out_path <- checkpoint_path(id)

	message("[", id, "] starting")
	tryCatch({
		df_lt <- dplyr::filter(data, leaf_type == lt)
		result <- oof_skill_by_bins(
			trait = tr,
			data = df_lt,
			covariates = COVARIATES,
			env_vars = ENV_VARS,
			hyper_grid = hyper_params[[lt]],
			standage_breaks = STANDAGE_BREAKS,
			probs = PROBS,
			v = N_FOLDS,
			repeats = 1L,
			group_col = GROUP_COL,
			base_seed = job_seed(tr, lt, rp),
			repeat_ids = rp
		) %>%
			mutate(
				leaf_type = lt,
				trait_label = recode(trait, !!!TRAIT_LABELS),
				variable_label = recode(variable, !!!ENV_LABELS),
				standage_mid = as.numeric(stringr::str_extract(
					as.character(standage_bin), "(?<=\\[)\\d+"
				)) + 5,
				job_id = id
			)

		payload <- list(
			status = "complete",
			analysis_version = ANALYSIS_VERSION,
			config_id = CONFIG_ID,
			input_signature = INPUT_SIGNATURE,
			settings_signature = SETTINGS_SIGNATURE,
			job_id = id,
			job = list(trait = tr, leaf_type = lt, repeat_id = rp),
			seed = job_seed(tr, lt, rp),
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
			job = list(trait = tr, leaf_type = lt, repeat_id = rp),
			seed = job_seed(tr, lt, rp),
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
	repeat_id = seq_len(N_REPEATS)
) %>%
	mutate(job_id = purrr::pmap_chr(
		list(trait, leaf_type, repeat_id), job_id
	))

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
		return(lapply(seq_len(nrow(remaining)), function(i) run_vecv_job(remaining[i, ])))
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
		"run_vecv_job", "job_id", "checkpoint_path", "error_path",
		"atomic_write_rds", "job_seed", "oof_skill_by_bins",
		"assign_group_folds", "fit_rf_model", "VEcv", "E1",
		"data", "hyper_params", "COVARIATES", "ENV_VARS", "PROBS",
		"STANDAGE_BREAKS", "N_FOLDS", "GROUP_COL", "BASE_SEED",
		"TRAITS_ALL", "LEAF_TYPES", "TRAIT_LABELS", "ENV_LABELS",
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
		NULL
	})

	job_list <- split(remaining, seq_len(nrow(remaining)))
	foreach::foreach(
		job = job_list,
		.packages = c("ranger", "dplyr", "tidyr", "purrr", "readr", "stringr"),
		.noexport = export_vars
	) %dopar% run_vecv_job(job)
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
		job_id, trait, leaf_type, repeat_id,
		status = if_else(complete_after, "complete", "incomplete")
	) %>%
	left_join(select(job_status, job_id, error), by = "job_id")
readr::write_csv(status_table, file.path(OUTPUT_DIR, "job_status.csv"))

completed_results <- payloads[complete_after] %>% map_dfr("result")
atomic_write_rds(completed_results, file.path(OUTPUT_DIR, "vecv_raw_partial.rds"))

if (!all(complete_after)) {
	failure_report <- status_table %>% filter(status != "complete")
	atomic_write_rds(failure_report, file.path(OUTPUT_DIR, "failure_report.rds"))
	stop(
		sprintf(
			"Grouped VEcv stopped with %d incomplete job(s). Completed checkpoints and the partial raw table are safe. Rerun the same command to resume.",
			sum(!complete_after)
		),
		call. = FALSE
	)
}

# Save all completed raw results before filtering or summarising.
vecv_raw_unfiltered <- completed_results
atomic_write_rds(vecv_raw_unfiltered, file.path(OUTPUT_DIR, "vecv_raw_unfiltered.rds"))

vecv_raw <- vecv_raw_unfiltered %>% filter(n >= MIN_BIN_N)
if (nrow(vecv_raw) == 0L) {
	stop("No VEcv cells remain after the minimum-bin-size filter.", call. = FALSE)
}
atomic_write_rds(vecv_raw, file.path(OUTPUT_DIR, "vecv_raw.rds"))

vecv_summary <- vecv_raw %>%
	group_by(
		trait, trait_label, leaf_type, variable, variable_label,
		env_group, standage_bin, standage_mid
	) %>%
	summarise(
		n_med = median(n, na.rm = TRUE),
		n_groups_med = median(n_groups, na.rm = TRUE),
		VEcv_med = median(VEcv, na.rm = TRUE),
		VEcv_lwr = quantile(VEcv, 0.025, na.rm = TRUE),
		VEcv_upr = quantile(VEcv, 0.975, na.rm = TRUE),
		E1_med = median(E1, na.rm = TRUE),
		E1_lwr = quantile(E1, 0.025, na.rm = TRUE),
		E1_upr = quantile(E1, 0.975, na.rm = TRUE),
		n_repeats = n_distinct(repeat_id),
		.groups = "drop"
	)
atomic_write_rds(vecv_summary, file.path(OUTPUT_DIR, "vecv_summary.rds"))

vecv_divergence_raw <- vecv_raw %>%
	select(
		trait, trait_label, leaf_type, variable, variable_label,
		standage_bin, standage_mid, repeat_id, env_group, VEcv
	) %>%
	pivot_wider(names_from = env_group, values_from = VEcv, names_prefix = "VEcv_") %>%
	mutate(
		delta_VEcv = VEcv_high - VEcv_low,
		abs_delta_VEcv = abs(delta_VEcv)
	)
atomic_write_rds(
	vecv_divergence_raw,
	file.path(OUTPUT_DIR, "vecv_divergence_raw.rds")
)

vecv_divergence <- vecv_divergence_raw %>%
	group_by(
		trait, trait_label, leaf_type, variable, variable_label,
		standage_bin, standage_mid
	) %>%
	summarise(
		delta_med = median(delta_VEcv, na.rm = TRUE),
		delta_lwr = quantile(delta_VEcv, 0.025, na.rm = TRUE),
		delta_upr = quantile(delta_VEcv, 0.975, na.rm = TRUE),
		abs_delta_med = median(abs_delta_VEcv, na.rm = TRUE),
		abs_delta_lwr = quantile(abs_delta_VEcv, 0.025, na.rm = TRUE),
		abs_delta_upr = quantile(abs_delta_VEcv, 0.975, na.rm = TRUE),
		n_repeats = n_distinct(repeat_id),
		.groups = "drop"
	) %>%
	mutate(
		direction_stable = delta_lwr > 0 | delta_upr < 0,
		# Backward-compatible alias used by the current supplementary plot code.
		sig_divergence = direction_stable
	)
atomic_write_rds(vecv_divergence, file.path(OUTPUT_DIR, "vecv_divergence.rds"))

run_metadata <- list(
	analysis_version = ANALYSIS_VERSION,
	config_id = CONFIG_ID,
	smoke_test = SMOKE_TEST,
	resampling_unit = GROUP_COL,
	n_folds = N_FOLDS,
	n_repeats = N_REPEATS,
	traits = TRAITS,
	environmental_variables = ENV_VARS,
	minimum_bin_n = MIN_BIN_N,
	base_seed = BASE_SEED,
	input_signature = INPUT_SIGNATURE,
	settings_signature = SETTINGS_SIGNATURE,
	script_md5 = unname(tools::md5sum("scripts/06_vecv.R")),
	functions_md5 = unname(tools::md5sum("scripts/functions.R")),
	plot_theme_md5 = unname(tools::md5sum("scripts/plot_theme.R")),
	completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
	session_info = utils::capture.output(sessionInfo())
)
atomic_write_rds(run_metadata, file.path(OUTPUT_DIR, "run_metadata.rds"))

message("\nGrouped VEcv complete.")
message("  Jobs completed: ", sum(complete_after), "/", nrow(jobs))
message("  Raw rows before filtering: ", nrow(vecv_raw_unfiltered))
message("  Raw rows retained: ", nrow(vecv_raw))
message("  Output directory: ", OUTPUT_DIR)
