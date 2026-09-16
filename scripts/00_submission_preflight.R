################################################################################
## Submission-hardening preflight and provenance freeze
##
## Read-only with respect to analytical outputs. Diagnostics are written to:
##   tables/diagnostics/submission_freeze/
##
## Run from the project root:
##   Rscript scripts/00_submission_preflight.R
################################################################################

rm(list = ls())

required_packages <- c("dplyr", "purrr", "readr", "tibble", "tidyr")
missing_packages <- required_packages[
	!vapply(required_packages, requireNamespace, quietly = TRUE,
				FUN.VALUE = logical(1))
]
if (length(missing_packages) > 0L) {
	stop("Missing required package(s): ", paste(missing_packages, collapse = ", "),
			 call. = FALSE)
}

library(dplyr)
library(purrr)
library(readr)
library(tibble)
library(tidyr)

source("scripts/functions.R")

OUTPUT_DIR <- "tables/diagnostics/submission_freeze"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

required_files <- c(
	"data_processed/fia_traits_clean.rds",
	"tables/perf_broadleaf.csv",
	"tables/perf_coniferous.csv",
	"tables/shap_importance.rds",
	"tables/shap_per_var.rds",
	"tables/shap_importance_ci.rds",
	"tables/shap_importance_cv_raw.rds",
	"scripts/00_submission_preflight.R",
	"scripts/functions.R",
	"scripts/plot_theme.R",
	"scripts/04_shap.R",
	"scripts/04_shap_CV_uncertainty.R",
	"scripts/04_shap_with_CV_uncertainty.R",
	"scripts/plot_fig2_shap_cv_uncertainty.R",
	"scripts/05_pdp.R",
	"scripts/06_vecv.R",
	"scripts/SUBMISSION_HARDENING.md"
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
	stop("Required submission file(s) missing:\n  ",
			 paste(missing_files, collapse = "\n  "), call. = FALSE)
}

git_value <- function(args) {
	paste(system2("git", args, stdout = TRUE, stderr = TRUE), collapse = "\n")
}

repository_state <- c(
	paste0("generated_at_utc: ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
	paste0("git_commit: ", git_value(c("rev-parse", "HEAD"))),
	paste0("git_branch: ", git_value(c("branch", "--show-current"))),
	"git_status:",
	git_value(c("status", "--short"))
)
writeLines(repository_state, file.path(OUTPUT_DIR, "repository_state.txt"))
writeLines(
	utils::capture.output(sessionInfo()),
	file.path(OUTPUT_DIR, "session_info.txt")
)

pipeline_manifest <- tribble(
	~stage, ~canonical_script, ~primary_outputs, ~resampling_unit, ~status,
	"Data preparation", "scripts/01_traits_prep.R",
	"data_processed/fia_traits_clean.rds", "inventory processing", "existing",
	"Environmental PCA", "scripts/02_environment_pca.R",
	"data_processed/fia_traits_clean.rds", "inventory processing", "existing",
	"Random forest fits", "scripts/03_rf_fit.R",
	"models/; tables/perf_*.csv", "legacy row split", "provenance review pending",
	"Figure 2 SHAP", "scripts/04_shap_CV_uncertainty.R",
	"tables/shap_*.rds", "PID-grouped repeated CV for uncertainty",
	"audited and complete",
	"Figure 2 plotting", "scripts/plot_fig2_shap_cv_uncertainty.R",
	"figures/main/fig2_shap_with_cv_uncertainty.*", "saved results only",
	"audited and complete",
	"Figure 3 PDP", "scripts/05_pdp.R",
	"tables/pdp_grouped_pid/<configuration>/", "PID cluster bootstrap",
	"hardened; smoke tested",
	"Figure 4 VEcv", "scripts/06_vecv.R",
	"tables/vecv_grouped_pid/<configuration>/", "PID-grouped repeated CV",
	"hardened; smoke tested"
)
write_csv(pipeline_manifest, file.path(OUTPUT_DIR, "pipeline_manifest.csv"))

data <- read_rds("data_processed/fia_traits_clean.rds")
forest_data <- data %>%
	mutate(
		leaf_type = case_when(
			biome_boreal_forests_or_taiga == 1 |
				biome_temperate_conifer_forests == 1 ~ "coniferous",
			biome_temperate_broadleaf_forests == 1 |
				biome_mediterranean_woodlands == 1 ~ "broadleaf",
			TRUE ~ NA_character_
		)
	) %>%
	filter(!is.na(leaf_type))

design_summary <- forest_data %>%
	group_by(leaf_type) %>%
	summarise(
		n_inventory_rows = n(),
		n_permanent_plots = n_distinct(PID),
		n_repeated_plots = sum(table(PID) > 1L),
		max_inventories_per_plot = max(table(PID)),
		.groups = "drop"
	)
write_csv(design_summary, file.path(OUTPUT_DIR, "design_summary.csv"))

shap_ci <- read_rds("tables/shap_importance_ci.rds")
shap_raw <- read_rds("tables/shap_importance_cv_raw.rds")
required_ci <- c(
	"trait", "leaf_type", "ratio_point", "ratio_lwr", "ratio_upr",
	"n_cv", "n_repeats", "n_folds"
)
required_raw <- c("trait", "leaf_type", "rep", "fold", "env_succ_ratio")
if (!all(required_ci %in% names(shap_ci)) ||
		!all(required_raw %in% names(shap_raw))) {
	stop("Figure 2 uncertainty tables do not match the expected schema.",
			 call. = FALSE)
}

fig2_cells <- shap_raw %>%
	count(trait, leaf_type, name = "n_cv")
figure2_freeze <- tibble(
	check = c(
		"trait_by_forest_cells", "raw_fold_estimates", "all_cells_have_50_estimates",
		"all_lower_intervals_above_one", "minimum_lower_interval",
		"minimum_point_ratio", "maximum_point_ratio"
	),
	value = c(
		nrow(shap_ci), nrow(shap_raw), all(fig2_cells$n_cv == 50L),
		all(shap_ci$ratio_lwr > 1), min(shap_ci$ratio_lwr),
		min(shap_ci$ratio_point), max(shap_ci$ratio_point)
	)
)
if (nrow(shap_ci) != 18L || nrow(shap_raw) != 900L ||
		!all(fig2_cells$n_cv == 50L) || !all(shap_ci$ratio_lwr > 1)) {
	stop("Figure 2 failed the submission-freeze completeness checks.", call. = FALSE)
}
write_csv(figure2_freeze, file.path(OUTPUT_DIR, "figure2_freeze_summary.csv"))

smoke_specs <- tribble(
	~analysis, ~status_file, ~metadata_file,
	"Figure 3 grouped PDP",
	"tables/pdp_grouped_pid/smoke_pid_b001_q25_seed42/job_status.csv",
	"tables/pdp_grouped_pid/smoke_pid_b001_q25_seed42/run_metadata.rds",
	"Figure 4 grouped VEcv",
	"tables/vecv_grouped_pid/smoke_pid_f02_r01_seed42/job_status.csv",
	"tables/vecv_grouped_pid/smoke_pid_f02_r01_seed42/run_metadata.rds",
	"Figure 3 grouped PDP (parallel)",
	"tables/pdp_grouped_pid/smoke_parallel_pid_b001_q25_seed42/job_status.csv",
	"tables/pdp_grouped_pid/smoke_parallel_pid_b001_q25_seed42/run_metadata.rds",
	"Figure 4 grouped VEcv (parallel)",
	"tables/vecv_grouped_pid/smoke_parallel_pid_f02_r01_seed42/job_status.csv",
	"tables/vecv_grouped_pid/smoke_parallel_pid_f02_r01_seed42/run_metadata.rds"
)
smoke_status <- pmap_dfr(smoke_specs, function(analysis, status_file, metadata_file) {
	if (!file.exists(status_file) || !file.exists(metadata_file)) {
		return(tibble(
			analysis = analysis, exists = FALSE, all_jobs_complete = FALSE,
			n_jobs = NA_integer_, n_complete = NA_integer_, resampling_unit = NA_character_
		))
	}
	status <- read_csv(status_file, show_col_types = FALSE)
	metadata <- read_rds(metadata_file)
	unit <- if (!is.null(metadata$resampling_unit)) {
		metadata$resampling_unit
	} else {
		metadata$bootstrap_unit
	}
	tibble(
		analysis = analysis,
		exists = TRUE,
		all_jobs_complete = nrow(status) > 0L && all(status$status == "complete"),
		n_jobs = nrow(status),
		n_complete = sum(status$status == "complete"),
		resampling_unit = unit
	)
})
if (!all(smoke_status$all_jobs_complete)) {
	stop("At least one grouped-resampling smoke test is incomplete.", call. = FALSE)
}
write_csv(smoke_status, file.path(OUTPUT_DIR, "smoke_test_status.csv"))

# Deterministic unit checks for the two grouped-resampling primitives and VEcv.
toy <- tibble(
	PID = rep(sprintf("plot%02d", 1:12), times = rep(c(1L, 2L, 3L), 4L)),
	value = seq_len(24)
)
fold_a <- assign_group_folds(toy, "PID", v = 3L, seed = 42L)
fold_b <- assign_group_folds(toy, "PID", v = 3L, seed = 42L)
folds_per_pid <- tapply(fold_a, toy$PID, function(x) length(unique(x)))

boot_a <- cluster_bootstrap(toy, "PID", seed = 42L)
boot_b <- cluster_bootstrap(toy, "PID", seed = 42L)
cluster_to_pid <- boot_a %>%
	group_by(.bootstrap_cluster) %>%
	summarise(n_pid = n_distinct(PID), .groups = "drop")
original_sizes <- count(toy, PID, name = "expected_rows")
bootstrap_sizes <- boot_a %>%
	count(.bootstrap_cluster, PID, name = "observed_rows") %>%
	left_join(original_sizes, by = "PID")

resampling_tests <- tibble(
	check = c(
		"grouped_folds_are_deterministic",
		"each_PID_occurs_in_one_fold",
		"cluster_bootstrap_is_deterministic",
		"each_bootstrap_cluster_is_one_PID",
		"bootstrap_retains_all_rows_of_each_drawn_PID",
		"VEcv_perfect_prediction_equals_one",
		"VEcv_mean_prediction_equals_zero"
	),
	passed = c(
		identical(fold_a, fold_b),
		all(folds_per_pid == 1L),
		identical(boot_a, boot_b),
		all(cluster_to_pid$n_pid == 1L),
		all(bootstrap_sizes$observed_rows == bootstrap_sizes$expected_rows),
		isTRUE(all.equal(VEcv(1:5, 1:5), 1)),
		isTRUE(all.equal(VEcv(1:5, rep(mean(1:5), 5)), 0))
	)
)
if (!all(resampling_tests$passed)) {
	stop("At least one grouped-resampling unit check failed.", call. = FALSE)
}
write_csv(
	resampling_tests,
	file.path(OUTPUT_DIR, "grouped_resampling_unit_checks.csv")
)

manifest_paths <- unique(c(
	required_files,
	"figures/main/fig2_shap_with_cv_uncertainty.png",
	"figures/main/fig2_shap_with_cv_uncertainty.pdf",
	list.files(
		"tables/shap_cv_checkpoints/normal_f10_r05_s050_n1500_seed42",
		full.names = TRUE, recursive = TRUE
	),
	list.files(
		"tables/pdp_grouped_pid/smoke_pid_b001_q25_seed42",
		full.names = TRUE, recursive = TRUE
	),
	list.files(
		"tables/vecv_grouped_pid/smoke_pid_f02_r01_seed42",
		full.names = TRUE, recursive = TRUE
	),
	list.files(
		"tables/pdp_grouped_pid/smoke_parallel_pid_b001_q25_seed42",
		full.names = TRUE, recursive = TRUE
	),
	list.files(
		"tables/vecv_grouped_pid/smoke_parallel_pid_f02_r01_seed42",
		full.names = TRUE, recursive = TRUE
	)
))
manifest_paths <- manifest_paths[file.exists(manifest_paths) & !dir.exists(manifest_paths)]
info <- file.info(manifest_paths)
file_manifest <- tibble(
	path = manifest_paths,
	size_bytes = info$size,
	modified_at = format(info$mtime, tz = "UTC", usetz = TRUE),
	md5 = unname(tools::md5sum(manifest_paths))
) %>% arrange(path)
write_csv(file_manifest, file.path(OUTPUT_DIR, "file_manifest.csv"))

message("Submission preflight passed.")
message("  Figure 2: 18/18 cells complete; 900 fold estimates; all lower intervals > 1.")
message("  Figure 3 and Figure 4 sequential and parallel smoke tests: complete.")
message("  Grouped-resampling and VEcv unit checks: complete.")
message("  Diagnostics: ", OUTPUT_DIR)
