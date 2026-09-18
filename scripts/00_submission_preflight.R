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

PDP_RUN_DIR <- "tables/pdp_grouped_pid/normal_pid_b100_q25_seed42"
VECV_RUN_DIR <- "tables/vecv_grouped_pid/normal_pid_f10_r30_seed42"
PDP_AUDIT_DIR <- "tables/diagnostics/submission_audit_grouped_pdp"
VECV_AUDIT_DIR <- "tables/diagnostics/submission_audit_grouped_vecv"

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
	"scripts/audit_grouped_pdp.R",
	"scripts/plot_fig3_grouped_pdp.R",
	"scripts/06_vecv.R",
	"scripts/audit_grouped_vecv.R",
	"scripts/plot_fig4_grouped_vecv.R",
	"scripts/revise_manuscript_submission.py",
	file.path(PDP_RUN_DIR, "run_metadata.rds"),
	file.path(PDP_RUN_DIR, "job_status.csv"),
	file.path(PDP_RUN_DIR, "pdp_raw.rds"),
	file.path(PDP_RUN_DIR, "pdp_stats.rds"),
	file.path(PDP_RUN_DIR, "pdp_summary.rds"),
	file.path(VECV_RUN_DIR, "run_metadata.rds"),
	file.path(VECV_RUN_DIR, "vecv_raw_unfiltered.rds"),
	file.path(VECV_RUN_DIR, "vecv_raw.rds"),
	file.path(VECV_RUN_DIR, "vecv_summary.rds"),
	file.path(VECV_RUN_DIR, "vecv_divergence.rds"),
	file.path(PDP_AUDIT_DIR, "audit_checks.csv"),
	file.path(VECV_AUDIT_DIR, "audit_checks.csv"),
	"figures/main/fig3_pdp_grouped_review.png",
	"figures/main/fig3_pdp_grouped_review.pdf",
	"figures/main/fig4_vecv_grouped_review.png",
	"figures/main/fig4_vecv_grouped_review.pdf",
	"Main_Manuscript_Trait_Succession_audited_revision.docx",
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
writeLines(repository_state, file.path(OUTPUT_DIR, "repository_state.local.txt"))
writeLines(
	utils::capture.output(sessionInfo()),
	file.path(OUTPUT_DIR, "session_info.local.txt")
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
	"audited and complete",
	"Figure 3 plotting", "scripts/plot_fig3_grouped_pdp.R",
	"figures/main/fig3_pdp_grouped_review.*", "saved results only",
	"audited and complete",
	"Figure 4 VEcv", "scripts/06_vecv.R",
	"tables/vecv_grouped_pid/<configuration>/", "PID-grouped repeated CV",
	"audited and complete",
	"Figure 4 plotting", "scripts/plot_fig4_grouped_vecv.R",
	"figures/main/fig4_vecv_grouped_review.*", "saved results only",
	"audited and complete"
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
	"trait", "leaf_type", "ratio_point", "ratio_cv_med", "ratio_lwr", "ratio_upr",
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
		"minimum_cv_median_ratio", "maximum_cv_median_ratio"
	),
	value = c(
		nrow(shap_ci), nrow(shap_raw), all(fig2_cells$n_cv == 50L),
		all(shap_ci$ratio_lwr > 1), min(shap_ci$ratio_lwr),
		min(shap_ci$ratio_cv_med), max(shap_ci$ratio_cv_med)
	)
)
if (nrow(shap_ci) != 18L || nrow(shap_raw) != 900L ||
		!all(fig2_cells$n_cv == 50L) || !all(shap_ci$ratio_lwr > 1)) {
	stop("Figure 2 failed the submission-freeze completeness checks.", call. = FALSE)
}
write_csv(figure2_freeze, file.path(OUTPUT_DIR, "figure2_freeze_summary.csv"))

smoke_specs <- tribble(
	~analysis, ~checkpoint_dir, ~metadata_file, ~expected_jobs,
	"Figure 3 grouped PDP",
	"tables/pdp_grouped_pid/smoke_pid_b001_q25_seed42/checkpoints",
	"tables/pdp_grouped_pid/smoke_pid_b001_q25_seed42/run_metadata.rds",
	2L,
	"Figure 4 grouped VEcv",
	"tables/vecv_grouped_pid/smoke_pid_f02_r01_seed42/checkpoints",
	"tables/vecv_grouped_pid/smoke_pid_f02_r01_seed42/run_metadata.rds",
	2L,
	"Figure 3 grouped PDP (parallel)",
	"tables/pdp_grouped_pid/smoke_parallel_pid_b001_q25_seed42/checkpoints",
	"tables/pdp_grouped_pid/smoke_parallel_pid_b001_q25_seed42/run_metadata.rds",
	2L,
	"Figure 4 grouped VEcv (parallel)",
	"tables/vecv_grouped_pid/smoke_parallel_pid_f02_r01_seed42/checkpoints",
	"tables/vecv_grouped_pid/smoke_parallel_pid_f02_r01_seed42/run_metadata.rds",
	2L
)
smoke_status <- pmap_dfr(
	smoke_specs,
	function(analysis, checkpoint_dir, metadata_file, expected_jobs) {
	if (!dir.exists(checkpoint_dir) || !file.exists(metadata_file)) {
		return(tibble(
			analysis = analysis, exists = FALSE, all_jobs_complete = FALSE,
			expected_jobs = expected_jobs, n_jobs = NA_integer_,
			n_complete = NA_integer_, resampling_unit = NA_character_
		))
	}
	checkpoint_files <- list.files(
		checkpoint_dir, pattern = "\\.rds$", full.names = TRUE
	)
	payloads <- lapply(checkpoint_files, function(path) {
		tryCatch(read_rds(path), error = function(e) NULL)
	})
	complete <- vapply(
		payloads,
		function(x) is.list(x) && identical(x$status, "complete") &&
			is.data.frame(x$result) && nrow(x$result) > 0L,
		FUN.VALUE = logical(1)
	)
	metadata <- read_rds(metadata_file)
	unit <- if (!is.null(metadata$resampling_unit)) {
		metadata$resampling_unit
	} else {
		metadata$bootstrap_unit
	}
	tibble(
		analysis = analysis,
		exists = TRUE,
		all_jobs_complete = length(payloads) == expected_jobs && all(complete),
		expected_jobs = expected_jobs,
		n_jobs = length(payloads),
		n_complete = sum(complete),
		resampling_unit = unit
	)
})
if (!all(smoke_status$all_jobs_complete)) {
	failed <- smoke_status %>%
		filter(!all_jobs_complete) %>%
		transmute(detail = sprintf(
			"%s: %s complete checkpoints (expected %s)",
			analysis, n_complete, expected_jobs
		)) %>%
		pull(detail)
	stop(
		"At least one grouped-resampling smoke test is incomplete:\n  ",
		paste(failed, collapse = "\n  "),
		call. = FALSE
	)
}
write_csv(smoke_status, file.path(OUTPUT_DIR, "smoke_test_status.csv"))

# Production-run gate. Full payload validation and reconstruction live in the
# dedicated audit scripts; here we require the expected checkpoint inventory,
# matching production metadata, and a fresh core audit without analytical
# failures. Narrative rows may deliberately reject a legacy manuscript claim.
production_specs <- tribble(
	~analysis, ~run_dir, ~audit_file, ~result_file, ~config_id, ~expected_jobs,
	"Figure 3 grouped PDP", PDP_RUN_DIR,
	file.path(PDP_AUDIT_DIR, "audit_checks.csv"),
	file.path(PDP_RUN_DIR, "pdp_summary.rds"),
	"normal_pid_b100_q25_seed42", 1800L,
	"Figure 4 grouped VEcv", VECV_RUN_DIR,
	file.path(VECV_AUDIT_DIR, "audit_checks.csv"),
	file.path(VECV_RUN_DIR, "vecv_summary.rds"),
	"normal_pid_f10_r30_seed42", 540L
)
production_status <- pmap_dfr(
	production_specs,
	function(analysis, run_dir, audit_file, result_file, config_id,
			 expected_jobs) {
		metadata <- read_rds(file.path(run_dir, "run_metadata.rds"))
		checkpoint_files <- list.files(
			file.path(run_dir, "checkpoints"),
			pattern = "\\.rds$", full.names = TRUE
		)
		audit <- read_csv(audit_file, show_col_types = FALSE)
		core_audit <- audit %>% filter(section != "narrative")
		metadata_ok <- identical(metadata$config_id, config_id) &&
			identical(metadata$smoke_test, FALSE)
		checkpoint_count_ok <- length(checkpoint_files) == expected_jobs
		core_audit_ok <- nrow(core_audit) > 0L &&
			!any(core_audit$status == "FAIL")
		audit_fresh <- file.info(audit_file)$mtime >= file.info(result_file)$mtime
		tibble(
			analysis = analysis,
			config_id = config_id,
			expected_jobs = expected_jobs,
			n_checkpoints = length(checkpoint_files),
			metadata_ok = metadata_ok,
			checkpoint_count_ok = checkpoint_count_ok,
			core_audit_ok = core_audit_ok,
			audit_fresh = audit_fresh,
			all_checks_pass = metadata_ok && checkpoint_count_ok &&
				core_audit_ok && audit_fresh
		)
	}
)
if (!all(production_status$all_checks_pass)) {
	failed <- production_status %>%
		filter(!all_checks_pass) %>%
		transmute(detail = sprintf(
			"%s: checkpoints=%d/%d; metadata=%s; core audit=%s; fresh=%s",
			analysis, n_checkpoints, expected_jobs, metadata_ok, core_audit_ok,
			audit_fresh
		)) %>%
		pull(detail)
	stop(
		"At least one grouped production run failed preflight:\n  ",
		paste(failed, collapse = "\n  "),
		call. = FALSE
	)
}
write_csv(
	production_status,
	file.path(OUTPUT_DIR, "production_run_status.csv")
)

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
message("  Figure 3: 1,800/1,800 production checkpoints; core audit passed.")
message("  Figure 4: 540/540 production checkpoints; core audit passed.")
message("  Figure 3 and Figure 4 sequential and parallel smoke tests: complete.")
message("  Grouped-resampling and VEcv unit checks: complete.")
message("  Diagnostics: ", OUTPUT_DIR)
