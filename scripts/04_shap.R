################################################################################
## Compatibility entry point for the canonical, checkpointed SHAP analysis.
##
## Normal run: Rscript scripts/04_shap.R
## Smoke test: Rscript scripts/04_shap.R --smoke-test
################################################################################

canonical_script <- "scripts/04_shap_CV_uncertainty.R"
if (!file.exists(canonical_script)) {
	stop("Canonical SHAP script not found: ", canonical_script, call. = FALSE)
}
message("Delegating to canonical SHAP workflow: ", canonical_script)
source(canonical_script, chdir = FALSE)
