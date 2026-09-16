################################################################################
## Deprecated compatibility entry point.
## The single maintained SHAP implementation is 04_shap_CV_uncertainty.R.
################################################################################

canonical_script <- "scripts/04_shap_CV_uncertainty.R"
if (!file.exists(canonical_script)) {
	stop("Canonical SHAP script not found: ", canonical_script, call. = FALSE)
}
message(
	"04_shap_with_CV_uncertainty.R is retained for compatibility; delegating to ",
	canonical_script
)
source(canonical_script, chdir = FALSE)
