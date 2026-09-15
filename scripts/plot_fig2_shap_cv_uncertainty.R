################################################################################
## Rebuild Figure 2 with repeated-CV uncertainty in panel a
##
## Run from the project root after 04_shap_with_CV_uncertainty.R completes:
##   Rscript scripts/plot_fig2_shap_cv_uncertainty.R
##
## This script only reads completed tables and creates a figure. It does not
## refit models or recompute SHAP values, and it does not overwrite fig2_shap.
################################################################################

rm(list = ls())

required_packages <- c("tidyverse", "patchwork")
missing_packages <- required_packages[
	!vapply(required_packages, requireNamespace, quietly = TRUE,
					FUN.VALUE = logical(1))
]
if (length(missing_packages) > 0L) {
	stop(
		"Missing required package(s): ", paste(missing_packages, collapse = ", "),
		call. = FALSE
	)
}

library(tidyverse)
library(patchwork)

theme_file <- "scripts/plot_theme.R"
input_files <- c(
	"tables/shap_per_var.rds",
	"tables/shap_importance.rds",
	"tables/shap_importance_ci.rds"
)

missing_files <- c(theme_file, input_files)[
	!file.exists(c(theme_file, input_files))
]
if (length(missing_files) > 0L) {
	stop(
		"Required file(s) not found:\n  ", paste(missing_files, collapse = "\n  "),
		"\nRun this script from the project root and confirm that the normal CV run finished.",
		call. = FALSE
	)
}

source(theme_file)

# Match the shortened labels used in the current main_plots.R.
names(COLS_ENV) <- stringr::str_remove(names(COLS_ENV), "\\s*PC$")
ENV_LABELS[] <- stringr::str_remove(ENV_LABELS, "\\s*PC$")
COLS_PREDICTORS <- c("Stand age" = COL_STANDAGE, COLS_ENV)
ALL_PREDICTOR_LABELS <- c("standage" = "Stand age", ENV_LABELS)

assert_columns <- function(x, required, description) {
	missing <- setdiff(required, names(x))
	if (length(missing) > 0L) {
		stop(
			description, " is missing column(s): ", paste(missing, collapse = ", "),
			call. = FALSE
		)
	}
}

shap_per_var <- readr::read_rds(input_files[[1L]])
shap_importance <- readr::read_rds(input_files[[2L]])
shap_importance_ci <- readr::read_rds(input_files[[3L]])

assert_columns(
	shap_per_var,
	c("trait", "leaf_type", "variable", "sum_abs_shap", "trait_label"),
	"tables/shap_per_var.rds"
)
assert_columns(
	shap_importance,
	c("trait", "leaf_type", "env_succ_ratio", "trait_label"),
	"tables/shap_importance.rds"
)
assert_columns(
	shap_importance_ci,
	c("trait", "leaf_type", "ratio_lwr", "ratio_upr"),
	"tables/shap_importance_ci.rds"
)

duplicate_ci <- shap_importance_ci %>%
	count(trait, leaf_type) %>%
	filter(n != 1L)
if (nrow(duplicate_ci) > 0L) {
	stop("The CV table must contain exactly one row per trait and leaf type.",
			 call. = FALSE)
}

ratio_plot_data <- shap_importance %>%
	left_join(
		shap_importance_ci %>%
			dplyr::select(trait, leaf_type, ratio_lwr, ratio_upr),
		by = c("trait", "leaf_type")
	)

if (anyNA(ratio_plot_data$ratio_lwr) || anyNA(ratio_plot_data$ratio_upr)) {
	missing_ci <- ratio_plot_data %>%
		filter(is.na(ratio_lwr) | is.na(ratio_upr)) %>%
		transmute(id = paste(leaf_type, trait, sep = "/")) %>%
		pull(id)
	stop(
		"Missing CV intervals for: ", paste(missing_ci, collapse = ", "),
		call. = FALSE
	)
}
if (any(!is.finite(ratio_plot_data$env_succ_ratio)) ||
		any(!is.finite(ratio_plot_data$ratio_lwr)) ||
		any(!is.finite(ratio_plot_data$ratio_upr)) ||
		any(ratio_plot_data$env_succ_ratio <= 0) ||
		any(ratio_plot_data$ratio_lwr <= 0) ||
		any(ratio_plot_data$ratio_upr <= 0) ||
		any(ratio_plot_data$ratio_lwr > ratio_plot_data$ratio_upr)) {
	stop("Ratio estimates or interval bounds are invalid for a log2 axis.",
			 call. = FALSE)
}

# Panel b data: unchanged from the current Figure 2.
trait_totals <- shap_per_var %>%
	group_by(trait, leaf_type) %>%
	summarise(total_shap = sum(sum_abs_shap), .groups = "drop")

shap_pct <- shap_per_var %>%
	left_join(trait_totals, by = c("trait", "leaf_type")) %>%
	mutate(pct = 100 * sum_abs_shap / total_shap)

trait_order_fig2 <- shap_pct %>%
	filter(variable == "standage") %>%
	group_by(trait_label) %>%
	summarise(mean_succ_pct = mean(pct), .groups = "drop") %>%
	arrange(desc(mean_succ_pct)) %>%
	pull(trait_label)

ratio_trait_order <- ratio_plot_data %>%
	group_by(trait_label) %>%
	summarise(mean_ratio = mean(env_succ_ratio), .groups = "drop") %>%
	arrange(mean_ratio) %>%
	pull(trait_label)

# Panel a: original held-out point estimates plus empirical repeated-CV
# 2.5th-97.5th percentile intervals. The CV median is deliberately not used as
# the point because the manuscript's reported point estimates remain unchanged.
fig2a <- ratio_plot_data %>%
	mutate(
		trait_label = factor(trait_label, levels = ratio_trait_order),
		leaf_type = tools::toTitleCase(leaf_type)
	) %>%
	ggplot(aes(x = env_succ_ratio, y = trait_label, colour = leaf_type)) +
	geom_vline(
		xintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.4
	) +
	geom_segment(
		aes(x = 1, xend = env_succ_ratio, yend = trait_label),
		linewidth = 0.5, alpha = 0.25
	) +
	geom_segment(
		aes(x = ratio_lwr, xend = ratio_upr, yend = trait_label),
		linewidth = 1.1, alpha = 0.8, lineend = "round"
	) +
	geom_point(size = 3.5) +
	scale_colour_manual(values = COLS_LEAFTYPE, name = NULL) +
	scale_x_continuous(
		trans = "log2",
		breaks = c(0.25, 0.5, 1, 2, 5, 10, 20, 50, 100),
		labels = function(x) paste0(x, "\u00d7")
	) +
	facet_wrap(~ leaf_type, ncol = 1) +
	labs(x = LAB_RATIO, y = NULL) +
	theme_succession(base_size = 9) +
	theme(
		legend.position = "none",
		panel.grid.major.y = element_blank()
	)

# Panel b: unchanged relative importance bars.
fig2b <- shap_pct %>%
	mutate(
		leaf_type = tools::toTitleCase(leaf_type),
		trait_label = factor(trait_label, levels = trait_order_fig2),
		variable_label = recode(variable, !!!ALL_PREDICTOR_LABELS),
		variable_label = factor(variable_label, levels = names(COLS_PREDICTORS))
	) %>%
	ggplot(aes(x = trait_label, y = pct, fill = variable_label)) +
	geom_col(colour = "white", linewidth = 0.25, width = 0.8) +
	facet_wrap(~ leaf_type, ncol = 1) +
	scale_fill_manual(
		name = NULL,
		values = COLS_PREDICTORS,
		breaks = names(COLS_PREDICTORS),
		guide = guide_legend(nrow = 2)
	) +
	scale_y_continuous(
		expand = expansion(mult = c(0, 0.02)),
		labels = function(x) paste0(x, "%")
	) +
	scale_x_discrete(guide = guide_axis(angle = 35)) +
	labs(x = NULL, y = "Relative importance (% of total |SHAP|)") +
	theme_succession(base_size = 9) +
	guides(fill = guide_legend(nrow = 6)) +
	theme(
		legend.position = "right",
		panel.grid.major.x = element_blank()
	)

fig2 <- (fig2a | fig2b) +
	plot_layout(widths = c(0.4, 0.6)) +
	theme(plot.margin = margin(4, 4, 4, 4)) +
	plot_annotation(tag_levels = "a")

output_dir <- "figures/main"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
png_path <- file.path(output_dir, "fig2_shap_with_cv_uncertainty.png")
pdf_path <- file.path(output_dir, "fig2_shap_with_cv_uncertainty.pdf")

ggsave(
	filename = png_path, plot = fig2,
	width = 240, height = 160, units = "mm", dpi = 400, bg = "white"
)
ggsave(
	filename = pdf_path, plot = fig2,
	width = 240, height = 160, units = "mm", device = grDevices::pdf,
	bg = "white"
)

message("Figure 2 with SHAP CV uncertainty complete.")
message("  PNG: ", png_path)
message("  PDF: ", pdf_path)
message("Panel a points are the original held-out estimates; horizontal lines are repeated-CV 95% empirical intervals.")