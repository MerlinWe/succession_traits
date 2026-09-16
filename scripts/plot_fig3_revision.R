################################################################################
## Figure 3 revision prototype
##
## Run from the Trait Succession project root:
##   source("scripts/plot_fig3_revision.R")
##
## This plotting-only script reads existing analysis tables. It does not refit
## models or overwrite the current figures/main/fig3_pdp.png.
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
	"tables/pdp_raw.rds",
	"tables/pdp_stats.rds",
	"tables/pdp_summary.rds",
	"tables/shap_per_var.rds"
)
missing_files <- c(theme_file, input_files)[
	!file.exists(c(theme_file, input_files))
]
if (length(missing_files) > 0L) {
	stop(
		"Required file(s) not found:\n  ", paste(missing_files, collapse = "\n  "),
		"\nRun this script from the project root.",
		call. = FALSE
	)
}

source(theme_file)

# Match the shortened labels used by the current main_plots.R.
names(COLS_ENV) <- stringr::str_remove(names(COLS_ENV), "\\s*PC$")
ENV_LABELS[] <- stringr::str_remove(ENV_LABELS, "\\s*PC$")

strip_pc <- function(df) {
	if ("variable_label" %in% names(df)) {
		df <- dplyr::mutate(
			df,
			variable_label = stringr::str_remove(
				as.character(variable_label), "\\s*PC$"
			)
		)
	}
	df
}

assert_columns <- function(x, required, description) {
	missing <- setdiff(required, names(x))
	if (length(missing) > 0L) {
		stop(
			description, " is missing column(s): ", paste(missing, collapse = ", "),
			call. = FALSE
		)
	}
}

pdp_raw <- readr::read_rds(input_files[[1L]]) %>% strip_pc()
pdp_stats <- readr::read_rds(input_files[[2L]]) %>% strip_pc()
pdp_summary <- readr::read_rds(input_files[[3L]]) %>% strip_pc()
shap_per_var <- readr::read_rds(input_files[[4L]]) %>% strip_pc()

assert_columns(
	pdp_raw,
	c("iteration", "leaf_type", "trait", "variable", "group", "standage", "yhat"),
	"tables/pdp_raw.rds"
)
assert_columns(
	pdp_stats,
	c("iteration", "leaf_type", "trait", "variable", "slope_diff"),
	"tables/pdp_stats.rds"
)
assert_columns(
	pdp_summary,
	c("leaf_type", "trait", "variable", "trait_label", "variable_label",
		"slope_median", "slope_lwr", "slope_upr"),
	"tables/pdp_summary.rds"
)
assert_columns(
	shap_per_var,
	c("trait", "leaf_type", "variable", "sum_abs_shap", "trait_label"),
	"tables/shap_per_var.rds"
)

# Preserve the trait order used in Figure 2 for cross-figure consistency.
trait_totals <- shap_per_var %>%
	group_by(trait, leaf_type) %>%
	summarise(total_shap = sum(sum_abs_shap), .groups = "drop")

trait_order_fig2 <- shap_per_var %>%
	left_join(trait_totals, by = c("trait", "leaf_type")) %>%
	mutate(pct = 100 * sum_abs_shap / total_shap) %>%
	filter(variable == "standage") %>%
	group_by(trait_label) %>%
	summarise(mean_succ_pct = mean(pct), .groups = "drop") %>%
	arrange(desc(mean_succ_pct)) %>%
	pull(trait_label)

# ==============================================================================
# Panel a: environmental trait differentiation in early vs later succession
# ==============================================================================

LATE_AGE <- 100L

early_late <- pdp_raw %>%
	group_by(iteration, leaf_type, trait, variable, group) %>%
	summarise(
		age_early = standage[which.min(abs(standage - 0))],
		yhat_early = yhat[which.min(abs(standage - 0))],
		yhat_late = yhat[which.min(abs(standage - LATE_AGE))],
		.groups = "drop"
	) %>%
	pivot_wider(
		names_from = group,
		values_from = c(age_early, yhat_early, yhat_late)
	) %>%
	mutate(
		delta_early = yhat_early_high - yhat_early_low,
		delta_late = yhat_late_high - yhat_late_low,
		abs_delta_early = abs(delta_early),
		abs_delta_late = abs(delta_late),
		# This paired quantity directly tests movement relative to the 1:1 line.
		gap_change = abs_delta_late - abs_delta_early
	)

early_late_summary <- early_late %>%
	group_by(leaf_type, trait, variable) %>%
	summarise(
		early_age_med = median(c(age_early_high, age_early_low), na.rm = TRUE),
		abs_early_med = median(abs_delta_early, na.rm = TRUE),
		abs_early_lwr = quantile(abs_delta_early, 0.025, na.rm = TRUE),
		abs_early_upr = quantile(abs_delta_early, 0.975, na.rm = TRUE),
		abs_late_med = median(abs_delta_late, na.rm = TRUE),
		abs_late_lwr = quantile(abs_delta_late, 0.025, na.rm = TRUE),
		abs_late_upr = quantile(abs_delta_late, 0.975, na.rm = TRUE),
		gap_change_med = median(gap_change, na.rm = TRUE),
		gap_change_lwr = quantile(gap_change, 0.025, na.rm = TRUE),
		gap_change_upr = quantile(gap_change, 0.975, na.rm = TRUE),
		gap_change_robust = gap_change_lwr > 0 | gap_change_upr < 0,
		gap_change_class = case_when(
			gap_change_lwr > 0 ~ "Widening",
			gap_change_upr < 0 ~ "Narrowing",
			TRUE ~ "Uncertain"
		),
		.groups = "drop"
	) %>%
	mutate(
		trait_label = recode(trait, !!!TRAIT_LABELS),
		variable_label = recode(variable, !!!ENV_LABELS),
		leaf_type = tools::toTitleCase(leaf_type),
		trait_label = factor(trait_label, levels = TRAIT_LABELS),
		variable_label = factor(variable_label, levels = names(COLS_ENV))
	)

# Report the actual early reference age so an age-0 label is not used if the
# partial-dependence grid starts later.
early_age_range <- range(early_late_summary$early_age_med, na.rm = TRUE)
if (all(abs(early_age_range) < 1e-8)) {
	early_axis_label <- expression(
		"Environmental trait differentiation at stand age 0 (" *
			"|" * Delta * "trait|" * ")"
	)
} else {
	early_axis_label <- expression(
		"Environmental trait differentiation in early succession (" *
			"|" * Delta * "trait|" * ")"
	)
	warning(sprintf(
		"The earliest modelled ages range from %.1f to %.1f years; panel a uses 'early succession' rather than 'age 0'.",
		early_age_range[[1L]], early_age_range[[2L]]
	), call. = FALSE)
}

ax_max_pa <- early_late_summary %>%
	summarise(across(
		c(abs_early_lwr, abs_early_upr, abs_late_lwr, abs_late_upr),
		~ max(., na.rm = TRUE)
	)) %>%
	unlist() %>%
	max() %>%
	{ ceiling(. * 10) / 10 }

if (!is.finite(ax_max_pa) || ax_max_pa <= 0) {
	stop("Could not determine a positive axis range for panel a.", call. = FALSE)
}

SHAPES_ENV_VARS <- c(
	"Temperature" = 16,
	"Soil water retention" = 17,
	"Precipitation" = 15,
	"Elevation" = 18,
	"Soil pH" = 8
)

region_labels <- tidyr::expand_grid(
	leaf_type = unique(early_late_summary$leaf_type),
	region = c("widen", "narrow")
) %>%
	mutate(
		x = if_else(region == "widen", 0.06, 0.94) * ax_max_pa,
		y = if_else(region == "widen", 0.90, 0.08) * ax_max_pa,
		hjust = if_else(region == "widen", 0, 1),
		label = if_else(
			region == "widen",
			"Environmental differences\nwiden (divergence)",
			"Environmental differences\nnarrow (convergence)"
		)
	)

panel_a <- ggplot(
	early_late_summary,
	aes(x = abs_early_med, y = abs_late_med)
) +
	geom_abline(slope = 1, intercept = 0, colour = "grey35", linewidth = 0.55) +
	geom_text(
		data = region_labels,
		aes(x = x, y = y, label = label, hjust = hjust),
		inherit.aes = FALSE,
		colour = "grey35",
		fontface = "italic",
		size = 2.8
	) +
	# Marginal bootstrap intervals. Non-zero caps keep short vertical intervals
	# visible rather than allowing the point symbol to hide them completely.
	geom_errorbar(
		aes(ymin = abs_late_lwr, ymax = abs_late_upr),
		colour = "grey55", linewidth = 0.4, width = 0.018 * ax_max_pa
	) +
	geom_errorbar(
		aes(
			xmin = abs_early_lwr, xmax = abs_early_upr,
			ymin = abs_late_med, ymax = abs_late_med
		),
		colour = "grey55", linewidth = 0.4,
		width = 0.018 * ax_max_pa, orientation = "y"
	) +
	geom_point(
		aes(
			colour = trait_label,
			shape = variable_label,
			alpha = gap_change_robust
		),
		size = 2.6,
		stroke = 0.5
	) +
	scale_colour_viridis_d(
		option = "turbo",
		name = NULL,
		guide = guide_legend(
			nrow = 3, order = 1,
			override.aes = list(size = 2.5, shape = 16, alpha = 1)
		)
	) +
	scale_shape_manual(
		values = SHAPES_ENV_VARS,
		name = NULL,
		guide = guide_legend(
			ncol = 2, order = 2,
			override.aes = list(size = 2.5, colour = "grey25", alpha = 1)
		)
	) +
	scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.25), guide = "none") +
	coord_equal(xlim = c(0, ax_max_pa), ylim = c(0, ax_max_pa), clip = "on") +
	facet_wrap(~ leaf_type, ncol = 2) +
	labs(
		x = early_axis_label,
		y = expression(
			"Environmental trait differentiation after 100 years (" *
				"|" * Delta * "trait|" * ")"
		)
	) +
	theme_succession(base_size = 9) +
	theme(
		legend.position = "bottom",
		legend.box = "horizontal",
		legend.spacing.x = unit(6, "mm"),
		legend.key.size = unit(3.5, "mm")
	)

# ==============================================================================
# Panel b: continuous environmental modification of successional trait change
# ==============================================================================

# Continuous bootstrap directional consistency: 0 = an even split between
# positive and negative differences; 1 = every bootstrap has the same sign.
slope_consistency <- pdp_stats %>%
	group_by(leaf_type, trait, variable) %>%
	summarise(
		p_positive = mean(slope_diff > 0, na.rm = TRUE),
		directional_consistency = 2 * abs(p_positive - 0.5),
		.groups = "drop"
	)

panel_b_data <- pdp_summary %>%
	left_join(
		slope_consistency,
		by = c("leaf_type", "trait", "variable")
	) %>%
	mutate(
		leaf_type = tools::toTitleCase(leaf_type),
		trait_label = factor(trait_label, levels = rev(trait_order_fig2)),
		variable_label = factor(
			variable_label,
			levels = c(
				"Temperature", "Soil water retention", "Precipitation",
				"Elevation", "Soil pH"
			)
		),
		# Traits were z-scored in 01_traits_prep.R. Multiplying by 100 expresses
		# the fitted difference in successional change per 100 years.
		delta_slope_100 = 100 * slope_median,
		abs_delta_slope_100 = abs(delta_slope_100)
	) %>%
	filter(
		!is.na(variable_label),
		is.finite(delta_slope_100),
		is.finite(directional_consistency)
	)

slope_limit <- max(abs(panel_b_data$delta_slope_100), na.rm = TRUE)
if (!is.finite(slope_limit) || slope_limit <= 0) {
	stop("Could not determine a positive effect-size range for panel b.",
			 call. = FALSE)
}

panel_b <- ggplot(
	panel_b_data,
	aes(x = variable_label, y = trait_label)
) +
	geom_point(
		aes(
			fill = delta_slope_100,
			size = abs_delta_slope_100,
			alpha = directional_consistency
		),
		shape = 21,
		colour = "grey25",
		stroke = 0.35
	) +
	scale_fill_gradient2(
		low = COL_NEG,
		mid = "grey96",
		high = COL_POS,
		midpoint = 0,
		limits = c(-slope_limit, slope_limit),
		name = "Difference in trait change per 100 years\n(upper - lower environmental quantile)",
		guide = guide_colourbar(
			order = 1,
			title.position = "top",
			barwidth = unit(55, "mm"),
			barheight = unit(3.5, "mm")
		)
	) +
	scale_size_continuous(
		range = c(1.5, 5.5),
		trans = "sqrt",
		guide = "none"
	) +
	scale_alpha_continuous(
		range = c(0.25, 1),
		limits = c(0, 1),
		breaks = c(0, 0.5, 1),
		labels = c("Low", "Moderate", "High"),
		name = "Bootstrap directional\nconsistency",
		guide = guide_legend(
			order = 2,
			title.position = "top",
			override.aes = list(size = 3.5, fill = "grey35")
		)
	) +
	scale_x_discrete(guide = guide_axis(angle = 35)) +
	facet_wrap(~ leaf_type, ncol = 2) +
	labs(
		title = "Environmental modification of successional trait change",
		x = NULL,
		y = NULL
	) +
	theme_succession(base_size = 9) +
	theme(
		panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
		legend.position = "bottom",
		plot.title = element_text(hjust = 0.5)
	)

fig3_revision <- panel_a / panel_b +
	plot_layout(heights = c(1.4, 1)) +
	plot_annotation(tag_levels = "a")

output_dir <- "figures/main"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
png_path <- file.path(output_dir, "fig3_pdp_revision.png")
pdf_path <- file.path(output_dir, "fig3_pdp_revision.pdf")

ggsave(
	filename = png_path,
	plot = fig3_revision,
	width = 180,
	height = 230,
	units = "mm",
	dpi = 400,
	bg = "white"
)
ggsave(
	filename = pdf_path,
	plot = fig3_revision,
	width = 180,
	height = 230,
	units = "mm",
	device = grDevices::pdf,
	bg = "white"
)

message("Figure 3 revision created without overwriting the current figure:")
message("  ", png_path)
message("  ", pdf_path)
message("Panel a fading now reflects the paired bootstrap interval for widening/narrowing relative to the 1:1 line.")
message("Panel b: fill and size encode signed effect magnitude; opacity encodes bootstrap directional consistency.")
