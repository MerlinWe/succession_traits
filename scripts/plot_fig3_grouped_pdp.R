################################################################################
## succession_traits: Figure 3 from the grouped PID-cluster PDP bootstrap
##
## Panel a asks whether environmental separation in modelled trait expression
## widens or narrows between stand ages 10 and 100 years.
## Panel b shows the magnitude of environmental modification of the modelled
## successional slope, without assigning a misleading "faster" label from the
## sign of the slope difference.
##
## Run from the project root:
##   Rscript scripts/audit_grouped_pdp.R
##   Rscript scripts/plot_fig3_grouped_pdp.R
################################################################################

rm(list = ls())

required_packages <- c("tidyverse", "patchwork")
missing_packages <- required_packages[
	!vapply(required_packages, requireNamespace, quietly = TRUE,
		FUN.VALUE = logical(1))
]
if (length(missing_packages) > 0L) {
	stop("Missing package(s): ", paste(missing_packages, collapse = ", "),
		call. = FALSE)
}
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(patchwork))
source("scripts/plot_theme.R")

RUN_DIR <- "tables/pdp_grouped_pid/normal_pid_b100_q25_seed42"
SUMMARY_PATH <- file.path(RUN_DIR, "pdp_summary.rds")
METADATA_PATH <- file.path(RUN_DIR, "run_metadata.rds")
OUT_DIR <- "figures/main"
DIAGNOSTIC_DIR <- "tables/diagnostics/submission_audit_grouped_pdp"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DIAGNOSTIC_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(SUMMARY_PATH) || !file.exists(METADATA_PATH)) {
	stop(
		"Grouped production output is missing. Run scripts/05_pdp.R first.",
		call. = FALSE
	)
}

metadata <- read_rds(METADATA_PATH)
if (!identical(metadata$config_id, "normal_pid_b100_q25_seed42") ||
	!identical(metadata$bootstrap_unit, "PID") ||
	!identical(metadata$n_boot, 100L) ||
	!identical(metadata$early_age, 10L) ||
	!identical(metadata$late_age, 100L)) {
	stop("Grouped PDP metadata do not match the production configuration.",
		call. = FALSE)
}

LEAF_LEVELS <- c("Broadleaf", "Coniferous")
ENV_LEVELS <- c(
	"Temperature", "Precipitation", "Soil water retention", "Elevation", "Soil pH"
)
ENV_SHAPES <- c(
	"Temperature" = 16,
	"Precipitation" = 15,
	"Soil water retention" = 17,
	"Elevation" = 18,
	"Soil pH" = 8
)
CHANGE_LEVELS <- c("Widening", "Uncertain", "Narrowing")
CHANGE_COLOURS <- c(
	"Widening" = "#D95F02",
	"Uncertain" = "grey70",
	"Narrowing" = "#4575B4"
)

pdp_summary <- read_rds(SUMMARY_PATH) %>%
	mutate(
		leaf_label = recode(
			leaf_type, broadleaf = "Broadleaf", coniferous = "Coniferous"
		),
		leaf_label = factor(leaf_label, levels = LEAF_LEVELS),
		variable_label = str_remove(variable_label, " PC$"),
		variable_label = factor(variable_label, levels = ENV_LEVELS),
		change_class = recode(
			gap_change_class,
			widening = "Widening",
			narrowing = "Narrowing",
			uncertain = "Uncertain"
		),
		change_class = factor(change_class, levels = CHANGE_LEVELS),
		modification_magnitude_100yr = 100 * abs(slope_median)
	) %>%
	filter(
		!is.na(leaf_label), !is.na(variable_label),
		is.finite(gap_early_med), is.finite(gap_late_med),
		is.finite(modification_magnitude_100yr)
	)

if (nrow(pdp_summary) != 90L || any(pdp_summary$n_boot != 100L)) {
	stop("Grouped PDP summary is incomplete.", call. = FALSE)
}

# -- Panel a: environmental separation at early and late succession -----------

panel_a_data <- pdp_summary
axis_max <- max(
	panel_a_data$gap_early_upr,
	panel_a_data$gap_late_upr,
	na.rm = TRUE
)
axis_max <- ceiling(axis_max * 10) / 10

region_labels <- tidyr::expand_grid(
	leaf_label = factor(LEAF_LEVELS, levels = LEAF_LEVELS),
	region = c("widen", "narrow")
) %>%
	mutate(
		x = if_else(region == "widen", 0.05, 0.95) * axis_max,
		y = if_else(region == "widen", 0.90, 0.08) * axis_max,
		hjust = if_else(region == "widen", 0, 1),
		label = if_else(
			region == "widen",
			"Environmental separation\nwidens",
			"Environmental separation\nnarrows"
		)
	)

panel_a <- ggplot(
	panel_a_data,
	aes(x = gap_early_med, y = gap_late_med)
) +
	geom_abline(slope = 1, intercept = 0, colour = "grey30", linewidth = 0.55) +
	geom_text(
		data = region_labels,
		aes(x = x, y = y, label = label, hjust = hjust),
		inherit.aes = FALSE,
		colour = "grey35", fontface = "italic", size = 2.8,
		lineheight = 0.95
	) +
	geom_errorbar(
		aes(ymin = gap_late_lwr, ymax = gap_late_upr),
		colour = "grey63", linewidth = 0.32, width = 0.018 * axis_max
	) +
	geom_errorbar(
		aes(
			xmin = gap_early_lwr, xmax = gap_early_upr,
			ymin = gap_late_med, ymax = gap_late_med
		),
		colour = "grey63", linewidth = 0.32,
		width = 0.018 * axis_max, orientation = "y"
	) +
	geom_point(
		aes(colour = change_class, shape = variable_label),
		size = 2.45, stroke = 0.45
	) +
	scale_colour_manual(
		values = CHANGE_COLOURS,
		name = "Change from age 10 to 100",
		drop = FALSE
	) +
	scale_shape_manual(
		values = ENV_SHAPES,
		name = "Environmental gradient",
		drop = FALSE
	) +
	coord_equal(
		xlim = c(0, axis_max), ylim = c(0, axis_max),
		expand = FALSE, clip = "on"
	) +
	facet_wrap(~ leaf_label, nrow = 1) +
	labs(
		title = "Environmental separation can widen or narrow through succession",
		subtitle = "Each point is one trait x environmental gradient; bars are 95% PID-cluster bootstrap intervals",
		x = expression(
			"Environmental trait separation at age 10 (" * abs(Delta * "trait") * ")"
		),
		y = expression(
			"Environmental trait separation at age 100 (" * abs(Delta * "trait") * ")"
		)
	) +
	theme_succession(base_size = 9) +
	theme(
		legend.position = "bottom",
		legend.box = "vertical",
		legend.box.just = "center",
		plot.title = element_text(hjust = 0.5),
		plot.subtitle = element_text(hjust = 0.5, colour = "grey35", size = 8),
		panel.grid.minor = element_blank()
	) +
	guides(
		colour = guide_legend(order = 1, nrow = 1),
		shape = guide_legend(order = 2, nrow = 1)
	)

# -- Panel b: magnitude of slope modification ---------------------------------

panel_b_data <- pdp_summary
panel_b_summary <- panel_b_data %>%
	group_by(leaf_label, variable_label) %>%
	summarise(
		median_magnitude = median(modification_magnitude_100yr),
		n_supported = sum(slope_direction_stable),
		n_traits = n(),
		.groups = "drop"
	)

panel_b_ymax <- max(panel_b_data$modification_magnitude_100yr, na.rm = TRUE)

panel_b <- ggplot(
	panel_b_data,
	aes(x = variable_label, y = modification_magnitude_100yr)
) +
	geom_point(
		aes(alpha = slope_direction_stable, colour = leaf_label),
		position = position_jitter(width = 0.10, height = 0, seed = 42),
		size = 2.2
	) +
	geom_point(
		data = panel_b_summary,
		aes(y = median_magnitude),
		shape = 23, size = 3.7, stroke = 0.55,
		fill = "black", colour = "white"
	) +
	geom_text(
		data = panel_b_summary,
		aes(
			x = variable_label,
			y = 1.10 * panel_b_ymax,
			label = paste0(n_supported, "/", n_traits)
		),
		inherit.aes = FALSE,
		size = 2.75, colour = "grey25"
	) +
	scale_colour_manual(values = COLS_LEAFTYPE, guide = "none") +
	scale_alpha_manual(
		values = c("TRUE" = 0.90, "FALSE" = 0.22),
		guide = "none"
	) +
	scale_x_discrete(guide = guide_axis(angle = 28)) +
	scale_y_continuous(expand = expansion(mult = c(0.02, 0.17))) +
	coord_cartesian(ylim = c(0, 1.14 * panel_b_ymax), clip = "off") +
	facet_wrap(~ leaf_label, nrow = 1) +
	labs(
		title = "Environmental modification of successional rates is stronger in broadleaf forests",
		subtitle = "Points are traits; diamonds are medians; labels give signed slope-difference intervals excluding zero",
		x = NULL,
		y = expression(
			"Slope modification magnitude (" * abs(Delta * "slope") *
				", trait SD per 100 years)"
		)
	) +
	theme_succession(base_size = 9) +
	theme(
		legend.position = "none",
		plot.title = element_text(hjust = 0.5),
		plot.subtitle = element_text(hjust = 0.5, colour = "grey35", size = 8),
		panel.grid.major.x = element_blank(),
		panel.grid.minor = element_blank()
	)

write_csv(panel_a_data,
	file.path(DIAGNOSTIC_DIR, "figure3_panel_a_data.csv"))
write_csv(panel_b_data,
	file.path(DIAGNOSTIC_DIR, "figure3_panel_b_data.csv"))
write_csv(panel_b_summary,
	file.path(DIAGNOSTIC_DIR, "figure3_panel_b_summary.csv"))

fig3_grouped <- panel_a / panel_b +
	plot_layout(heights = c(1.35, 0.95)) +
	plot_annotation(tag_levels = "a")

png_review <- file.path(OUT_DIR, "fig3_pdp_grouped_review.png")
pdf_review <- file.path(OUT_DIR, "fig3_pdp_grouped_review.pdf")
ggsave(
	png_review, fig3_grouped,
	width = 180, height = 220, units = "mm", dpi = 400, bg = "white"
)
ggsave(
	pdf_review, fig3_grouped,
	width = 180, height = 220, units = "mm", device = cairo_pdf, bg = "white"
)

message("Grouped Figure 3 review version saved:")
message("  ", png_review)
message("  ", pdf_review)
