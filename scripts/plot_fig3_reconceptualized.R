################################################################################
## Figure 3 reconceptualized prototype
##
## Main message:
## Environmental context modifies successional rates, but does not generate a
## single universal pattern of functional convergence or divergence.
##
## Run from the Trait Succession project root:
##   Rscript scripts/plot_fig3_reconceptualized.R
##
## Outputs:
##   figures/main/fig3_pdp_reconceptualized.png
##   figures/main/fig3_pdp_reconceptualized.pdf
################################################################################

rm(list = ls())

required_packages <- c("tidyverse", "patchwork")
missing_packages <- required_packages[
	!vapply(required_packages, requireNamespace, quietly = TRUE,
				FUN.VALUE = logical(1))
]
if (length(missing_packages) > 0L) {
	stop("Missing required package(s): ", paste(missing_packages, collapse = ", "),
			 call. = FALSE)
}

library(tidyverse)
library(patchwork)

theme_file <- "scripts/plot_theme.R"
input_files <- c(
	"tables/pdp_raw.rds",
	"tables/pdp_summary.rds"
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

# Match the shortened labels in the current main_plots.R.
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

pdp_raw <- readr::read_rds(input_files[[1L]]) %>% strip_pc()
pdp_summary <- readr::read_rds(input_files[[2L]]) %>% strip_pc()

required_raw <- c(
	"iteration", "leaf_type", "trait", "variable", "group", "standage", "yhat"
)
required_summary <- c(
	"leaf_type", "trait", "variable", "trait_label", "variable_label",
	"slope_median", "slope_lwr", "slope_upr", "slope_robust"
)
if (!all(required_raw %in% names(pdp_raw))) {
	stop("tables/pdp_raw.rds lacks required columns.", call. = FALSE)
}
if (!all(required_summary %in% names(pdp_summary))) {
	stop("tables/pdp_summary.rds lacks required columns.", call. = FALSE)
}

LEAF_LEVELS <- c("Broadleaf", "Coniferous")
ENV_LEVELS <- c(
	"Temperature", "Precipitation", "Soil water retention", "Elevation", "Soil pH"
)

# ==============================================================================
# Panel a: change in environmental trait differentiation through succession
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
		gap_early = abs(yhat_early_high - yhat_early_low),
		gap_late = abs(yhat_late_high - yhat_late_low),
		gap_change = gap_late - gap_early
	)

gap_summary <- early_late %>%
	group_by(leaf_type, trait, variable) %>%
	summarise(
		early_age_med = median(c(age_early_high, age_early_low), na.rm = TRUE),
		gap_early_med = median(gap_early, na.rm = TRUE),
		gap_early_lwr = quantile(gap_early, 0.025, na.rm = TRUE),
		gap_early_upr = quantile(gap_early, 0.975, na.rm = TRUE),
		gap_late_med = median(gap_late, na.rm = TRUE),
		gap_late_lwr = quantile(gap_late, 0.025, na.rm = TRUE),
		gap_late_upr = quantile(gap_late, 0.975, na.rm = TRUE),
		gap_change_lwr = quantile(gap_change, 0.025, na.rm = TRUE),
		gap_change_upr = quantile(gap_change, 0.975, na.rm = TRUE),
		change_class = case_when(
			gap_change_lwr > 0 ~ "Robust widening",
			gap_change_upr < 0 ~ "Robust narrowing",
			TRUE ~ "Uncertain change"
		),
		.groups = "drop"
	) %>%
	mutate(
		leaf_type = factor(tools::toTitleCase(leaf_type), levels = LEAF_LEVELS),
		change_class = factor(
			change_class,
			levels = c("Robust widening", "Uncertain change", "Robust narrowing")
		)
	)

early_age_range <- range(gap_summary$early_age_med, na.rm = TRUE)
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
}

axis_max <- gap_summary %>%
	summarise(across(
		c(gap_early_lwr, gap_early_upr, gap_late_lwr, gap_late_upr),
		~ max(., na.rm = TRUE)
	)) %>%
	unlist() %>%
	max() %>%
	{ ceiling(. * 10) / 10 }

CHANGE_COLOURS <- c(
	"Robust widening" = "#D95F02",
	"Uncertain change" = "grey70",
	"Robust narrowing" = "#4575B4"
)

region_labels <- tidyr::expand_grid(
	leaf_type = factor(LEAF_LEVELS, levels = LEAF_LEVELS),
	region = c("widen", "narrow")
) %>%
	mutate(
		x = if_else(region == "widen", 0.06, 0.94) * axis_max,
		y = if_else(region == "widen", 0.90, 0.08) * axis_max,
		hjust = if_else(region == "widen", 0, 1),
		label = if_else(
			region == "widen",
			"Environmental differences\nwiden",
			"Environmental differences\nnarrow"
		)
	)

panel_a <- ggplot(gap_summary, aes(gap_early_med, gap_late_med)) +
	geom_abline(slope = 1, intercept = 0, colour = "grey30", linewidth = 0.6) +
	geom_text(
		data = region_labels,
		aes(x = x, y = y, label = label, hjust = hjust),
		inherit.aes = FALSE,
		colour = "grey35",
		fontface = "italic",
		size = 2.8
	) +
	geom_errorbar(
		aes(ymin = gap_late_lwr, ymax = gap_late_upr),
		colour = "grey65", linewidth = 0.3, width = 0.018 * axis_max
	) +
	geom_errorbar(
		aes(
			xmin = gap_early_lwr, xmax = gap_early_upr,
			ymin = gap_late_med, ymax = gap_late_med
		),
		colour = "grey65", linewidth = 0.3,
		width = 0.018 * axis_max, orientation = "y"
	) +
	geom_point(aes(colour = change_class), size = 2.7) +
	scale_colour_manual(
		values = CHANGE_COLOURS,
		name = "Change through succession",
		drop = FALSE
	) +
	coord_equal(xlim = c(0, axis_max), ylim = c(0, axis_max)) +
	facet_wrap(~ leaf_type, ncol = 2) +
	labs(
		title = "Do environmental differences widen or narrow through succession?",
		x = early_axis_label,
		y = expression(
			"Environmental trait differentiation after 100 years (" *
				"|" * Delta * "trait|" * ")"
		)
	) +
	theme_succession(base_size = 9) +
	theme(
		legend.position = "bottom",
		plot.title = element_text(hjust = 0.5),
		panel.grid.minor = element_blank()
	)

# ==============================================================================
# Panel b: how strongly each environmental axis modifies successional rates
# ==============================================================================

slope_data <- pdp_summary %>%
	mutate(
		leaf_type = factor(tools::toTitleCase(leaf_type), levels = LEAF_LEVELS),
		variable_label = factor(variable_label, levels = ENV_LEVELS),
		# Traits are z-scored; this is the absolute difference between the two
		# signed slopes, expressed as trait SD per 100 years.
		modification_magnitude = 100 * abs(slope_median)
	) %>%
	filter(!is.na(variable_label), is.finite(modification_magnitude))

slope_group_summary <- slope_data %>%
	group_by(leaf_type, variable_label) %>%
	summarise(
		median_magnitude = median(modification_magnitude),
		n_robust = sum(slope_robust, na.rm = TRUE),
		n_traits = n(),
		.groups = "drop"
	)

y_max <- max(slope_data$modification_magnitude, na.rm = TRUE)

panel_b <- ggplot(
	slope_data,
	aes(x = variable_label, y = modification_magnitude)
) +
	geom_point(
		aes(colour = leaf_type, alpha = slope_robust),
		position = position_jitter(width = 0.09, height = 0, seed = 42),
		size = 2.2
	) +
	geom_point(
		data = slope_group_summary,
		aes(y = median_magnitude),
		shape = 23,
		size = 3.5,
		stroke = 0.6,
		fill = "black",
		colour = "white"
	) +
	geom_text(
		data = slope_group_summary,
		aes(
			x = variable_label,
			y = 1.10 * y_max,
			label = paste0(n_robust, "/", n_traits)
		),
		inherit.aes = FALSE,
		size = 2.8,
		colour = "grey30"
	) +
	scale_colour_manual(values = COLS_LEAFTYPE, guide = "none") +
	scale_alpha_manual(values = c("TRUE" = 0.85, "FALSE" = 0.25), guide = "none") +
	scale_x_discrete(guide = guide_axis(angle = 30)) +
	scale_y_continuous(
		expand = expansion(mult = c(0.02, 0.16))
	) +
	coord_cartesian(ylim = c(0, 1.13 * y_max), clip = "off") +
	facet_wrap(~ leaf_type, ncol = 2) +
	labs(
		title = "Where does environmental context most strongly alter trait change?",
		subtitle = "Points are traits; diamonds are medians; labels give robust differences out of nine traits",
		x = NULL,
		y = expression(
			"Modification magnitude (" * "|" * Delta * "slope|" *
				", trait SD per 100 years)"
		)
	) +
	theme_succession(base_size = 9) +
	theme(
		legend.position = "none",
		plot.title = element_text(hjust = 0.5),
		plot.subtitle = element_text(hjust = 0.5, colour = "grey35"),
		panel.grid.major.x = element_blank()
	)

fig3_reconceptualized <- panel_a / panel_b +
	plot_layout(heights = c(1.35, 1)) +
	plot_annotation(tag_levels = "a")

output_dir <- "figures/main"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
png_path <- file.path(output_dir, "fig3_pdp_reconceptualized.png")
pdf_path <- file.path(output_dir, "fig3_pdp_reconceptualized.pdf")

ggsave(
	png_path,
	fig3_reconceptualized,
	width = 180,
	height = 210,
	units = "mm",
	dpi = 400,
	bg = "white"
)
ggsave(
	pdf_path,
	fig3_reconceptualized,
	width = 180,
	height = 210,
	units = "mm",
	device = grDevices::pdf,
	bg = "white"
)

message("Reconceptualized Figure 3 saved without overwriting existing figures:")
message("  ", png_path)
message("  ", pdf_path)
