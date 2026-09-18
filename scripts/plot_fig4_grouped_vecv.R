################################################################################
## succession_traits: review figure for grouped VEcv results
##
## Panel a: predictive skill through succession, averaged with equal weight
##          across traits after averaging their ten environmental strata.
## Panel b: absolute high-low environmental difference in predictive skill,
##          averaged with equal weight across traits for each gradient.
##
## Ribbons are empirical 95% repeated-partition intervals. They quantify
## sensitivity to grouped CV partitions, not population sampling uncertainty.
## The common age window (bin midpoints 15-125 years) prevents sparse edge bins
## from driving the visual summary.
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

RUN_DIR <- "tables/vecv_grouped_pid/normal_pid_f10_r30_seed42"
RAW_PATH <- file.path(RUN_DIR, "vecv_raw.rds")
METADATA_PATH <- file.path(RUN_DIR, "run_metadata.rds")
OUT_DIR <- "figures/main"
DIAGNOSTIC_DIR <- "tables/diagnostics/submission_audit_grouped_vecv"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DIAGNOSTIC_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(RAW_PATH) || !file.exists(METADATA_PATH)) {
	stop(
		"Grouped production output is missing. Run scripts/06_vecv.R on the server first.",
		call. = FALSE
	)
}

metadata <- read_rds(METADATA_PATH)
if (!identical(metadata$config_id, "normal_pid_f10_r30_seed42") ||
	!identical(metadata$resampling_unit, "PID") ||
	!identical(metadata$n_repeats, 30L)) {
	stop("Grouped VEcv metadata do not match the production configuration.",
		call. = FALSE)
}

q025 <- function(x) unname(quantile(x, 0.025, na.rm = TRUE, type = 8))
q975 <- function(x) unname(quantile(x, 0.975, na.rm = TRUE, type = 8))

vecv_raw <- read_rds(RAW_PATH) %>%
	filter(between(standage_mid, 15, 125)) %>%
	mutate(
		leaf_label = recode(leaf_type, !!!LEAFTYPE_LABELS),
		variable_label = recode(variable, !!!ENV_LABELS) %>%
			str_remove(" PC$")
	)

# Each trait contributes equally. Environmental strata are averaged within a
# trait before the nine trait means are averaged within each repeat.
panel_a_repeat <- vecv_raw %>%
	group_by(leaf_type, leaf_label, trait, repeat_id, standage_mid) %>%
	summarise(
		VEcv = mean(VEcv),
		n_context_cells = n(),
		.groups = "drop"
	) %>%
	filter(n_context_cells == 10L) %>%
	group_by(leaf_type, leaf_label, repeat_id, standage_mid) %>%
	summarise(VEcv = mean(VEcv), n_traits = n_distinct(trait), .groups = "drop") %>%
	filter(n_traits == 9L)

panel_a_data <- panel_a_repeat %>%
	group_by(leaf_type, leaf_label, standage_mid) %>%
	summarise(
		VEcv_med = median(VEcv),
		VEcv_lwr = q025(VEcv),
		VEcv_upr = q975(VEcv),
		.groups = "drop"
	)

# Pair upper and lower environmental strata within the same trait and repeat,
# then take absolute differences before aggregating across traits. This avoids
# cancellation among traits with opposing signed differences.
panel_b_repeat <- vecv_raw %>%
	select(
		trait, leaf_type, leaf_label, variable, variable_label,
		standage_mid, repeat_id, env_group, VEcv
	) %>%
	pivot_wider(names_from = env_group, values_from = VEcv) %>%
	filter(is.finite(high), is.finite(low)) %>%
	mutate(abs_delta_VEcv = abs(high - low)) %>%
	group_by(
		leaf_type, leaf_label, variable, variable_label,
		repeat_id, standage_mid
	) %>%
	summarise(
		abs_delta_VEcv = mean(abs_delta_VEcv),
		n_traits = n_distinct(trait),
		.groups = "drop"
	) %>%
	filter(n_traits == 9L)

panel_b_data <- panel_b_repeat %>%
	group_by(leaf_type, leaf_label, variable, variable_label, standage_mid) %>%
	summarise(
		abs_delta_med = median(abs_delta_VEcv),
		abs_delta_lwr = q025(abs_delta_VEcv),
		abs_delta_upr = q975(abs_delta_VEcv),
		.groups = "drop"
	) %>%
	mutate(
		variable_label = factor(
			variable_label,
			levels = c(
				"Temperature", "Precipitation", "Soil water retention",
				"Elevation", "Soil pH"
			)
		),
		leaf_label = factor(leaf_label, levels = names(COLS_LEAFTYPE))
	)

write_csv(panel_a_data,
	file.path(DIAGNOSTIC_DIR, "figure4_panel_a_data.csv"))
write_csv(panel_b_data,
	file.path(DIAGNOSTIC_DIR, "figure4_panel_b_data.csv"))

panel_a <- ggplot(
	panel_a_data,
	aes(x = standage_mid, colour = leaf_label, fill = leaf_label)
) +
	geom_ribbon(
		aes(ymin = VEcv_lwr, ymax = VEcv_upr),
		alpha = 0.18, colour = NA
	) +
	geom_line(aes(y = VEcv_med), linewidth = 0.9) +
	geom_point(aes(y = VEcv_med), size = 1.4) +
	facet_wrap(~ leaf_label, nrow = 1) +
	scale_colour_manual(values = COLS_LEAFTYPE, guide = "none") +
	scale_fill_manual(values = COLS_LEAFTYPE, guide = "none") +
	scale_x_continuous(
		breaks = c(25, 50, 75, 100, 125),
		expand = expansion(mult = c(0.02, 0.03))
	) +
	scale_y_continuous(
		limits = c(0, 0.55), breaks = seq(0, 0.5, 0.1),
		expand = expansion(mult = c(0, 0.03))
	) +
	labs(
		x = NULL,
		y = "Trait predictability (VEcv)",
		subtitle = "Mean predictive skill across nine traits and environmental strata"
	) +
	theme_succession(base_size = 9) +
	theme(
		plot.subtitle = element_text(size = 8.5, colour = "grey30"),
		panel.grid.minor = element_blank()
	)

panel_b <- ggplot(
	panel_b_data,
	aes(x = standage_mid, colour = leaf_label, fill = leaf_label)
) +
	geom_hline(
		yintercept = 0, linetype = "dashed", colour = "grey55",
		linewidth = 0.35
	) +
	geom_ribbon(
		aes(ymin = abs_delta_lwr, ymax = abs_delta_upr),
		alpha = 0.18, colour = NA
	) +
	geom_line(aes(y = abs_delta_med), linewidth = 0.7) +
	geom_point(aes(y = abs_delta_med), size = 0.9) +
	facet_grid(variable_label ~ leaf_label, drop = TRUE) +
	scale_colour_manual(values = COLS_LEAFTYPE, guide = "none") +
	scale_fill_manual(values = COLS_LEAFTYPE, guide = "none") +
	scale_x_continuous(
		breaks = c(25, 50, 75, 100, 125),
		expand = expansion(mult = c(0.02, 0.03))
	) +
	scale_y_continuous(
		limits = c(0, 0.40), breaks = c(0, 0.1, 0.2, 0.3, 0.4),
		expand = expansion(mult = c(0, 0.03))
	) +
	labs(
		x = LAB_STANDAGE,
		y = expression("Environmental predictability gap (mean " *
			abs(Delta * "VEcv") * ")"),
		subtitle = "Absolute upper-lower quantile difference; zero indicates identical predictive skill"
	) +
	theme_succession(base_size = 8.5) +
	theme(
		plot.subtitle = element_text(size = 8.2, colour = "grey30"),
		strip.text.y = element_text(angle = 0, face = "bold", size = 8),
		strip.text.x = element_text(face = "bold", size = 8.5),
		panel.spacing.y = unit(1.4, "mm"),
		panel.grid.minor = element_blank()
	)

fig4_grouped <- panel_a / panel_b +
	plot_layout(heights = c(0.9, 2.1)) +
	plot_annotation(tag_levels = "a")

png_path <- file.path(OUT_DIR, "fig4_vecv_grouped_review.png")
pdf_path <- file.path(OUT_DIR, "fig4_vecv_grouped_review.pdf")
ggsave(
	png_path, fig4_grouped, width = 180, height = 210,
	units = "mm", dpi = 400, bg = "white"
)
ggsave(
	pdf_path, fig4_grouped, width = 180, height = 210,
	units = "mm", device = cairo_pdf, bg = "white"
)

message("Grouped Figure 4 review version saved:")
message("  ", png_path)
message("  ", pdf_path)
