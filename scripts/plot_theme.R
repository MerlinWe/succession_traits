################################################################################
## succession_traits: plot_theme.R
## Shared aesthetics, colour palettes, and label mappings.
## Source at the top of main_plots.R and supp_plots.R.
## Do not source from pipeline scripts (03–06).
##
## Contents:
##   1. ggplot2 base theme
##   2. Colour palettes (colour-blind safe)
##   3. Label mappings
##   4. Shape mappings
##   5. Convenience save function
##
## Author: M. Weiss @ Maynard Lab UCL / ETH Zürich
################################################################################

library(tidyverse)

# ══════════════════════════════════════════════════════════════════════════════
# 1. Base ggplot theme
# ══════════════════════════════════════════════════════════════════════════════
# All publication figures use this theme as a base.
# Individual figures may override specific elements as needed.

theme_succession <- function(base_size = 9, ...) {
	theme_bw(base_size = base_size) %+replace%
		theme(
			# Strips
			strip.text = element_text(face = "bold", size = base_size,
																margin = margin(t = 4, b = 4, l = 4, r = 4)),
			strip.background   = element_rect(fill = "white", colour = "black",
																				linewidth = 0.5),
			# Grid
			panel.grid.minor   = element_blank(),
			panel.grid.major   = element_line(colour = "grey92", linewidth = 0.3),
			panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
			# Axes
			axis.text          = element_text(size = base_size - 1),
			axis.title         = element_text(size = base_size),
			# Legend
			legend.background  = element_rect(fill = "white", colour = NA),
			legend.key         = element_rect(fill = "white", colour = NA),
			legend.text        = element_text(size = base_size - 1),
			legend.title       = element_text(size = base_size, face = "bold"),
			legend.position    = "bottom",
			# Plot
			plot.title         = element_text(face = "bold", size = base_size,
																				hjust = 0),
			plot.subtitle      = element_text(size = base_size - 1, hjust = 0,
																				colour = "grey40"),
			plot.margin        = margin(5, 8, 5, 8),
			...
		)
}

# ══════════════════════════════════════════════════════════════════════════════
# 2. Colour palettes (colour-blind safe)
# ══════════════════════════════════════════════════════════════════════════════

# ── Forest types ──────────────────────────────────────────────────────────────
# Established in pipeline; kept consistent across all figures
COLS_LEAFTYPE <- c(
	"Broadleaf"  = "#228B22",   # forest green
	"Coniferous" = "#D95F02"    # burnt orange
)

# ── Environmental quantile groups ─────────────────────────────────────────────
COLS_ENVGROUP <- c(
	"Lower environmental quantile" = "#4575B4",  # cool blue
	"Upper environmental quantile" = "#D73027"   # warm red
)

# ── Environmental predictors (Okabe-Ito palette — colour-blind safe) ──────────
# Full 8-colour Okabe-Ito: black, orange, sky blue, green, yellow, blue, red, grey
# We use 5 colours for the 5 environmental predictors
COLS_ENV <- c(
	"Temperature PC"          = "#E69F00",  # orange
	"Precipitation PC"        = "#56B4E9",  # sky blue
	"Soil water retention PC" = "#009E73",  # green
	"Elevation"               = "#CC79A7",  # mauve/pink
	"Soil pH"                 = "#0072B2"   # blue
)

# ── Stand age (successional predictor) ───────────────────────────────────────
COL_STANDAGE <- "#8B4513"  # saddle brown — distinct from environmental palette

# ── All predictors combined (for stacked bar) ────────────────────────────────
COLS_PREDICTORS <- c(
	"Stand age"               = COL_STANDAGE,
	COLS_ENV
)

# ── Divergence direction (signed Δ) ──────────────────────────────────────────
COL_POS <- "#D73027"   # warm red  — high env has steeper/higher trajectory
COL_NEG <- "#4575B4"   # cool blue — low env has steeper/higher trajectory
COL_NS  <- "grey80"    # non-robust differences

# ══════════════════════════════════════════════════════════════════════════════
# 3. Label mappings
# ══════════════════════════════════════════════════════════════════════════════

TRAIT_LABELS <- c(
	"bark_thickness"     = "Bark Thickness",
	"conduit_diam"       = "Conduit Diameter",
	"height"             = "Tree Height",
	"leaf_density"       = "Leaf Density",
	"leaf_k"             = "Leaf Potassium",
	"root_depth"         = "Root Depth",
	"seed_dry_mass"      = "Seed Dry Mass",
	"shade_tolerance"    = "Shade Tolerance",
	"specific_leaf_area" = "Specific Leaf Area"
)

ENV_LABELS <- c(
	"temp_pc"   = "Temperature PC",
	"soil_pc"   = "Soil water retention PC",
	"rain_pc"   = "Precipitation PC",
	"elevation" = "Elevation",
	"soil_ph"   = "Soil pH"
)

ALL_PREDICTOR_LABELS <- c(
	"standage" = "Stand age",
	ENV_LABELS
)

LEAFTYPE_LABELS <- c(
	"broadleaf"  = "Broadleaf",
	"coniferous" = "Coniferous"
)

# Functional group labels for trait grouping in figures
TRAIT_GROUPS <- tibble::tribble(
	~trait,              ~group,
	"specific_leaf_area", "Leaf economics",
	"leaf_k",             "Leaf economics",
	"leaf_density",       "Leaf economics",
	"conduit_diam",       "Hydraulics",
	"bark_thickness",     "Bark & structure",
	"height",             "Architecture",
	"root_depth",         "Belowground",
	"seed_dry_mass",      "Reproduction",
	"shade_tolerance",    "Light competition"
)

# ══════════════════════════════════════════════════════════════════════════════
# 4. Shape mappings
# ══════════════════════════════════════════════════════════════════════════════
# 9 traits need 9 distinguishable filled shapes.
# Using filled shapes (21–25) allows fill + border colour independently.
# First 5 shapes repeat with a second set to cover 9 traits.

SHAPES_TRAITS <- c(21, 22, 23, 24, 25, 21, 22, 23, 24)
names(SHAPES_TRAITS) <- TRAIT_LABELS

# Forest types
SHAPES_LEAFTYPE <- c("Broadleaf" = 16, "Coniferous" = 17)

# ══════════════════════════════════════════════════════════════════════════════
# 5. Convenience save function
# ══════════════════════════════════════════════════════════════════════════════
# Wraps ggsave with consistent defaults.
# dpi = 400 for main figures, 300 for supplementary.

save_fig <- function(plot, filename, dir = "figures/main",
										 width = 180, height = 120,
										 units = "mm", dpi = 400) {
	path <- file.path(dir, filename)
	ggsave(
		filename = path,
		plot     = plot,
		width    = width,
		height   = height,
		units    = units,
		dpi      = dpi,
		bg       = "white"
	)
	message(sprintf("  ✓ Saved: %s", path))
}

# ══════════════════════════════════════════════════════════════════════════════
# 6. Shared axis / legend labels (used across multiple figures)
# ══════════════════════════════════════════════════════════════════════════════

LAB_STANDAGE   <- "Stand age (years)"
LAB_VECV       <- "VEcv (cross-validated variance explained)"
LAB_DELTA_VECV <- expression(Delta * "VEcv (upper \u2212 lower quantile)")
LAB_SHAP       <- "Mean |SHAP| importance"
LAB_SLOPE_DIFF <- expression(Delta * " slope (upper \u2212 lower quantile)")
LAB_INT_DIFF   <- expression(Delta * " intercept (upper \u2212 lower quantile)")
LAB_RATIO      <- "Environmental / successional importance (ratio, log\u2082 scale)"

message("plot_theme.R loaded.")