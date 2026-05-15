################################################################################
## succession_traits: main_plots.R
## Builds all main publication figures from pre-computed analysis outputs.
## Reads .rds files from tables/ and writes figures to figures/main/.
##
## Figures produced:
##   fig1_map.png        — Data figure (plot density + stand age maps)
##   fig2_shap.png       — RQ1: SHAP stacked bar ordered by stand age importance
##   fig3_pdp.png        — RQ2: SHAP-weighted ellipse + direction dotplot
##   fig4_vecv.png       — RQ3: VEcv trajectory subset + ΔVEcv by env variable
##
## Run after all pipeline scripts (01–06) have completed successfully.
## Requires: tables/shap_per_var.rds, shap_importance.rds, shap_values.rds,
##           pdp_stats_full (from pdp_stats.rds + pdp_summary.rds),
##           vecv_raw.rds, vecv_summary.rds, vecv_divergence.rds
##
## Author: M. Weiss @ Maynard Lab UCL / ETH Zürich
################################################################################

rm(list = ls())

# ── Libraries ─────────────────────────────────────────────────────────────────
library(tidyverse)
library(patchwork)
library(ggbeeswarm)

source("scripts/plot_theme.R")

# ── Output directory ──────────────────────────────────────────────────────────
DIR_MAIN <- "figures/main"
dir.create(DIR_MAIN, recursive = TRUE, showWarnings = FALSE)

# ── Load all tables ───────────────────────────────────────────────────────────
shap_per_var    <- read_rds("tables/shap_per_var.rds")
shap_importance <- read_rds("tables/shap_importance.rds")
shap_values     <- read_rds("tables/shap_values.rds")
pdp_raw         <- read_rds("tables/pdp_raw.rds")
pdp_stats       <- read_rds("tables/pdp_stats.rds")
pdp_summary     <- read_rds("tables/pdp_summary.rds")
vecv_raw        <- read_rds("tables/vecv_raw.rds")
vecv_summary    <- read_rds("tables/vecv_summary.rds")
vecv_divergence <- read_rds("tables/vecv_divergence.rds")

# ── Analytical constants ──────────────────────────────────────────────────────
ENV_VARS   <- c("temp_pc", "soil_pc", "rain_pc", "elevation", "soil_ph")
LEAF_TYPES <- c("broadleaf", "coniferous")
TRAITS     <- names(TRAIT_LABELS)


################################################################################
## Figure 1 — Data map
## Layout:
##   Row 1: Plot sampling intensity (a) | Median stand age (b)
##   Row 2: Stand age ridge (full width, c)
##   Row 3: Forest type map + stand age by forest type density (d) |
##          Temperature PC × elevation environmental space (e)
################################################################################

message("Building Figure 1 (data map)...")

library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(scales)
library(ggridges)
library(patchwork)
library(grid)

# ── Load data ─────────────────────────────────────────────────────────────────
dat_map <- read_rds("data_processed/fia_traits_clean.rds") %>%
	dplyr::select(PID_rep, standage, LAT, LON, temp_pc, elevation,
								biome_boreal_forests_or_taiga,
								biome_temperate_conifer_forests,
								biome_temperate_broadleaf_forests,
								biome_mediterranean_woodlands) %>%
	filter(complete.cases(.)) %>%
	mutate(
		leaf_type = case_when(
			biome_boreal_forests_or_taiga   == 1 |
				biome_temperate_conifer_forests == 1 ~ "Coniferous",
			biome_temperate_broadleaf_forests == 1 |
				biome_mediterranean_woodlands   == 1 ~ "Broadleaf",
			TRUE ~ NA_character_
		)
	) %>%
	filter(!is.na(leaf_type))

message(sprintf("  Map data: %d plots (%d broadleaf / %d coniferous)",
								nrow(dat_map),
								sum(dat_map$leaf_type == "Broadleaf"),
								sum(dat_map$leaf_type == "Coniferous")))

# ── Projections and base geography ────────────────────────────────────────────
crs_conus <- 5070L
crs_ak    <- 3338L

usa_wgs <- ne_countries(
	scale = "medium", country = "united states of america", returnclass = "sf"
)

conus_bbox_wgs <- st_as_sfc(
	st_bbox(c(xmin = -125, xmax = -66.5, ymin = 24.5, ymax = 49.5),
					crs = st_crs(4326))
)
usa_conus <- st_transform(st_intersection(usa_wgs, conus_bbox_wgs), crs_conus)

ak_bbox_wgs <- st_as_sfc(
	st_bbox(c(xmin = -180, xmax = -130, ymin = 50, ymax = 72),
					crs = st_crs(4326))
)
usa_ak <- st_transform(st_intersection(usa_wgs, ak_bbox_wgs), crs_ak)

pts_wgs   <- st_as_sf(dat_map, coords = c("LON", "LAT"), crs = 4326L)
pts_conus <- pts_wgs %>% st_transform(crs_conus) %>%
	st_intersection(st_transform(conus_bbox_wgs, crs_conus))
pts_ak    <- pts_wgs %>% st_transform(crs_ak) %>%
	st_intersection(st_transform(ak_bbox_wgs, crs_ak))

# ── Hex aggregation ───────────────────────────────────────────────────────────
hex_aggregate <- function(points_sf, region_sf,
													cell_km = 100, min_n_for_age = 3) {
	cellsize_m <- units::drop_units(
		units::set_units(units::set_units(cell_km, "km"), "m")
	)
	hex <- st_make_grid(region_sf, cellsize = cellsize_m,
											square = FALSE, what = "polygons") %>%
		st_sf(geometry = .) %>%
		st_make_valid() %>%
		st_intersection(st_make_valid(st_union(region_sf))) %>%
		st_collection_extract("POLYGON") %>%
		filter(st_is_valid(.)) %>%
		mutate(.area = as.numeric(st_area(.))) %>%
		filter(.area > 0) %>%
		dplyr::select(-.area) %>%
		mutate(id = row_number())
	
	j <- st_join(points_sf, hex, join = st_within, left = FALSE)
	
	agg <- j %>%
		st_drop_geometry() %>%
		group_by(id) %>%
		summarise(
			n              = n(),
			median_age     = median(standage,  na.rm = TRUE),
			median_temp_pc = median(temp_pc,   na.rm = TRUE),
			median_elev    = median(elevation, na.rm = TRUE),
			# dominant leaf type per hex — majority vote
			leaf_type      = {
				lt <- leaf_type
				names(which.max(table(lt[!is.na(lt)])))
			},
			.groups = "drop"
		)
	
	hex %>%
		left_join(agg, by = "id") %>%
		mutate(
			median_age_shown = if_else(!is.na(n) & n >= min_n_for_age,
																 median_age, NA_real_),
			temp_pc_shown    = if_else(!is.na(n) & n >= min_n_for_age,
																 median_temp_pc, NA_real_),
			elev_shown       = if_else(!is.na(n) & n >= min_n_for_age,
																 median_elev, NA_real_),
			leaf_type_shown  = if_else(!is.na(n) & n >= min_n_for_age,
																 leaf_type, NA_character_)
		)
}

message("  Building hex grids...")
conus_hex <- hex_aggregate(pts_conus, usa_conus, cell_km = 100, min_n_for_age = 3)
ak_hex    <- hex_aggregate(pts_ak,    usa_ak,    cell_km = 150, min_n_for_age = 3)

# ── Shared theme and scales ───────────────────────────────────────────────────
theme_map_panel <- function(base_size = 8) {
	theme_void(base_size = base_size) +
		theme(
			legend.position   = "top",
			legend.title      = element_text(size = base_size, face = "bold"),
			legend.text       = element_text(size = base_size - 1),
			legend.key.height = unit(2.5, "mm"),
			legend.key.width  = unit(8,   "mm"),
			plot.margin       = margin(1, 1, 1, 1)
		)
}

scale_count <- scale_fill_viridis_c(
	"Plots per hex",
	trans    = "sqrt",
	breaks   = pretty_breaks(4),
	na.value = "grey85",
	option   = "viridis"
)

scale_age <- scale_fill_gradientn(
	"Median stand age",
	colours  = c("#2C1654", "#7B2D8B", "#D95F02", "#FFC107", "#FFFDE7"),
	na.value = "grey85",
	breaks   = c(40, 60, 80, 100)
)

# Leaf type colour scale — consistent with COLS_LEAFTYPE from plot_theme.R
scale_leaftype <- scale_fill_manual(
	name     = "Forest type",
	values   = c("Broadleaf"  = "#228B22",
							 "Coniferous" = "#D95F02"),
	na.value = "grey85",
	guide    = guide_legend(
		direction      = "horizontal",
		title.position = "top",
		override.aes   = list(colour = NA)
	)
)

# ── Panel builders ──────────────────────────────────
build_panel <- function(fill_var, hex_conus, hex_ak, fill_scale) {
	
	p_main <- ggplot() +
		geom_sf(data = hex_conus,
						aes(fill = .data[[fill_var]]), colour = NA) +
		geom_sf(data = usa_conus,
						fill = NA, colour = "black", linewidth = 0.2) +
		fill_scale +
		theme_map_panel()
	
	p_ak <- ggplot() +
		geom_sf(data = hex_ak,
						aes(fill = .data[[fill_var]]), colour = NA) +
		geom_sf(data = usa_ak,
						fill = NA, colour = "black", linewidth = 0.15) +
		fill_scale +
		theme_void(6) +
		theme(legend.position = "none", plot.background = element_blank())
	
	p_main +
		inset_element(p_ak,
									left = 0.0, bottom = 0.0, right = 0.22, top = 0.28,
									align_to = "plot", ignore_tag = TRUE)
}

# Categorical panel builder — used for forest type map
build_cat_panel <- function(fill_var, hex_conus, hex_ak, fill_scale) {
	
	p_main <- ggplot() +
		geom_sf(data = hex_conus %>% filter(!is.na(.data[[fill_var]])),
						aes(fill = .data[[fill_var]]), colour = NA) +
		geom_sf(data = hex_conus %>% filter(is.na(.data[[fill_var]])),
						fill = "grey85", colour = NA) +
		geom_sf(data = usa_conus,
						fill = NA, colour = "black", linewidth = 0.2) +
		fill_scale +
		theme_map_panel()
	
	p_ak <- ggplot() +
		geom_sf(data = ak_hex %>% filter(!is.na(.data[[fill_var]])),
						aes(fill = .data[[fill_var]]), colour = NA) +
		geom_sf(data = ak_hex %>% filter(is.na(.data[[fill_var]])),
						fill = "grey85", colour = NA) +
		geom_sf(data = usa_ak,
						fill = NA, colour = "black", linewidth = 0.15) +
		fill_scale +
		theme_void(6) +
		theme(legend.position = "none", plot.background = element_blank())
	
	p_main +
		inset_element(p_ak,
									left = 0.0, bottom = 0.0, right = 0.22, top = 0.28,
									align_to = "plot", ignore_tag = TRUE)
}

# ── Build map panels ──────────────────────────────────────────────────────────
message("  Building map panels...")

p_counts   <- build_panel("n",               conus_hex, ak_hex, scale_count)
p_age      <- build_panel("median_age_shown", conus_hex, ak_hex, scale_age)
p_leaftype <- build_cat_panel("leaf_type_shown", conus_hex, ak_hex,
															scale_leaftype)

# ── Stand age ridge — full width, coloured by stand age gradient ──────────────
age_vals <- dat_map$standage
dens     <- density(age_vals, na.rm = TRUE, bw = "nrd0", n = 512)
df_dens  <- tibble(age = dens$x, dens = dens$y) %>%
	mutate(height = scales::rescale(dens, to = c(0, 1)))

p_ridge <- ggplot(df_dens, aes(x = age, y = 0, height = height, fill = age)) +
	geom_ridgeline_gradient(colour = "grey30", linewidth = 0.25, alpha = 0.95) +
	scale_fill_gradientn(
		colours = c("#2C1654", "#7B2D8B", "#D95F02", "#FFC107", "#FFFDE7"),
		guide   = "none"
	) +
	scale_x_continuous(limits = c(0, 150), breaks = c(0, 50, 100, 150)) +
	scale_y_continuous(expand = expansion(mult = c(0.1, 0.15))) +
	labs(x = "Stand age (years)", y = NULL) +
	theme_succession(base_size = 9) +
	theme(
		panel.grid     = element_blank(),
		panel.grid.major.x = element_blank(),
		panel.grid.minor.x = element_blank(),
		panel.grid.major.y = element_blank(),
		panel.grid.minor.y = element_blank(),
		axis.text.y    = element_blank(),
		axis.ticks.y   = element_blank(),
		panel.border   = element_blank(),
		plot.margin    = margin(0, 2, 1, 2)
	)

# ── Stand age by forest type density ─────────────────────────────────────────
# Separate density curves per forest type — shows successional range overlap
p_age_by_lt <- dat_map %>%
	ggplot(aes(x = standage, colour = leaf_type, fill = leaf_type)) +
	geom_density(alpha = 0.25, linewidth = 0.6) +
	scale_colour_manual(values = COLS_LEAFTYPE, name = NULL) +
	scale_fill_manual(  values = COLS_LEAFTYPE, name = NULL) +
	scale_x_continuous(limits = c(0, 150), breaks = c(0, 50, 100, 150)) +
	labs(x = "Stand age (years)", y = NULL) +
	theme_succession(base_size = 8) +
	theme(
		legend.position = "none",
		panel.grid.major.x = element_blank(),
		panel.grid.minor.x = element_blank(),
		panel.grid.major.y = element_blank(),
		panel.grid.minor.y = element_blank(),
		axis.text.y    = element_blank(),
		axis.ticks.y   = element_blank(),
		panel.border    = element_blank(),
		plot.margin     = margin(0, 2, 2, 2)
	)

# ── Environmental space: temperature PC × elevation, coloured by forest type ──
# 2D density contours in feature space — no geographic boundary issues
# Thin point layer underneath for texture, contours on top for clarity
p_env_space <- dat_map %>%
	ggplot(aes(x = temp_pc, y = elevation, colour = leaf_type,
						 fill = leaf_type)) +
	geom_point(alpha = 0.04, size = 0.3) +
	stat_density_2d(
		aes(fill = leaf_type),
		geom     = "polygon",
		alpha    = 0.25,
		colour   = NA,
		contour_var = "ndensity",   # normalised so both types visible equally
		bins     = 6
	) +
	stat_density_2d(
		aes(colour = leaf_type),
		geom      = "contour",
		linewidth = 0.4,
		contour_var = "ndensity",
		bins      = 6
	) +
	scale_colour_manual(values = COLS_LEAFTYPE, name = "Forest type") +
	scale_fill_manual(  values = COLS_LEAFTYPE, name = "Forest type") +
	labs(
		x = "Temperature PC",
		y = "Elevation (m)"
	) +
	theme_succession(base_size = 8) +
	theme(
		legend.position = "bottom",
		legend.key.size = unit(3, "mm")
	)

# ── Assemble ──────────────────────────────────────────────────────────────────
# Row 1: plot density | stand age map
# Row 2: stand age ridge (full width)
# Row 3: forest type map + age density | env space scatter

row_top <- (p_counts | p_age) +
	plot_layout(ncol = 2)

row_mid <- p_ridge

col_left_bot <- (p_leaftype / p_age_by_lt) +
	plot_layout(heights = c(1, 0.35))

col_right_bot <- p_env_space

row_bot <- (col_left_bot | col_right_bot) +
	plot_layout(ncol = 2)

fig1 <- (row_top / row_mid / row_bot) +
	plot_layout(heights = c(1, 0.18, 1.15)) +
	plot_annotation(
		tag_levels = NULL,
		theme      = theme(plot.margin = margin(2, 2, 2, 2))
	)

save_fig(fig1, "fig1_map.png",
				 width  = 180,
				 height = 230,
				 dpi    = 400)

message("  ✓ Figure 1 saved")

# ══════════════════════════════════════════════════════════════════════════════
# Figure 2 — RQ1: Relative importance of environmental vs successional filtering
# ══════════════════════════════════════════════════════════════════════════════
# Two-panel figure:
#   Panel A: Lollipop showing the env/succ importance ratio per trait on a
#            log2 scale. Directly answers "how much does environment dominate
#            over succession?" with a single quotable number per trait.
#            Dashed line at 1× = equal importance.
#   Panel B: Stacked bar showing % of total |SHAP| per predictor for each
#            trait. Answers "which environmental predictor drives that
#            dominance?" — the mechanistic follow-up to Panel A.
#
# Key design decisions:
#   - Traits ordered LEFT → RIGHT in Panel B by INCREASING stand age proportion
#     matching the top-to-bottom ordering in Panel A (most successional first)
#     so readers can cross-reference between panels immediately
#   - Stand age (brown) always rendered as the bottom segment in Panel B so
#     its varying height is the primary visual signal
#   - Panel widths 40/60 — lollipop needs less horizontal space than stacked bar

message("Building Figure 2 (SHAP stacked bar)...")

# Compute stand age proportion per trait × leaf type for ordering
standage_order <- shap_per_var %>%
	filter(variable == "standage") %>%
	group_by(trait_label) %>%
	summarise(mean_succ_prop = mean(sum_abs_shap / {
		shap_per_var %>%
			group_by(trait, leaf_type) %>%
			summarise(total = sum(sum_abs_shap), .groups = "drop") %>%
			group_by(trait) %>%
			summarise(total = mean(total)) %>%
			filter(trait == trait[1]) %>%
			pull(total)
	}), .groups = "drop")

# Cleaner approach: compute proportion directly
trait_totals <- shap_per_var %>%
	group_by(trait, leaf_type) %>%
	summarise(total_shap = sum(sum_abs_shap), .groups = "drop")

shap_pct <- shap_per_var %>%
	left_join(trait_totals, by = c("trait", "leaf_type")) %>%
	mutate(pct = 100 * sum_abs_shap / total_shap)

# Order traits by mean stand age proportion (ascending = succession first)
trait_order_fig2 <- shap_pct %>%
	filter(variable == "standage") %>%
	group_by(trait_label) %>%
	summarise(mean_succ_pct = mean(pct), .groups = "drop") %>%
	arrange(desc(mean_succ_pct)) %>%   # most successional first (left)
	pull(trait_label)

# ── Panel A: ratio lollipop ───────────────────────────────────────────────────
fig2a <- shap_importance %>%
	mutate(
		trait_label = factor(
			trait_label,
			levels = shap_importance %>%
				group_by(trait_label) %>%
				summarise(m = mean(env_succ_ratio), .groups = "drop") %>%
				arrange(m) %>%
				pull(trait_label)
		),
		leaf_type = tools::toTitleCase(leaf_type)
	) %>%
	ggplot(aes(x = env_succ_ratio, y = trait_label, colour = leaf_type)) +
	geom_vline(xintercept = 1, linetype = "dashed",
						 colour = "grey50", linewidth = 0.4) +
	geom_segment(aes(x = 1, xend = env_succ_ratio,
									 y = trait_label, yend = trait_label),
							 linewidth = 0.6, alpha = 0.5) +
	geom_point(size = 3.5) +
	scale_colour_manual(values = COLS_LEAFTYPE, name = NULL) +
	scale_x_continuous(
		trans  = "log2",
		breaks = c(0.5, 1, 2, 5, 10, 20, 30),
		labels = function(x) paste0(x, "\u00d7")
	) +
	facet_wrap(~ leaf_type, ncol = 1) +
	labs(x = LAB_RATIO, y = NULL) +
	theme_succession(base_size = 9) +
	theme(
		legend.position    = "none",
		panel.grid.major.y = element_blank()
	)

# ── Panel B: stacked bar ──────────────────────────────────────────────────────
fig2b <- shap_pct %>%
	mutate(
		leaf_type      = tools::toTitleCase(leaf_type),
		trait_label    = factor(trait_label, levels = trait_order_fig2),
		variable_label = recode(variable, !!!ALL_PREDICTOR_LABELS),
		variable_label = factor(variable_label, levels = names(COLS_PREDICTORS))
	) %>%
	ggplot(aes(x = trait_label, y = pct, fill = variable_label)) +
	geom_col(colour = "white", linewidth = 0.25, width = 0.8) +
	facet_wrap(~ leaf_type, ncol = 1) +
	scale_fill_manual(
		name   = NULL,
		values = COLS_PREDICTORS,
		breaks = names(COLS_PREDICTORS),
		guide  = guide_legend(nrow = 2)
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
		legend.position    = "right",
		panel.grid.major.x = element_blank()
	)

# ── Combine ───────────────────────────────────────────────────────────────────
fig2 <- (fig2a | fig2b) +
	plot_layout(widths = c(0.4, 0.6)) +
	theme(plot.margin = margin(4, 4, 4, 4))
fig2 <- fig2 + plot_annotation(tag_levels = "a")

save_fig(fig2, "fig2_shap.png", width = 240, height = 160)

# ══════════════════════════════════════════════════════════════════════════════
# Figure 3 — RQ2: Environmental modulation of successional trait trajectories
# ══════════════════════════════════════════════════════════════════════════════
# Two-panel figure:
#   Panel A: Scatter of |Δ trait| at stand age 0 (x) vs |Δ trait| at stand age
#            100 (y). One point per trait × environmental variable combination.
#            Diagonal y = x is the reference: points on the diagonal = gap
#            unchanged through succession; above = gap widened (diverging);
#            below = gap narrowed (converging). Colour = trait (viridis_d),
#            shape = environmental variable. Faded points = non-robust.
#   Panel B: Direction and robustness dotplot across all trait × env var combos.
#            Dot colour = direction (red = high env steeper, blue = low env
#            steeper), open circle = non-robust, size = |Δslope|.

message("Building Figure 3 (PDP early-late scatter + direction dotplot)...")

# ── Panel A data: |Δ trait| at early vs late succession ──────────────────────
LATE_AGE <- 100L   # stand age (years) used as "late succession" reference point

# Extract predicted yhat at standage ≈ 0 and standage ≈ LATE_AGE per bootstrap
# iteration, then compute high − low gap at each time point.
early_late <- pdp_raw %>%
	group_by(iteration, leaf_type, trait, variable, group) %>%
	summarise(
		yhat_early = yhat[which.min(abs(standage - 0))],
		yhat_late  = yhat[which.min(abs(standage - LATE_AGE))],
		.groups    = "drop"
	) %>%
	pivot_wider(
		names_from  = group,
		values_from = c(yhat_early, yhat_late)
	) %>%
	mutate(
		delta_early     = yhat_early_high - yhat_early_low,
		delta_late      = yhat_late_high  - yhat_late_low,
		# Absolute values computed per iteration — correct basis for CI
		abs_delta_early = abs(delta_early),
		abs_delta_late  = abs(delta_late)
	)

# Summarise across bootstrap iterations.
# CIs on the absolute scale are computed from per-iteration absolute values,
# which correctly handles the case where the signed CI crosses zero.
# Robustness criterion retains the original signed convention (CI of the signed
# delta_late excludes zero), consistent with the dotplot in Panel B.
early_late_summary <- early_late %>%
	group_by(leaf_type, trait, variable) %>%
	summarise(
		abs_early_med = median(abs_delta_early,          na.rm = TRUE),
		abs_early_lwr = quantile(abs_delta_early, 0.025, na.rm = TRUE),
		abs_early_upr = quantile(abs_delta_early, 0.975, na.rm = TRUE),
		abs_late_med  = median(abs_delta_late,           na.rm = TRUE),
		abs_late_lwr  = quantile(abs_delta_late,  0.025, na.rm = TRUE),
		abs_late_upr  = quantile(abs_delta_late,  0.975, na.rm = TRUE),
		robust        = (quantile(delta_late, 0.025, na.rm = TRUE) > 0 |
										 	quantile(delta_late, 0.975, na.rm = TRUE) < 0),
		.groups = "drop"
	) %>%
	mutate(
		trait_label    = recode(trait,    !!!TRAIT_LABELS),
		variable_label = recode(variable, !!!ENV_LABELS),
		leaf_type      = tools::toTitleCase(leaf_type),
		# Factor levels: traits follow TRAIT_LABELS order; env vars follow COLS_ENV
		trait_label    = factor(trait_label,    levels = TRAIT_LABELS),
		variable_label = factor(variable_label, levels = names(COLS_ENV))
	)

# Shared square axis limit — ensures diagonal is visually 45°
ax_max_pa <- early_late_summary %>%
	summarise(across(c(abs_early_lwr, abs_early_upr,
										 abs_late_lwr,  abs_late_upr),
									 ~ max(., na.rm = TRUE))) %>%
	unlist() %>% max() %>%
	{ ceiling(. * 10) / 10 }

# Solid shape mapping for 5 environmental variables
# Chosen to be clearly distinguishable at small sizes
SHAPES_ENV_VARS <- c(
	"Temperature PC"          = 16,   # filled circle
	"Soil water retention PC" = 17,   # filled triangle up
	"Precipitation PC"        = 15,   # filled square
	"Elevation"               = 18,   # filled diamond
	"Soil pH"                 = 8     # asterisk
)

panel_a <- early_late_summary %>%
	ggplot(aes(x = abs_early_med, y = abs_late_med)) +
	
	# Diagonal y = x: gap unchanged through succession
	geom_abline(slope = 1, intercept = 0,
							colour = "grey35", linewidth = 0.55) +
	
	# Dashed reference lines at zero
	geom_hline(yintercept = 0, linetype = "dashed",
						 colour = "grey60", linewidth = 0.3) +
	geom_vline(xintercept = 0, linetype = "dashed",
						 colour = "grey60", linewidth = 0.3) +
	
	# Vertical CI whiskers (uncertainty on late-succession gap)
	geom_errorbar(
		aes(ymin = abs_late_lwr, ymax = abs_late_upr),
		colour = "grey65", linewidth = 0.3, width = 0
	) +
	
	# Horizontal CI whiskers (uncertainty on early-succession gap)
	# orientation = "y" replaces deprecated geom_errorbarh (ggplot2 >= 4.0)
	geom_errorbar(
		aes(xmin = abs_early_lwr, xmax = abs_early_upr,
				ymin = abs_late_med,  ymax = abs_late_med),
		colour = "grey65", linewidth = 0.3,
		width = 0, orientation = "y"
	) +
	
	# Points: colour = trait, shape = environmental variable
	# Alpha fades non-robust points rather than removing them
	geom_point(
		aes(colour = trait_label,
				shape  = variable_label,
				alpha  = robust),
		size   = 2.8,
		stroke = 0.5
	) +
	
	# Trait colours: viridis_d (turbo) — 9 well-separated hues
	# Consistent with scale_fill_viridis_d used in the previous panel A
	scale_colour_viridis_d(
		option = "turbo",
		name   = NULL,
		guide  = guide_legend(
			nrow         = 3,
			order        = 1,
			override.aes = list(size = 2.5, shape = 16, alpha = 1)
		)
	) +
	
	# Environmental variable shapes
	scale_shape_manual(
		values = SHAPES_ENV_VARS,
		name   = NULL,
		guide  = guide_legend(
			nrow         = 2,
			order        = 2,
			override.aes = list(size = 2.5, colour = "grey25", alpha = 1)
		)
	) +
	
	# Non-robust points faded — no legend entry needed
	scale_alpha_manual(
		values = c("TRUE" = 1, "FALSE" = 0.2),
		guide  = "none"
	) +
	
	# Square coordinate system so diagonal is 45°
	coord_equal(
		xlim = c(0, ax_max_pa),
		ylim = c(0, ax_max_pa)
	) +
	
	facet_wrap(~ leaf_type, ncol = 2) +
	
	labs(
		x = expression(
			"|" * Delta * "| trait value at stand age 0 (upper \u2212 lower quantile)"
		),
		y = expression(
			"|" * Delta * "| trait value at stand age 100 (upper \u2212 lower quantile)"
		)
	) +
	theme_succession(base_size = 9) +
	theme(
		legend.position  = "bottom",
		legend.box       = "horizontal",    # trait colours and env shapes side-by-side
		legend.spacing.x = unit(6, "mm"),   # gap between the two legend blocks
		legend.key.size  = unit(3.5, "mm")
	)

# ── Panel B: Direction and robustness dotplot ─────────────────────────────────
# Rows = traits (ordered as in fig2 — most successional at top)
# Columns = environmental variables
# Dot colour = direction of Δslope:
#   red  (COL_POS) = high-env group has steeper successional trajectory
#   blue (COL_NEG) = low-env group has steeper successional trajectory
# Open circle = non-robust (95% bootstrap CI includes zero)
# Dot size = |median Δslope|
# Facet columns = leaf type

panel_b_data <- pdp_summary %>%
	mutate(
		leaf_type      = tools::toTitleCase(leaf_type),
		trait_label    = factor(trait_label, levels = rev(trait_order_fig2)),
		variable_label = factor(variable_label,
														levels = c("Temperature PC",
																			 "Soil water retention PC",
																			 "Precipitation PC",
																			 "Elevation",
																			 "Soil pH")),
		direction      = case_when(
			slope_robust & slope_median > 0 ~ "positive",
			slope_robust & slope_median < 0 ~ "negative",
			TRUE                            ~ "non-robust"
		),
		abs_slope = abs(slope_median)
	) %>%
	filter(!is.na(variable_label))

panel_b <- ggplot(panel_b_data,
									aes(x = variable_label, y = trait_label)) +
	geom_point(
		aes(size = abs_slope, colour = direction, shape = direction),
		stroke = 0.6
	) +
	scale_colour_manual(
		values = c("positive"   = COL_POS,
							 "negative"   = COL_NEG,
							 "non-robust" = COL_NS),
		labels = c("positive"   = "High env steeper",
							 "negative"   = "Low env steeper",
							 "non-robust" = "Non-robust"),
		name   = "Direction"
	) +
	scale_shape_manual(
		values = c("positive"   = 16,
							 "negative"   = 16,
							 "non-robust" = 1),
		guide  = "none"
	) +
	scale_size_continuous(
		name  = "|Δ slope|",
		range = c(1, 6),
		guide = guide_legend(override.aes = list(colour = "grey40", shape = 16))
	) +
	scale_x_discrete(guide = guide_axis(angle = 35)) +
	facet_wrap(~ leaf_type, ncol = 2) +
	labs(x = NULL, y = NULL) +
	theme_succession(base_size = 9) +
	theme(
		panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
		legend.position  = "bottom"
	)

# ── Combine panels A and B ────────────────────────────────────────────────────
# Heights: panel A is square (coord_equal) so will naturally be taller than
# panel B. Heights = c(1.4, 1) gives panel A slightly more room while keeping
# both readable. Figure height increased to 230mm to accommodate legends from
# both panels without cropping.

fig3 <- panel_a / panel_b +
	plot_layout(heights = c(1.4, 1))

save_fig(fig3, "fig3_pdp.png", width = 180, height = 230)

# ══════════════════════════════════════════════════════════════════════════════
# Figure 4 — RQ3: Divergence in trait predictability across environmental
# gradients
# ══════════════════════════════════════════════════════════════════════════════
# Single-focus figure showing ΔVEcv (upper minus lower environmental quantile)
# as a function of stand age, averaged across traits, one line per environmental
# variable, split by forest type.
#
# The overall predictability increase (VEcv rises from ~0.55 to ~0.81 in
# broadleaf through succession) is reported in the results and shown in full
# in Figure S-8. This figure focuses on the divergence finding — that
# environmental context creates persistent differences in predictability
# throughout succession — which is the novel RQ3 contribution.
#
# Key design decisions:
#   - Larger panels than previous two-panel version — divergence is the
#     sole visual argument
#   - Dashed reference line at ΔVEcv = 0 (no divergence)
#   - Broadleaf and coniferous as facet columns for direct comparison
#   - COLS_ENV palette (Okabe-Ito) for environmental variables
#   - Loess smoothing (span = 0.5) reduces noise in averaged trajectories
#   - Ribbons show mean of per-trait CIs — gives sense of cross-trait
#     variability in divergence

message("Building Figure 4 (VEcv divergence)...")
delta_avg <- vecv_divergence %>%
	filter(!is.na(delta_med)) %>%
	group_by(leaf_type, variable_label, standage_mid) %>%
	summarise(
		delta_med = mean(delta_med, na.rm = TRUE),
		delta_lwr = mean(delta_lwr, na.rm = TRUE),
		delta_upr = mean(delta_upr, na.rm = TRUE),
		.groups   = "drop"
	) %>%
	mutate(
		leaf_type      = tools::toTitleCase(leaf_type),
		variable_label = factor(variable_label, levels = names(COLS_ENV))
	)
fig4 <- ggplot(delta_avg,
							 aes(x = standage_mid, y = delta_med,
							 		colour = variable_label, fill = variable_label)) +
	# Uncertainty ribbon
	geom_smooth(aes(ymin = delta_lwr, ymax = delta_upr),
							method = "loess", span = 0.5,
							stat = "smooth", alpha = 0.12, linewidth = 0) +
	# Central line
	geom_smooth(aes(y = delta_med),
							method = "loess", span = 0.5,
							se = FALSE, linewidth = 0.9) +
	# Reference line
	geom_hline(yintercept = 0, linetype = "dashed",
						 colour = "grey40", linewidth = 0.5) +
	# Annotation showing overall direction
	annotate("label", x = 140, y = Inf, hjust = 1, vjust = 1.5,
					 label    = "Upper quantile\nmore predictable",
					 size     = 2.5, colour = "grey50",
					 fill     = "white", label.size = 0.3,
					 label.padding = unit(0.15, "lines")) +
	annotate("label", x = 140, y = -Inf, hjust = 1, vjust = -0.5,
					 label    = "Lower quantile\nmore predictable",
					 size     = 2.5, colour = "grey50",
					 fill     = "white", label.size = 0.3,
					 label.padding = unit(0.15, "lines")) +
	facet_wrap(~ leaf_type, ncol = 2) +
	scale_colour_manual(values = COLS_ENV, name = NULL) +
	scale_fill_manual(  values = COLS_ENV, name = NULL) +
	scale_x_continuous(breaks = c(0, 50, 100, 150),
										 expand = expansion(mult = c(0.02, 0.02))) +
	scale_y_continuous(labels = scales::label_number(accuracy = 0.01)) +
	labs(
		x = LAB_STANDAGE,
		y = expression(Delta * "VEcv (upper \u2212 lower environmental quantile, mean across traits)")
	) +
	theme_succession(base_size = 10) +
	theme(
		legend.position  = "bottom",
		legend.key.size  = unit(4, "mm"),
		strip.text       = element_text(face = "bold", size = 10),
		panel.grid.major = element_line(colour = "grey92", linewidth = 0.3)
	)

save_fig(fig4, "fig4_vecv.png", width = 180, height = 120)

# ══════════════════════════════════════════════════════════════════════════════
# Summary
# ══════════════════════════════════════════════════════════════════════════════

message("\n── main_plots.R complete ───────────────────────────────────────────")
message("Figures saved to figures/main/:")
list.files(DIR_MAIN) %>% walk(~ message(sprintf("  %s", .x)))