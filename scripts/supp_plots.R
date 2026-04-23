################################################################################
## succession_traits: supp_plots.R
## Builds all supplementary figures from pre-computed analysis outputs.
## Reads .rds files from tables/ and writes figures to figures/supplementary/.
##
## Supplementary figure order (follows main text cross-reference sequence):
##   S-1: PCA scree + varimax contributions       (methods: environmental PCA)
##   S-2: Environmental strata spatial + density  (methods: stratification)
##   S-3: Quantile sensitivity analysis           (methods: sensitivity)
##   S-4: SHAP vs linear correlation              (methods: ML justification)
##   S-5: SHAP beeswarm plots                     (RQ1 results)
##   S-6: SHAP dependence plots                   (RQ1 results)
##   S-7: PDP fit lines                           (RQ2 results)
##   S-8: Full VEcv trajectories all 9 traits     (RQ3 results)
##   S-9: Full ΔVEcv per trait                    (RQ3 results)
##
## Run after all pipeline scripts (01–06) and quantile_sensitivity.R.
##
## Author: M. Weiss @ Maynard Lab UCL / ETH Zürich
################################################################################

rm(list = ls())

# ── Libraries ─────────────────────────────────────────────────────────────────
library(tidyverse)
library(patchwork)
library(ggbeeswarm)
library(ggh4x)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(scales)
library(ggrepel)

source("scripts/plot_theme.R")
source("scripts/functions.R")

# ── Output directories ────────────────────────────────────────────────────────
DIR_SUPP     <- "figures/supplementary"
DIR_BEESWARM <- file.path(DIR_SUPP, "shap_beeswarm")
DIR_DEPEND   <- file.path(DIR_SUPP, "shap_dependence")
DIR_PDP      <- file.path(DIR_SUPP, "pdp")
DIR_VECV     <- file.path(DIR_SUPP, "vecv")
DIR_PCA      <- file.path(DIR_SUPP, "pca")

walk(c(DIR_SUPP, DIR_BEESWARM, DIR_DEPEND, DIR_PDP, DIR_VECV, DIR_PCA),
		 ~ dir.create(.x, recursive = TRUE, showWarnings = FALSE))

# ── Load tables ───────────────────────────────────────────────────────────────
shap_values      <- read_rds("tables/shap_values.rds")
shap_per_var     <- read_rds("tables/shap_per_var.rds")
pdp_raw          <- read_rds("tables/pdp_raw.rds")
vecv_raw         <- read_rds("tables/vecv_raw.rds")
vecv_divergence  <- read_rds("tables/vecv_divergence.rds")
data_clean       <- read_rds("data_processed/fia_traits_clean.rds")

# Load sensitivity outputs if available
sens_available <- file.exists("tables/sensitivity_pdp.rds") &&
	file.exists("tables/sensitivity_vecv.rds")
if (sens_available) {
	sensitivity_pdp  <- read_rds("tables/sensitivity_pdp.rds")
	sensitivity_vecv <- read_rds("tables/sensitivity_vecv.rds")
}

# ── Analytical constants ──────────────────────────────────────────────────────
ENV_VARS   <- c("temp_pc", "soil_pc", "rain_pc", "elevation", "soil_ph")
LEAF_TYPES <- c("broadleaf", "coniferous")
TRAITS     <- names(TRAIT_LABELS)
COVARIATES <- c("standage", "temp_pc", "soil_pc", "rain_pc", "elevation", "soil_ph")

THRESHOLDS <- list(
	"10/90" = c(0.10, 0.90),
	"15/85" = c(0.15, 0.85),
	"20/80" = c(0.20, 0.80),
	"25/75" = c(0.25, 0.75),
	"30/70" = c(0.30, 0.70)
)

COLS_STRATA <- c(
	"Lower 25%"  = "#4575B4",
	"Middle 50%" = "grey85",
	"Upper 25%"  = "#D73027"
)

strata_label <- function(x, q) {
	case_when(
		x <= q[1] ~ "Lower 25%",
		x >= q[2] ~ "Upper 25%",
		TRUE       ~ "Middle 50%"
	)
}

save_supp <- function(plot, filename, dir = DIR_SUPP,
											width = 180, height = 140, dpi = 300) {
	save_fig(plot, filename, dir = dir,
					 width = width, height = height, dpi = dpi)
}


# ══════════════════════════════════════════════════════════════════════════════
# S-1 — PCA scree + varimax contributions
# ══════════════════════════════════════════════════════════════════════════════
message("S-1: PCA figure — already in figures/supplementary/pca/.")
message("     To rebuild, re-run 02_environment_pca.R.")


# ══════════════════════════════════════════════════════════════════════════════
# S-2 — Environmental strata: spatial maps + feature distributions
# ══════════════════════════════════════════════════════════════════════════════
# For each of the five environmental predictors: left panel shows the spatial
# distribution of lower/middle/upper 25% strata across CONUS hexagonal plots;
# right panel shows kernel density estimates for lower and upper strata split
# by forest type. Provides geographic and numerical context for the
# stratification used throughout the PDP and VEcv analyses.

message("Building S-2 (environmental strata maps + distributions)...")

crs_conus <- 5070L
usa_wgs <- ne_countries(
	scale = "medium", country = "united states of america", returnclass = "sf"
)
conus_bbox_wgs <- st_as_sfc(
	st_bbox(c(xmin = -125, xmax = -66.5, ymin = 24.5, ymax = 49.5),
					crs = st_crs(4326))
)
usa_conus <- st_transform(st_intersection(usa_wgs, conus_bbox_wgs), crs_conus)

data_lt <- data_clean %>%
	mutate(
		leaf_type = case_when(
			biome_boreal_forests_or_taiga    == 1 |
				biome_temperate_conifer_forests  == 1 ~ "Coniferous",
			biome_temperate_broadleaf_forests == 1 |
				biome_mediterranean_woodlands    == 1 ~ "Broadleaf",
			TRUE ~ NA_character_
		)
	) %>%
	filter(!is.na(leaf_type))

pts_wgs <- st_as_sf(
	data_lt %>% dplyr::select(PID_rep, LAT, LON, all_of(ENV_VARS)),
	coords = c("LON", "LAT"), crs = 4326L
)
pts_conus <- pts_wgs %>%
	st_transform(crs_conus) %>%
	st_intersection(st_transform(conus_bbox_wgs, crs_conus))

cellsize_m <- units::drop_units(
	units::set_units(units::set_units(100, "km"), "m")
)
hex <- st_make_grid(usa_conus, cellsize = cellsize_m,
										square = FALSE, what = "polygons") %>%
	st_sf(geometry = .) %>%
	st_make_valid() %>%
	st_intersection(st_make_valid(st_union(usa_conus))) %>%
	st_collection_extract("POLYGON") %>%
	filter(st_is_valid(.)) %>%
	mutate(.area = as.numeric(st_area(.))) %>%
	filter(.area > 0) %>%
	dplyr::select(-.area) %>%
	mutate(id = row_number())

j <- st_join(pts_conus, hex, join = st_within, left = FALSE)

hex_agg <- j %>%
	st_drop_geometry() %>%
	group_by(id) %>%
	summarise(across(all_of(ENV_VARS), ~ median(.x, na.rm = TRUE)),
						n = n(), .groups = "drop")

hex_vals <- hex %>%
	left_join(hex_agg, by = "id") %>%
	filter(!is.na(n) & n >= 3)

env_qs <- map(ENV_VARS, function(env) {
	quantile(data_lt[[env]], c(0.25, 0.75), na.rm = TRUE)
}) %>% set_names(ENV_VARS)

hex_strata <- hex_vals %>%
	mutate(across(all_of(ENV_VARS),
								~ strata_label(.x, env_qs[[cur_column()]]),
								.names = "{.col}_strata"))

theme_map_s2 <- function() {
	theme_void(base_size = 8) +
		theme(
			legend.position      = "top",
			legend.title         = element_blank(),
			legend.text          = element_text(size = 6),
			legend.key.height    = unit(2, "mm"),
			legend.key.width     = unit(5, "mm"),
			legend.margin        = margin(0, 0, 0, 0),
			legend.box.spacing   = unit(0, "mm"),
			plot.margin          = margin(0, 0, 0, 0),
			plot.subtitle        = element_text(size = 7, margin = margin(0,0,1,0))
		)
}

build_s2_row <- function(env_var, env_label) {
	
	strata_col <- paste0(env_var, "_strata")
	qs         <- env_qs[[env_var]]
	
	p_map <- ggplot() +
		geom_sf(data = hex_strata,
						aes(fill = .data[[strata_col]]), colour = NA) +
		geom_sf(data = usa_conus, fill = NA, colour = "black",
						linewidth = 0.15) +
		scale_fill_manual(
			values   = COLS_STRATA,
			na.value = "grey90",
			name     = NULL,
			breaks   = c("Lower 25%", "Middle 50%", "Upper 25%"),
			guide    = guide_legend(direction    = "horizontal",
															override.aes = list(colour = NA))
		) +
		coord_sf(expand = FALSE) +          # remove padding around map extent
		theme_map_s2() +
		labs(subtitle = env_label) +
		theme(plot.margin = margin(2, 0, 0, 0))   # zero margin on map panel
	
	dens_data <- data_lt %>%
		dplyr::select(leaf_type, val = all_of(env_var)) %>%
		mutate(stratum = strata_label(val, qs)) %>%
		filter(stratum != "Middle 50%") %>%
		mutate(stratum = factor(stratum, levels = c("Lower 25%", "Upper 25%")))
	
	p_dens <- ggplot(dens_data,
									 aes(x = val, fill = stratum, colour = stratum,
									 		linetype = leaf_type)) +
		geom_density(alpha = 0.35, linewidth = 0.5) +
		scale_fill_manual(
			values = COLS_STRATA[c("Lower 25%", "Upper 25%")],
			name   = "Stratum"
		) +
		scale_colour_manual(
			values = COLS_STRATA[c("Lower 25%", "Upper 25%")],
			name   = "Stratum"
		) +
		scale_linetype_manual(
			values = c("Broadleaf" = "solid", "Coniferous" = "dashed"),
			name   = "Forest type"
		) +
		labs(x = env_label, y = "Density") +
		theme_succession(base_size = 8) +
		theme(
			legend.position       = "right",
			legend.key.size       = unit(3, "mm"),
			legend.margin         = margin(0, 0, 0, 4),
			legend.box.spacing    = unit(2, "mm"),
			panel.grid.minor      = element_blank(),
			plot.margin           = margin(0, 2, 0, 2)
		)
	
	(p_map | p_dens) +
		plot_layout(widths = c(0.65, 0.35))
}

s2_rows <- map2(ENV_VARS, ENV_LABELS, build_s2_row)

s2 <- wrap_plots(s2_rows, ncol = 1) +
	theme(plot.margin = margin(0, 1, 0, 1))

s2 <- s2 +
	plot_annotation(tag_levels = "a",
									theme = theme(plot.margin = margin(2, 2, 2, 2)))

save_supp(s2, "S2_strata_spatial_distributions.png",
					width = 220, height = 300, dpi = 300)
message("  ✓ S-2 saved")


# ══════════════════════════════════════════════════════════════════════════════
# S-3 — Quantile sensitivity analysis
# ══════════════════════════════════════════════════════════════════════════════
message("S-3: Quantile sensitivity — checking for existing file...")
sens_src <- "figures/supplementary/S3_quantile_sensitivity.png"
sens_dst <- file.path(DIR_SUPP, "S3_quantile_sensitivity.png")
if (file.exists(sens_src)) {
	message("  ✓ S-3 already present")
} else {
	message("  ⚠ S-3 not found — run quantile_sensitivity.R")
}


# ══════════════════════════════════════════════════════════════════════════════
# S-4 — SHAP importance vs linear trait-feature correlation
# ══════════════════════════════════════════════════════════════════════════════

message("Building S-4 (SHAP vs correlation)...")

splits_bl <- read_rds("models/splits_broadleaf.rds")
splits_co <- read_rds("models/splits_coniferous.rds")

test_data <- bind_rows(
	splits_bl$test %>% mutate(leaf_type = "broadleaf"),
	splits_co$test %>% mutate(leaf_type = "coniferous")
)

cor_data <- map_dfr(LEAF_TYPES, function(lt) {
	df_lt <- test_data %>% filter(leaf_type == lt)
	map_dfr(TRAITS, function(tr) {
		map_dfr(COVARIATES, function(cov) {
			r <- cor(df_lt[[cov]], df_lt[[tr]], use = "pairwise.complete.obs")
			tibble(trait = tr, variable = cov, leaf_type = lt, correlation = r)
		})
	})
}) %>%
	mutate(variable_label = recode(variable, !!!ALL_PREDICTOR_LABELS))

s4_data <- shap_per_var %>%
	dplyr::select(trait, leaf_type, variable, mean_abs_shap, variable_label) %>%
	left_join(
		cor_data %>% dplyr::select(-variable_label),
		by = c("trait", "leaf_type", "variable")
	) %>%
	filter(!is.na(correlation)) %>%
	mutate(
		leaf_type   = str_to_title(leaf_type),
		trait_label = recode(trait, !!!TRAIT_LABELS)
	)

s4 <- ggplot(s4_data,
						 aes(x = correlation, y = mean_abs_shap,
						 		colour = variable_label, shape = leaf_type)) +
	geom_hline(yintercept = 0, linetype = "dashed",
						 colour = "grey70", linewidth = 0.3) +
	geom_vline(xintercept = 0, linetype = "dashed",
						 colour = "grey70", linewidth = 0.3) +
	geom_point(size = 2.5, alpha = 0.85, stroke = 0.4) +
	ggrepel::geom_text_repel(
		data = s4_data %>%
			group_by(leaf_type) %>%
			slice_max(mean_abs_shap, n = 6),
		aes(label = paste0(trait_label, "\n(", variable_label, ")")),
		size = 2.2, show.legend = FALSE,
		max.overlaps = 12, segment.colour = "grey60"
	) +
	scale_colour_manual(values = COLS_PREDICTORS, name = "Predictor") +
	scale_shape_manual(values  = SHAPES_LEAFTYPE,  name = "Forest type") +
	facet_wrap(~ leaf_type, ncol = 2) +
	labs(
		x = "Pearson correlation (feature value vs trait CWM)",
		y = "Mean |SHAP| importance"
	) +
	theme_succession(base_size = 10) +
	theme(
		legend.position = "bottom"
	) +
	guides(shape = "none")

save_supp(s4, "S4_shap_vs_correlation.png",
					width = 180, height = 120)
message("  ✓ S-4 saved")

# ══════════════════════════════════════════════════════════════════════════════
# S-5 — SHAP beeswarm plots (one per trait × leaf type)
# ══════════════════════════════════════════════════════════════════════════════

message("Building S-5 (SHAP beeswarm compound figures)...")

MAX_POINTS <- 2000L  # reduce per-panel points since 9 panels per page

for (lt in LEAF_TYPES) {
	
	plot_list <- map(TRAITS, function(tr) {
		
		var_order <- shap_values %>%
			filter(trait == tr, leaf_type == lt) %>%
			group_by(variable) %>%
			summarise(mean_abs = mean(abs(shap_value), na.rm = TRUE),
								.groups = "drop") %>%
			arrange(mean_abs) %>%
			pull(variable)
		
		df_sub <- shap_values %>%
			filter(trait == tr, leaf_type == lt) %>%
			group_by(variable) %>%
			mutate(value_col = symmetric_scale(feature_value)) %>%
			ungroup() %>%
			mutate(variable = factor(variable, levels = var_order)) %>%
			{ if (nrow(.) > MAX_POINTS) slice_sample(., n = MAX_POINTS) else . }
		
		ggplot(df_sub,
					 aes(x = variable, y = shap_value, colour = value_col)) +
			ggbeeswarm::geom_quasirandom(alpha = 0.35, size = 0.6) +
			geom_hline(yintercept = 0, linetype = "dashed",
								 colour = "grey50", linewidth = 0.3) +
			coord_flip() +
			scale_x_discrete(labels = ALL_PREDICTOR_LABELS) +
			scale_colour_viridis_c(
				name   = "Feature value",
				limits = c(-1, 1),
				breaks = c(-1, 0, 1),
				labels = c("Low", "Mid", "High"),
				guide  = guide_colourbar(barheight = unit(2, "mm"),
																 barwidth  = unit(15, "mm"))
			) +
			labs(
				x        = NULL,
				y        = "SHAP value",
				subtitle = TRAIT_LABELS[[tr]]
			) +
			theme_succession(base_size = 7) +
			theme(
				panel.grid.major.y = element_blank(),
				axis.text.y        = element_text(size = 6),
				legend.position    = "none",
				plot.subtitle      = element_text(face = "bold", size = 7,
																					hjust = 0.5)
			)
	})
	
	# Extract legend from one panel to share across the compound figure
	legend_panel <- ggplot(
		shap_values %>%
			filter(trait == TRAITS[1], leaf_type == lt) %>%
			mutate(value_col = symmetric_scale(feature_value)) %>%
			slice(1:100),
		aes(x = variable, y = shap_value, colour = value_col)
	) +
		ggbeeswarm::geom_quasirandom() +
		scale_colour_viridis_c(
			name   = "Feature value",
			limits = c(-1, 1),
			breaks = c(-1, 0, 1),
			labels = c("Low", "Mid", "High")
		) +
		theme_succession() +
		theme(legend.position = "bottom")
	
	shared_legend <- cowplot::get_legend(legend_panel)
	
	# Assemble 3×3 grid + shared legend at bottom
	compound <- wrap_plots(plot_list, ncol = 3, nrow = 3) +
		plot_annotation(
			title    = paste0(str_to_title(lt), " forests"),
			tag_levels = "a",
			theme    = theme(
				plot.title  = element_text(face = "bold", size = 9, hjust = 0.5),
				plot.margin = margin(4, 4, 4, 4)
			)
		)
	
	# Add shared legend below using cowplot
	final <- cowplot::plot_grid(
		compound,
		shared_legend,
		ncol    = 1,
		rel_heights = c(1, 0.05)
	)
	
	save_supp(final,
						sprintf("S5_beeswarm_%s.png", lt),
						dir    = DIR_BEESWARM,
						width  = 220,
						height = 260,
						dpi    = 300)
	message(sprintf("  ✓ S-5 beeswarm %s saved", lt))
}


# ══════════════════════════════════════════════════════════════════════════════
# S-6 — SHAP dependence plots (one per variable × trait)
# ══════════════════════════════════════════════════════════════════════════════

message("Building S-6 (SHAP dependence compound figures)...")

for (tr in TRAITS) {
	
	plot_list <- map(COVARIATES, function(var) {
		
		shap_values %>%
			filter(trait == tr, variable == var) %>%
			mutate(leaf_type = str_to_title(leaf_type)) %>%
			ggplot(aes(x = feature_value, y = shap_value,
								 colour = leaf_type, shape = leaf_type)) +
			geom_point(alpha = 0.15, size = 0.7) +
			geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
									se = FALSE, linewidth = 0.8) +
			geom_hline(yintercept = 0, linetype = "dashed",
								 colour = "grey50", linewidth = 0.3) +
			scale_colour_manual(values = COLS_LEAFTYPE, name = NULL) +
			scale_shape_manual(values  = SHAPES_LEAFTYPE, name = NULL) +
			labs(
				x        = ALL_PREDICTOR_LABELS[[var]],
				y        = "SHAP value",
				subtitle = ALL_PREDICTOR_LABELS[[var]]
			) +
			theme_succession(base_size = 7) +
			theme(
				legend.position  = "none",
				panel.grid.minor = element_blank(),
				plot.subtitle    = element_text(face = "bold", size = 7,
																				hjust = 0.5)
			)
	})
	
	# Shared legend
	legend_panel <- shap_values %>%
		filter(trait == tr, variable == COVARIATES[1]) %>%
		mutate(leaf_type = str_to_title(leaf_type)) %>%
		ggplot(aes(x = feature_value, y = shap_value,
							 colour = leaf_type, shape = leaf_type)) +
		geom_point() +
		scale_colour_manual(values = COLS_LEAFTYPE, name = "Forest type") +
		scale_shape_manual(values  = SHAPES_LEAFTYPE, name = "Forest type") +
		theme_succession() +
		theme(legend.position = "bottom")
	
	shared_legend <- cowplot::get_legend(legend_panel)
	
	# 2×3 grid (6 predictors: standage + 5 env vars)
	compound <- wrap_plots(plot_list, ncol = 3, nrow = 2) +
		plot_annotation(
			title      = TRAIT_LABELS[[tr]],
			tag_levels = "a"
		)
	
	final <- cowplot::plot_grid(
		compound,
		shared_legend,
		ncol        = 1,
		rel_heights = c(1, 0.06)
	)
	
	save_supp(final,
						sprintf("S6_dependence_%s.png", tr),
						dir    = DIR_DEPEND,
						width  = 200,
						height = 160,
						dpi    = 300)
	message(sprintf("  ✓ S-6 dependence %s saved", tr))
}


# ══════════════════════════════════════════════════════════════════════════════
# S-7 — PDP fit lines (raw bootstrap curves)
# ══════════════════════════════════════════════════════════════════════════════

message("Building S-7 (PDP fit lines)...")

s7 <- pdp_raw %>%
	mutate(
		leaf_type      = str_to_title(leaf_type),
		trait_label    = recode(trait,    !!!TRAIT_LABELS),
		variable_label = recode(variable, !!!ENV_LABELS),
		trait_label    = factor(trait_label,    levels = TRAIT_LABELS),
		variable_label = factor(variable_label, levels = ENV_LABELS)
	) %>%
	ggplot(aes(x = standage, y = yhat)) +
	geom_point(aes(fill = group), alpha = 0.03, size = 0.3,
						 shape = 21, stroke = 0.1, colour = "black") +
	geom_smooth(aes(colour = group), method = "lm",
							se = FALSE, linewidth = 0.5, formula = y ~ x) +
	ggh4x::facet_nested(
		rows = vars(leaf_type, trait_label),
		cols = vars(variable_label),
		scales = "free_y"
	) +
	scale_fill_manual(
		values = c("low" = "lightcyan", "high" = "lightcoral"),
		labels = c("low" = "Lower quantile", "high" = "Upper quantile"),
		name   = NULL
	) +
	scale_colour_manual(
		values = c("low" = "navy", "high" = "darkred"),
		labels = c("low" = "Lower quantile", "high" = "Upper quantile"),
		name   = NULL
	) +
	labs(x = LAB_STANDAGE, y = "Predicted trait expression") +
	theme_succession(base_size = 7) +
	theme(legend.position = "top")

save_supp(s7, "S7_pdp_lines.png",
					dir = DIR_PDP, width = 260, height = 360)
message("  ✓ S-7 saved")


# ══════════════════════════════════════════════════════════════════════════════
# S-8 — Full VEcv trajectories (all 9 traits)
# ══════════════════════════════════════════════════════════════════════════════

message("Building S-8 (full VEcv trajectories)...")

vecv_avg_full <- vecv_raw %>%
	group_by(trait, trait_label, leaf_type, env_group,
					 standage_bin, standage_mid, repeat_id) %>%
	summarise(VEcv = mean(VEcv, na.rm = TRUE), .groups = "drop") %>%
	group_by(trait, trait_label, leaf_type, env_group,
					 standage_bin, standage_mid) %>%
	summarise(
		VEcv_med = median(VEcv, na.rm = TRUE),
		VEcv_lwr = quantile(VEcv, 0.025, na.rm = TRUE),
		VEcv_upr = quantile(VEcv, 0.975, na.rm = TRUE),
		.groups  = "drop"
	) %>%
	mutate(
		env_group   = factor(env_group,
												 levels = c("low", "high"),
												 labels = c("Lower environmental quantile",
												 					 "Upper environmental quantile")),
		leaf_type   = str_to_title(leaf_type),
		trait_label = factor(trait_label, levels = TRAIT_LABELS)
	)

s8 <- ggplot(vecv_avg_full,
						 aes(x = standage_mid, colour = env_group, fill = env_group)) +
	geom_ribbon(aes(ymin = VEcv_lwr, ymax = VEcv_upr),
							alpha = 0.2, colour = NA) +
	geom_line(aes(y = VEcv_med), linewidth = 0.7) +
	geom_hline(yintercept = 0, linetype = "dashed",
						 colour = "grey50", linewidth = 0.3) +
	facet_grid(leaf_type ~ trait_label, scales = "free_y") +
	scale_colour_manual(values = COLS_ENVGROUP, name = NULL) +
	scale_fill_manual(  values = COLS_ENVGROUP, name = NULL) +
	scale_x_continuous(breaks = c(0, 50, 100, 150)) +
	scale_y_continuous(limits = c(0, NA)) +
	labs(
		x        = LAB_STANDAGE,
		y        = LAB_VECV,
	) +
	theme_succession(base_size = 8) +
	theme(
		legend.position = "top",
		strip.text      = element_text(face = "bold", size = 7)
	)

save_supp(s8, "S8_vecv_full.png",
					dir = DIR_VECV, width = 260, height = 140)
message("  ✓ S-8 saved")


# ══════════════════════════════════════════════════════════════════════════════
# S-9 — Full ΔVEcv per trait × environmental variable
# ══════════════════════════════════════════════════════════════════════════════

message("Building S-9 (full ΔVEcv per trait)...")

s9 <- vecv_divergence %>%
	filter(!is.na(delta_med)) %>%
	mutate(
		leaf_type      = str_to_title(leaf_type),
		trait_label    = factor(trait_label,    levels = TRAIT_LABELS),
		variable_label = factor(variable_label, levels = ENV_LABELS)
	) %>%
	ggplot(aes(x = standage_mid, y = delta_med,
						 colour = variable_label, fill = variable_label)) +
	geom_ribbon(aes(ymin = delta_lwr, ymax = delta_upr),
							alpha = 0.15, colour = NA) +
	geom_line(linewidth = 0.6) +
	geom_hline(yintercept = 0, linetype = "dashed",
						 colour = "grey40", linewidth = 0.35) +
	facet_grid(leaf_type ~ trait_label, scales = "free_y") +
	scale_colour_manual(values = COLS_ENV, name = NULL) +
	scale_fill_manual(  values = COLS_ENV, name = NULL) +
	scale_x_continuous(breaks = c(0, 50, 100, 150)) +
	labs(x = LAB_STANDAGE, y = LAB_DELTA_VECV) +
	theme_succession(base_size = 8) +
	theme(
		legend.position = "bottom",
		strip.text      = element_text(face = "bold", size = 7)
	)

save_supp(s9, "S9_vecv_divergence_full.png",
					dir = DIR_VECV, width = 260, height = 140)
message("  ✓ S-9 saved")


# ══════════════════════════════════════════════════════════════════════════════
# Summary
# ══════════════════════════════════════════════════════════════════════════════

message("\n── supp_plots.R complete ───────────────────────────────────────────")
message("Supplementary figures saved to figures/supplementary/:")
list.files(DIR_SUPP, recursive = FALSE, pattern = "\\.png$") %>%
	sort() %>%
	walk(~ message(sprintf("  %s", .x)))