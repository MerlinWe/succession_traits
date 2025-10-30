# not my best work 

library(ggplot2)
library(ggbump)
library(dplyr)
library(viridis)
library(MetBrewer)
library(showtext)
library(scales)
library(patchwork)

importance_rank <- importance_rank %>%
	mutate(
		trait = factor(
			trait,
			levels = c(
				"height",
				"bark_thickness",
				"conduit_diam",
				"root_depth",
				"specific_leaf_area",
				"leaf_density",
				"leaf_k",
				"seed_dry_mass",
				"shade_tolerance"
			)
		)
	)

trait_labels <- c(
	"height"             = "Height",
	"bark_thickness"     = "Bark thickness",
	"conduit_diam"       = "Conduit diameter",
	"root_depth"         = "Root depth",
	"specific_leaf_area" = "Specific leaf area (SLA)",
	"leaf_density"       = "Leaf density",
	"leaf_k"             = "Leaf potassium (K)",
	"seed_dry_mass"      = "Seed dry mass",
	"shade_tolerance"    = "Shade tolerance"
)

feature_labels <- c(
	"standage"  = "Stand age",
	"temp_pc"   = "Temperature",
	"rain_pc"   = "Precipitation",
	"soil_pc"   = "Soil W.R.",
	"soil_ph"   = "Soil pH",
	"elevation" = "Elevation"
)

feature_colors <- c(
	"standage"  = "#8B5A2B",  # brown-green
	"temp_pc"   = "#E6550D",  # warm orange
	"rain_pc"   = "#3182BD",  # blue
	"soil_pc"   = "#7B9C4E",  # olive green
	"soil_ph"   = "#756BB1",  # purple
	"elevation" = "#969696"   # gray
)


highlight_vars <- c("standage", "temp_pc")

p_bump <- ggplot(
	importance_rank,
	aes(
		x = trait,
		y = rank,
		group = interaction(leaf_type, variable),
		color = variable
	)
) +
	geom_bump(linewidth = 0.8, alpha = 0.3, smooth = 8) +
	geom_bump(
		data = ~ .x %>% filter(variable %in% highlight_vars),
		linewidth = 1.4, alpha = 0.95, smooth = 8
	) +
	geom_point(size = 1.6, alpha = 0.8) +
	geom_text(
		aes(label = recode(variable, !!!feature_labels)),
		x = 9.6, hjust = 0,
		color = "gray30",
		family = "Roboto Condensed", size = 3.2,
		data = importance_rank %>%
			group_by(leaf_type, variable) %>%
			slice_max(order_by = as.numeric(trait)) %>%
			ungroup()
	) +
	scale_color_manual(values = feature_colors) +
	scale_y_reverse(breaks = 1:6, expand = c(0.01, 0.01)) +
	facet_wrap(
		~ leaf_type, ncol = 1, scales = "free_x",
		labeller = as_labeller(c("broadleaf" = "Broadleaf", "coniferous" = "Coniferous"))
	) +
	scale_x_discrete(labels = trait_labels) +
	coord_cartesian(clip = "off") +
	theme_minimal(base_family = "Roboto Condensed", base_size = 12) +
	labs(
		y = "Relative Rank Importance\n(1 = highest mean |SHAP| importance)",
		x = NULL
	) +
	theme(
		legend.position = "none",
		panel.grid = element_blank(),
		plot.title.position = "plot",
		plot.margin = margin(5, 10, 5, 5),
		axis.text.x = element_text(angle = 35, hjust = 1, size = 10),
		axis.text.y = element_text(size = 11, margin = margin(r = 2)),
		axis.title.y = element_text(size = 11, margin = margin(r = 2)),
		strip.text = element_text(face = "bold", size = 11, family = "Roboto Condensed"),
		plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
		plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40")
	)

# ---------- panel B ----------

shap_slopes <- shap_by_leaf %>%
	group_by(leaf_type, trait, variable) %>%
	summarise(
		slope     = tryCatch(coef(lm(shap_value ~ feature_value))[2], error = function(e) NA),
		mean_shap = mean(shap_value, na.rm = TRUE),
		mean_abs  = mean(abs(shap_value), na.rm = TRUE),
		.groups = "drop"
	) %>%
	filter(!is.na(slope))

shap_slopes <- shap_slopes %>%
	mutate(
		arrow = case_when(
			slope >  0.06  ~ "↗",   # strong positive
			slope >  0.015 ~ "↑",   # weak positive
			slope < -0.06  ~ "↙",   # strong negative
			slope < -0.015 ~ "↓",   # weak negative
			TRUE           ~ "•"    # near-zero slope
		),
		alpha = rescale(mean_abs, to = c(0.4, 1))
	)

var_labels <- c(
	"standage"   = "Stand age",
	"temp_pc"    = "Temperature",
	"rain_pc"    = "Precipitation",
	"soil_pc"    = "Soil W.R.",
	"elevation"  = "Elevation",
	"soil_ph"    = "Soil pH"
)

trait_labels <- c(
	"bark_thickness"     = "Bark thickness",
	"conduit_diam"       = "Conduit diameter",
	"height"             = "Height",
	"leaf_density"       = "Leaf density",
	"leaf_k"             = "Leaf potassium",
	"root_depth"         = "Root depth",
	"seed_dry_mass"      = "Seed dry mass",
	"shade_tolerance"    = "Shade tolerance",
	"specific_leaf_area" = "Specific leaf area"
)

shap_slopes <- shap_slopes %>%
	mutate(
		variable_num = match(variable, names(var_labels)),
		trait_num = match(trait, names(trait_labels))
	)

p_arrow <- ggplot(shap_slopes, aes(x = variable_num, y = trait_num)) +
	geom_tile(aes(fill = mean_shap), color = "white", linewidth = 0.3) +
	geom_text(aes(label = arrow, alpha = alpha),
						size = ifelse(shap_slopes$arrow == "•", 3.5, 5)) +
	scale_fill_gradient2(
		low = "#4575B4", mid = "white", high = "#D73027",
		midpoint = 0, name = "Mean SHAP\n(sign)"
	) +
	scale_alpha_identity() +
	scale_x_continuous(
		breaks = 1:length(var_labels),
		labels = var_labels,
		expand = expansion(mult = c(0.05, 0.05))
	) +
	scale_y_continuous(
		breaks = 1:length(trait_labels),
		labels = trait_labels,
		expand = expansion(mult = c(0.05, 0.05))
	) +
	facet_wrap(
		~ leaf_type,
		labeller = as_labeller(c("broadleaf" = "Broadleaf", "coniferous" = "Coniferous"))
	) +
	coord_fixed() +
	labs(x = NULL, y = NULL) +
	theme_bw(base_family = "Roboto Condensed", base_size = 12) +
	theme(
		panel.grid = element_blank(),
		axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
		axis.text.y = element_text(size = 11, face = "bold"),
		axis.title.y = element_text(size = 11, margin = margin(r = 2)),
		strip.text = element_text(face = "bold", size = 11, family = "Roboto Condensed"),
		panel.border = element_rect(color = "grey60"),
		legend.position = "right",
		legend.key.height = unit(22, "pt"),
		plot.margin = margin(5, 10, 5, 5)
	)


p_final <- p_bump / p_arrow +
	plot_layout(heights = c(1.2, 0.8)) + 
	theme(plot.margin = margin(0, 0, 0, 0)) +
	plot_annotation(
		tag_levels = "a",
		theme = theme(
			plot.tag = element_text(
				face = "bold",
				size = 14,
				family = "Roboto Condensed",
				margin = margin(b = 1)
			)
		)
	)

ggsave(
	filename = "plots/fig2.eps",
	plot = p_final,
	device = cairo_ps,
	width = 200, height = 240, units = "mm",
	fallback_resolution = 600,
	limitsize = FALSE,
	bg = "transparent"
)


# Trait‐specific importance and directionality of environmental predictors for broadleaf and coniferous tree species.
# (a) Relative rank importance of six environmental predictors across nine functional traits, based on mean absolute SHAP values from random forest models. Lower ranks indicate stronger average predictive importance (1 = highest). Lines connect ranks across traits, separately for broadleaf (top) and coniferous (bottom) species, with colours denoting individual predictors. The ordering of traits along the x‐axis follows the analytical sequence of trait models.
# (b) Direction and strength of SHAP associations between predictors (x‐axis) and traits (y‐axis) for each leaf type. Tile colour shows the mean signed SHAP value (red = positive effect, blue = negative effect, white = weak or inconsistent effect). Arrow symbols indicate the dominant direction of the SHAP–feature relationship: upward (↑) for positive, downward (↓) for negative, and dot (•) for weak or inconsistent responses.




# --- Custom feature colours ---
feature_colors <- c(
	"standage"  = "#8B5A2B",  # brown-green
	"temp_pc"   = "#E6550D",  # warm orange
	"rain_pc"   = "#3182BD",  # blue
	"soil_pc"   = "#7B9C4E",  # olive green
	"soil_ph"   = "#756BB1",  # purple
	"elevation" = "#969696"   # gray
)

# --- Clean variable labels ---
var_labels <- c(
	"standage"   = "Stand age",
	"temp_pc"    = "Temperature",
	"rain_pc"    = "Precipitation",
	"soil_pc"    = "Soil water retention",
	"elevation"  = "Elevation",
	"soil_ph"    = "Soil pH"
)

trait_labels <- c(
	"bark_thickness"     = "Bark thickness",
	"conduit_diam"       = "Conduit diameter",
	"height"             = "Height",
	"leaf_density"       = "Leaf density",
	"leaf_k"             = "Leaf potassium",
	"root_depth"         = "Root depth",
	"seed_dry_mass"      = "Seed mass",
	"shade_tolerance"    = "Shade tolerance",
	"specific_leaf_area" = "SLA"
)

# --- Compute normalised SHAP importance ---
shap_importance <- shap_by_leaf %>%
	group_by(trait, leaf_type, variable) %>%
	summarise(mean_abs = mean(abs(shap_value), na.rm = TRUE), .groups = "drop") %>%
	group_by(trait, leaf_type) %>%
	mutate(prop_importance = mean_abs / sum(mean_abs) * 100) %>%
	ungroup()

# --- Build plot ---
p_bar <- ggplot(
	shap_importance,
	aes(
		x = factor(trait, levels = names(trait_labels)),
		y = prop_importance,
		fill = variable
	)
) +
	geom_col(position = "stack", width = 0.75, color = "white", linewidth = 0.2) +
	facet_wrap(
		~leaf_type, ncol = 1,
		labeller = as_labeller(c("broadleaf"="Broadleaf", "coniferous"="Coniferous"))
	) +
	scale_x_discrete(labels = trait_labels) +
	scale_fill_manual(
		name = "Predictor",
		values = feature_colors,
		labels = var_labels
	) +
	labs(
		y = "Relative importance (% of total |SHAP|)",
		x = NULL
	) +
	theme_bw(base_size = 11) +
	theme(
		panel.grid.minor = element_blank(),
		panel.grid.major.x = element_blank(),
		axis.text.x = element_text(angle = 35, hjust = 1, size = 9),
		axis.text.y = element_text(size = 9),
		axis.title.y = element_text(size = 10, face = "bold"),
		strip.text = element_text(face = "bold"),
		legend.position = "right",
		legend.key.height = unit(4, "mm"),
		legend.title = element_text(size = 9, face = "bold"),
		legend.text = element_text(size = 8)
	)


ggplot(shap_importance,
			 aes(x = variable, y = prop_importance, fill = variable)) +
	geom_col(position = position_dodge(width = 0.8), width = 0.7) +
	facet_wrap(~trait + leaf_type, ncol = 6, scales = "free_y") +
	
	scale_x_discrete(labels = trait_labels) +
	scale_fill_manual(
		name = "Predictor",
		values = feature_colors,
		labels = var_labels
	) +
	labs(
		y = "Relative importance (% of total |SHAP|)",
		x = NULL
	) +
	theme_bw(base_size = 11) +
	theme(
		panel.grid.minor = element_blank(),
		panel.grid.major.x = element_blank(),
		axis.text.x = element_text(angle = 35, hjust = 1, size = 9),
		axis.text.y = element_text(size = 9),
		axis.title.y = element_text(size = 10, face = "bold"),
		strip.text = element_text(face = "bold"),
		legend.position = "right",
		legend.key.height = unit(4, "mm"),
		legend.title = element_text(size = 9, face = "bold"),
		legend.text = element_text(size = 8)
	)
