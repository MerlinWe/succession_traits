
feature_colors <- c(
	"standage"  = "#8C6D53",  # warm oak brown — distinct from elevation
	"temp_pc"   = "#D08A61",  # soft terracotta (temperature)
	"rain_pc"   = "#6C91B2",  # misty blue (precipitation)
	"soil_pc"   = "#7D8F70",  # mossy green-gray (soil)
	"soil_ph"   = "#8C84A1",  # muted lavender-gray (pH)
	"elevation" = "#A8B6B0"   # cool pale slate (elevation)
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

var_labels <- c(
	"standage"   = "Stand age",
	"temp_pc"    = "Temperature",
	"rain_pc"    = "Precipitation",
	"soil_pc"    = "Soil water retention",
	"elevation"  = "Elevation",
	"soil_ph"    = "Soil pH"
)

# =========== Panel A ===========

feature_patterns <- c(
	"standage"  = "none",       # simple flat colour
	"temp_pc"   = "stripe",     # diagonal hatching
	"rain_pc"   = "crosshatch", # grid pattern
	"soil_pc"   = "wave",       # subtle wavy pattern
	"soil_ph"   = "circle",     # dotted texture
	"elevation" = "none"        # pale flat fill
)

shap_importance <- shap_by_leaf %>%
	group_by(trait, leaf_type, variable) %>%
	summarise(mean_abs = mean(abs(shap_value), na.rm = TRUE), .groups = "drop") %>%
	group_by(trait, leaf_type) %>%
	mutate(prop_importance = mean_abs / sum(mean_abs) * 100) %>%
	ungroup()

ggplot(
	shap_importance,
	aes(
		x = factor(trait, levels = names(trait_labels)),
		y = prop_importance,
		fill = variable,
		pattern = variable
	)
) +
	geom_col_pattern(
		position = "stack",
		width = 0.75,
		color = "white",
		linewidth = 0.2,
		pattern_fill = "grey60",    # lighter pattern colour
		pattern_alpha = 0.4,        # semi-transparent overlay
		pattern_density = 0.03,     # sparser pattern
		pattern_spacing = 0.06,     # more space between lines
		pattern_angle = 45
	) +
	facet_wrap(
		~leaf_type, ncol = 1,
		labeller = as_labeller(c("broadleaf"="Broadleaf", "coniferous"="Coniferous"))
	) +
	scale_x_discrete(labels = trait_labels) +
	scale_fill_manual(values = feature_colors, labels = var_labels, name = "Predictor") +
	scale_pattern_manual(values = feature_patterns, labels = var_labels, name = "Predictor") +
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

# =========== Panel B ===========

# Scale slope and SHAP for visual encoding
shap_slopes_fixed <- shap_by_leaf %>%
	group_by(leaf_type, trait, variable) %>%
	summarise(
		slope     = tryCatch(coef(lm(shap_value ~ feature_value))[2], error = function(e) NA),
		mean_shap = mean(shap_value, na.rm = TRUE),
		mean_abs  = mean(abs(shap_value), na.rm = TRUE),
		.groups = "drop"
	) %>%
	filter(!is.na(slope)) %>%
	mutate(
		variable_num = match(variable, names(var_labels)),
		trait_num    = match(trait, names(trait_labels)),
		# scale slope and SHAP to ellipse parameters
		size   = rescale(mean_abs, to = c(0.15, 0.35)),   # overall size inside cell
		angle  = slope * 400,                             # tilt based on slope
		ratio  = pmax(0.3, 1 - abs(slope) * 10)           # shape elongation
	)

# Plot
ggplot(shap_slopes_fixed, aes(x = variable_num, y = trait_num)) +
	geom_tile(fill = "gray97", color = "white", linewidth = 0.3) +
	
	geom_ellipse(
		aes(
			x0 = variable_num,
			y0 = trait_num,
			a  = size, 
			b  = size * ratio, 
			angle = angle,
			fill = mean_shap
		),
		color = "grey70",
		alpha = 1,
		linewidth = 0.1
	) +
	
	scale_fill_gradient2(
		low = "#4575B4", mid = "white", high = "#D73027",
		midpoint = 0, name = "Mean SHAP\n(sign)"
	) +
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
		~leaf_type,
		labeller = as_labeller(c("broadleaf" = "Broadleaf", "coniferous" = "Coniferous"))
	) +
	coord_fixed() +
	theme_bw(base_family = "Helvetica", base_size = 11) +
	labs(x = NULL, y = NULL) +
	theme(
		panel.grid = element_blank(),
		axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
		axis.text.y = element_text(size = 11, face = "bold"),
		strip.text = element_text(face = "bold", size = 11, family = "Helvetica"),
		panel.border = element_rect(color = "grey60"),
		legend.position = "right",
		legend.key.height = unit(20, "pt"),
		plot.margin = margin(5, 10, 5, 5)
	)



# ---- matrix of combinations for panel b interpretation ----


make_ellipse <- function(x0, y0, a, b, angle, n = 200) {
	t <- seq(0, 2 * pi, length.out = n)
	x <- a * cos(t)
	y <- b * sin(t)
	xr <- x * cos(angle) - y * sin(angle) + x0
	yr <- x * sin(angle) + y * cos(angle) + y0
	tibble(x = xr, y = yr)
}

legend_data <- expand.grid(
	color_cat = c("Negative (blue)", "Neutral (white)", "Positive (red)"),
	slope_cat = c("Downward (↘)", "Flat (◯)", "Upward (↗)"),
	stringsAsFactors = FALSE
) %>%
	mutate(
		mean_shap = case_when(
			color_cat == "Positive (red)" ~  0.002,
			color_cat == "Neutral (white)" ~  0,
			TRUE                           ~ -0.002
		),
		slope = case_when(
			slope_cat == "Upward (↗)"   ~  0.15,
			slope_cat == "Downward (↘)" ~ -0.15,
			TRUE                        ~  0
		),
		radius_x = 0.35,
		radius_y = 0.15,
		tilt = slope * 3
	)

ellipses <- legend_data %>%
	rowwise() %>%
	mutate(ellipse = list(make_ellipse(0, 0, radius_x, radius_y, tilt))) %>%
	unnest(ellipse)

p_legend <- ggplot(ellipses, aes(x, y)) +
	geom_polygon(aes(fill = mean_shap), color = "grey60", alpha = 0.9) +
	facet_grid(color_cat ~ slope_cat, switch = "both") +
	scale_fill_gradient2(
		low = "#4575B4", mid = "white", high = "#D73027",
		midpoint = 0, name = "Mean SHAP\n(sign)"
	) +
	coord_equal() +
	theme_bw(base_family = "Helvetica", base_size = 10) +
	theme(
		panel.grid = element_blank(),
		axis.text = element_blank(),
		axis.ticks = element_blank(),
		axis.title = element_blank(),
		panel.border = element_rect(color = "grey60", linewidth = 0.4),
		strip.text.x = element_text(size = 8, face = "bold"),
		strip.text.y = element_text(size = 8, face = "bold", angle = 0),
		strip.background = element_rect(fill = "grey95", color = NA),
		legend.position = "right",
		plot.margin = margin(2, 5, 2, 5)
	) +
	labs(
		title = "Interpretation of colour × tilt combinations",
		subtitle = "Colour = sign of mean SHAP; Tilt = slope of SHAP–feature relationship"
	)

text_matrix <- tribble(
	~Colour, ~Tilt, ~Meaning,
	"🔴 Red (positive mean SHAP)", "↗ Upward", "Increasing positive effect across predictor range",
	"🔴 Red", "↘ Downward", "Positive at low values, saturating or reversing at high values",
	"🔴 Red", "◯ Flat", "Consistently positive, no gradient change",
	"⚪ Neutral", "↗ Upward", "Weak positive effect (context-dependent)",
	"⚪ Neutral", "↘ Downward", "Weak negative effect (context-dependent)",
	"⚪ Neutral", "◯ Flat", "No consistent relationship (irrelevant predictor)",
	"🔵 Blue (negative mean SHAP)", "↗ Upward", "Negative at low range, diminishing with predictor",
	"🔵 Blue", "↘ Downward", "Increasingly negative with predictor (strong filtering)",
	"🔵 Blue", "◯ Flat", "Uniformly negative across range"
)

table_theme <- ttheme_default(
	core = list(fg_params = list(cex = 0.8, fontfamily = "Helvetica")),
	colhead = list(fg_params = list(fontface = "bold", cex = 0.8, fontfamily = "Helvetica")),
	padding = unit(c(2, 4), "mm")
)

p_table <- tableGrob(text_matrix, rows = NULL, theme = table_theme)

# ---- 3. Combine ----
final_plot <- p_legend / as.ggplot(p_table) +
	plot_layout(heights = c(1, 1.1))

library(ggplot2)
library(dplyr)
library(tibble)
library(patchwork)
library(gridExtra)
library(grid)							
							
							