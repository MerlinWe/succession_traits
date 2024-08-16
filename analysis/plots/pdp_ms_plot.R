library(tidyverse)

pdp_data <- read_csv("/Volumes/ritd-ag-project-rd01pr-dmayn10/merlin/traits_output/tables/pdp_data.csv")

# Quantile labels
quantile_labels <- list(
	temp_pc = c("Cold Temperatures (Lower 25%)", "Warm Temperatures (Upper 25%)"),
	rain_pc = c("Low Precipitation (Lower 25%)", "High Precipitation (Upper 25%)")
)

# Clean labels for traits
trait_labels <- c(
	wood_density = "Wood\nDensity",
	bark_thickness = "Bark\nThickness",
	conduit_diam = "Conduit\nDiameter",
	leaf_n = "Leaf\nNitrogen",
	specific_leaf_area = "Specific\nLeaf Area",
	seed_dry_mass = "Seed\nDry Mass",
	shade_tolerance = "Shade\nTolerance",
	height = "Tree\nHeight"
)

# Function for pdp plotting
create_pdp_plot <- function(data, levels, colors, labels, x_lab, y_lab, show_y_strip_labels) {
	data %>%
		mutate(group = factor(group, levels = levels),
					 variable = factor(variable, levels = c("Temperature (PC)", "Precipitation (PC)"), labels = c("Temperature", "Precipitation"))) %>%
		ggplot(aes(x = standage, y = yhat, color = group, shape = group)) +
		geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, linewidth = 1) +
		scale_color_manual(values = setNames(colors, labels)) +
		scale_fill_manual(values = setNames(colors, labels)) +
		labs(x = x_lab, y = y_lab, color = "Quantile", shape = "Quantile", fill = "Quantile") +
		theme_bw(base_line_size = .3, base_rect_size = .5) +
		theme(
			legend.position = "top",
			legend.box = "horizontal",
			legend.title = element_blank(),
			legend.text = element_text(size = 10),
			text = element_text(family = "sans", size = 12),
			strip.background.x = element_rect(fill = "white", color = "black", linewidth = .75),
			strip.background.y = if (show_y_strip_labels) element_rect(fill = "white", color = "black", linewidth = .75) else element_blank(),
			strip.text.y = if (show_y_strip_labels) element_text() else element_blank(),
			strip.text.x = element_text(face = "bold")
		) +
		guides(
			color = guide_legend(nrow = 2, ncol = 2),
			shape = guide_legend(nrow = 2, ncol = 2)
		) +
		facet_grid(trait ~ variable, scales = "free", labeller = labeller(trait = trait_labels))
}

# Generate the PDP plot for temperature and precipitation
pdp_plot <- create_pdp_plot(
	data = pdp_data %>% filter(variable %in% c("Temperature (PC)", "Precipitation (PC)")),
	levels = c(quantile_labels$temp_pc, quantile_labels$rain_pc),
	colors = c("darkslateblue", "darkred", "tan4", "mediumseagreen"),
	labels = c(quantile_labels$temp_pc, quantile_labels$rain_pc),
	x_lab = "Standage (years)",
	y_lab = "Predicted Trait Value (log-scaled)",
	show_y_strip_labels = TRUE
)

# Save the plot
ggsave(filename = paste0("/Users/serpent/Desktop/pdp_ms_plot.png"),
			 plot = pdp_plot, 
			 bg = "white",
			 width = 200, 
			 height = 220, 
			 units = "mm", 
			 dpi = 600)