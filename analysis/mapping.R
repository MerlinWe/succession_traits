##########  MSc Diss. Forest Succession: Map ###### 

rm(list = ls()) # make sure environment is clean 
set.seed(42)    # set seed for reproducibility

# Load necessary libraries
library(sf)
library(hexbin)
library(viridis)
library(dggridR)
library(rnaturalearth)
library(ggspatial)
library(cowplot)
library(tidyverse)

# Check which device is running; set paths conditionally
node_name <- Sys.info()["nodename"]
path_in <- ifelse(node_name == "threadeast", 
									"/home/merlin/RDS_drive/merlin/data/fia_traits", 
									"/Volumes/ritd-ag-project-rd01pr-dmayn10/merlin/data/fia_traits")

path_out <- ifelse(node_name == "threadeast", 
									 "/home/merlin/traits_output", 
									 "/Users/serpent/Documents/MSc/Thesis/Code/analysis")

# Read data and prep as in rf analysis
data <- read_csv(paste0(path_in, "/plotlevel_data_2024-07-15.csv")) %>%
		select(PID_rep, standage, biome, managed, state, LAT, LON) %>%
	filter(complete.cases(.)) %>%
		filter(standage < quantile(standage, 0.9) |managed == 0) %>%
	filter(managed == 0) %>%
	select(-managed, -standage) 
data <- data %>% filter(biome %in% (
		data %>%
			count(biome) %>%
			filter(n > 100) %>%
			pull(biome)))

data <- data %>% dplyr::select(biome, state, LAT, LON) 
points_sf <- st_as_sf(data, coords = c("LON", "LAT"), crs = 4326)	

# First: USA mainland
world <- ne_countries(scale = "medium", returnclass = "sf")
usa_main <- ne_countries(scale = "medium", returnclass = "sf", country = "united states of america")
usa_main <- st_crop(usa_main, st_bbox(c(xmin = -125, xmax = -66.5, ymin = 24.5, ymax = 49.5), crs = st_crs(usa_main)))

points_sf <- st_as_sf(data, coords = c("LON", "LAT"), crs = 4326)
points_sf <- st_intersection(points_sf, usa_main)
us_data <- st_drop_geometry(points_sf) %>% 
	mutate(LAT = st_coordinates(points_sf)[,2], 
				 LON = st_coordinates(points_sf)[,1])

biome_colors <- c("Mediterranean woodlands" = "aquamarine4", 
									"Temperate broadleaf forests" = "darkgreen",
									"Temperate conifer forests" = "saddlebrown",
									"Temperate grasslands" = "yellow4",
									"Xeric shrublands" = "orange2")

biome_stats <- us_data %>%
	group_by(biome) %>%
	summarise(
		total_points = n(),
		mean_points_per_hex = mean(hexbin(LON, LAT, xbins = 80)@count))

legend_labels <- paste0(
	names(biome_colors), "\n",
	"n plots: ", biome_stats$total_points[match(names(biome_colors), biome_stats$biome)], "\n",
	"Mean plots per hexagon: ", round(biome_stats$mean_points_per_hex[match(names(biome_colors), biome_stats$biome)], 2))

usa <- ggplot() +
	geom_sf(data = usa_main, fill = "antiquewhite", color = "black", alpha = .5, linewidth = .3) +
	stat_bin_hex(data = us_data, aes(x = LON, y = LAT, fill = biome), 
							 color = "black", linewidth = .05, bins = 80, alpha = .8) +
	scale_fill_manual(values = biome_colors, name = "Biome", labels = legend_labels) +
	theme_bw(base_line_size = .3, base_rect_size = .7) +
	annotation_scale(location = "tr", width_hint = 0.7) +
	annotation_north_arrow(location = "tl", which_north = "true", style = north_arrow_fancy_orienteering) +
	coord_sf(xlim = c(-135, -61), ylim = c(24, 52)) +
	labs(x = NULL, y = NULL) +
	guides(color = guide_legend(override.aes = list(alpha = 1))) +
	theme(panel.background = element_rect(fill = "aliceblue"),
				legend.key.spacing.y =  unit(.5, "lines"),
				text = element_text(family = "sans", size = 6),
				legend.position = c(.999, 0.002),
				legend.justification = c("right", "bottom"),
				legend.title = element_blank(),
				legend.background = element_rect(fill = "white", colour = "black", linewidth = .3),
				plot.margin = margin(5, 5, 5, 5),
				legend.text = element_text(hjust = 0))

## Alaska 
alaska <- ne_countries(scale = "medium", returnclass = "sf", country = "united states of america")
alaska <- st_crop(alaska, st_bbox(c(xmin = -180, xmax = -130, ymin = 50, ymax = 72), crs = st_crs(usa_main)))

# Transform the Data
points_sf <- st_as_sf(data, coords = c("LON", "LAT"), crs = 4326)
points_sf <- st_intersection(points_sf, alaska)
al_data <- st_drop_geometry(points_sf) %>% 
	mutate(LAT = st_coordinates(points_sf)[,2], 
				 LON = st_coordinates(points_sf)[,1])

biome_colors <- c("Boreal forests or taiga" = "brown", 
									"Temperate conifer forests" = "saddlebrown",
									"Tundra" = "azure3")

biome_stats <- al_data %>%
	group_by(biome) %>%
	summarise(
		total_points = n(),
		mean_points_per_hex = mean(hexbin(LON, LAT, xbins = 50)@count))

legend_labels <- paste0(
	names(biome_colors), "\n",
	"n plots: ", biome_stats$total_points[match(names(biome_colors), biome_stats$biome)], "\n",
	"Mean plots per hexagon: ", round(biome_stats$mean_points_per_hex[match(names(biome_colors), biome_stats$biome)], 2))

# Plot using geom_point and geom_density2d
alaska <- ggplot() +
	geom_sf(data = alaska, fill = "antiquewhite", color = "black", alpha = .5, linewidth = .2) +
	stat_bin_hex(data = al_data, aes(x = LON, y = LAT, fill = biome), 
							 color = "black", linewidth = .05, bins = 50, alpha = .8) +
	scale_fill_manual(values = biome_colors, name = "Biome", labels = legend_labels) +
	theme_bw(base_line_size = .3, base_rect_size = .7) +
	coord_sf(xlim = c(-170, -131.5), ylim = c(54, 71)) +
	labs(x = NULL, y = NULL) +
	guides(fill = guide_legend(ncol = 1, override.aes = list(alpha = 1))) +
	theme(panel.background = element_rect(fill = "aliceblue"),
				legend.key.spacing.y =  unit(.3, "lines"),
				text = element_text(family = "sans", size = 5),
				legend.position = "top",
				legend.justification = "left",
				legend.title = element_blank(),
				legend.background = element_rect(fill = "white", colour = "black", linewidth = .3),
				plot.margin = margin(0,0,0,0),
				legend.box.margin = margin(0,0,-11.3,0),
				axis.title = element_blank(),
				axis.text = element_blank(),
				axis.ticks = element_blank(),
				legend.key.size = unit(4, "mm"),
				plot.background = element_rect(fill = "transparent", colour = "transparent"),
				legend.text = element_text(hjust = 0))

# Combine plots
full <- ggdraw() +
	draw_plot(usa) +
	draw_plot(alaska, x = -.1068, y = 0.125, width = 0.47, height = 0.47)

ggsave(filename = paste0(path_out, "/plots/map.png"),
			 plot = full, 
			 bg = "white",
			 width = 200, 
			 height = 120, 
			 units = "mm", 
			 dpi = 600)
