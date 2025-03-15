##########  MSc Diss. Forest Succession: Map ###### 

rm(list = ls()) # make sure environment is clean 
set.seed(42)    # set seed for reproducibility

# Load necessary libraries
library(sf)
library(viridis)
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
									 "/Users/serpent/Documents/MSc/Thesis/Code/output/plots")

# Read data and prep as in RF analysis
data <- read_csv(paste0(path_in, "/plotlevel_data_2024-07-15.csv")) %>%
	select(PID_rep, standage, biome, managed, state, LAT, LON) %>%
	filter(complete.cases(.)) %>%
	filter(standage < quantile(standage, 0.9) | managed == 0) %>%
	filter(managed == 0) %>%
	select(-managed, -standage) 

data <- data %>% filter(biome %in% (
	data %>%
		count(biome) %>%
		filter(n > 100) %>%
		pull(biome)))

# Get spatial points
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

usa <- ggplot() +
	geom_sf(data = usa_main, fill = "beige", color = "black", alpha = .5, linewidth = .5) +
	stat_density_2d_filled(data = us_data, aes(x = LON, y = LAT, fill = after_stat(level)),
												 contour = TRUE, breaks = c(5,10, 50, 100, 200, 300, 400, 500),
												 geom = "polygon", contour_var = "count", alpha = .7) +
	geom_sf(data = usa_main, fill = NA, color = "black", alpha = .5, linewidth = .5) +
	scale_fill_manual(values = viridis::viridis(6), name = "Number of Plots:",
										labels = custom_labels <- c("5-10", "11-50", "51-100", "101-200", "201-300", "301-400", "401-500")) +
	theme_bw(base_line_size = .3, base_rect_size = .7) +
	annotation_scale(location = "tr", width_hint = 0.3) +
	annotation_north_arrow(location = "tl", which_north = "true", style = north_arrow_fancy_orienteering) +
	coord_sf(xlim = c(-131, -61), ylim = c(24, 52)) +
	labs(x=NULL, y=NULL) +
	theme(panel.background = element_rect(fill = "aliceblue"),
				legend.key.spacing.y =  unit(.5, "lines"),
				text = element_text(family = "sans", size = 8),
				axis.text = element_text(family = "sans", size = 8),
				legend.position = c(.999, 0.002),
				legend.justification = c("right", "bottom"),
				legend.background = element_rect(fill = "white", colour = "black", linewidth = .3),
				legend.text = element_text(hjust = 0),
				plot.margin = margin(t=1, r=10, b=1, l=1))

## Alaska 
alaska <- ne_countries(scale = "medium", returnclass = "sf", country = "united states of america")
alaska <- st_crop(alaska, st_bbox(c(xmin = -180, xmax = -130, ymin = 50, ymax = 72), crs = st_crs(usa_main)))

# Transform the Data
points_sf <- st_as_sf(data, coords = c("LON", "LAT"), crs = 4326)
points_sf <- st_intersection(points_sf, alaska)
al_data <- st_drop_geometry(points_sf) %>% 
	mutate(LAT = st_coordinates(points_sf)[,2], 
				 LON = st_coordinates(points_sf)[,1])

# Plot using geom_point and geom_density2d
alaska <- ggplot() +
	geom_sf(data = alaska, fill = "beige", color = "black", alpha = .5, linewidth = .4) +
	stat_density_2d_filled(data = al_data, aes(x = LON, y = LAT, fill = after_stat(level)),
												 contour = TRUE, breaks = c(5,10, 50, 100), geom = "polygon", contour_var = "count", alpha = .7) +
	geom_sf(data = alaska, fill = NA, color = "black", alpha = .5, linewidth = .4) +
	scale_fill_manual(values = viridis::viridis(6), name = "Number of Plots",
										labels = custom_labels <- c("1-10", "11-50", "51-100")) +
	theme_bw(base_line_size = .3, base_rect_size = .7) +
	coord_sf(xlim = c(-170, -131.5), ylim = c(54, 71)) +
	labs(x=NULL, y=NULL) +
	guides(fill = guide_legend(ncol = 1, override.aes = list(alpha = 1))) +
	theme(panel.background = element_rect(fill = "aliceblue"),
				text = element_text(family = "sans", size = 5),
				plot.margin = margin(0,0,0,0),
				legend.position = "none",
				axis.title = element_blank(),
				axis.text = element_blank(),
				axis.ticks = element_blank(),
				plot.background = element_rect(fill = "transparent", colour = "transparent"))

# Combine plots
full <- ggdraw() +
	draw_plot(usa) +
	draw_plot(alaska, x = -.015, y = 0.1135, width = 0.28, height = 0.28)

ggsave(filename = paste0(path_out, "/fig1.png"),
			 plot = full, 
			 bg = "white",
			 width = 200, 
			 height = 120, 
			 units = "mm", 
			 dpi = 600)
