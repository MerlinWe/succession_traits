##### Descriptive Tables for appendix #####

library(gt)
library(ggplot2)
library(dplyr)
library(magick)
library(gridExtra)

## >>>>>>>>>>>>>>> Trait descriptives <<<<<<<<<<<<<<<

traits <- data %>% 
	select(wood_density, bark_thickness, conduit_diam, leaf_n, specific_leaf_area, seed_dry_mass, shade_tolerance, height)

# Function to create histogram for a variable
create_histogram <- function(traits, var, path) {
	p <- ggplot(traits, aes_string(x = var)) +
		geom_histogram(binwidth = 0.1, fill = "forestgreen", color = "black") +
		theme_void() +
		theme(axis.title = element_blank(),
					axis.text = element_blank(),
					axis.ticks = element_blank())
	ggsave(path, plot = p, width = 200, height = 130, units = "mm", dpi = 600)
}

# Directory to save histograms
histogram_dir <- "/Users/serpent/Documents/MSc/Thesis/Code/analysis/descriptives/histograms/"
dir.create(histogram_dir, showWarnings = FALSE)

# Create a list of histograms and save them in the specified directory
histogram_paths <- sapply(names(traits), function(var) {
	temp_path <- file.path(histogram_dir, paste0(var, "_hist.png"))
	create_histogram(traits, var, temp_path)
	temp_path
})

names(histogram_paths) <- names(traits)

# Summarize the traits manually
summary_stats <- data.frame(
	Variable = names(traits),
	Mean = round(sapply(traits, mean), 2),
	Min = round(sapply(traits, min), 2),
	Max = round(sapply(traits, max), 2),
	Unit = "log-scale")

# Function to format variable names
format_variable_name <- function(name) {
	name <- gsub("_", " ", name)
	name <- tolower(name)
	name <- sub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", name, perl = TRUE)
	return(name)
}

summary_stats <- summary_stats %>%
	mutate(Trait = sapply(Variable, format_variable_name)) %>%
	mutate(Trait = ifelse(Trait == "Conduit diam", "Conduit diameter", Trait))

# Create the gt table
table <- gt(summary_stats) %>%
	tab_header(
		title = "Functional Tree Traits Used") %>%
	text_transform(
		locations = cells_body(columns = "Variable"),
		fn = function(x) {
			lapply(x, function(var) {
				path <- histogram_paths[[var]]
				paste0('<img src="', path, '" style="width: 200px; height: 60px;" />')
			})
		}) %>%
	cols_label(
		Variable = "Distribution",
		Mean = "Mean",
		Min = "Min",
		Max = "Max",
		Unit = "Unit",
		Trait = "Trait") %>%
	cols_move_to_start(columns = "Trait") %>%
	tab_options(
		table.border.top.style = "solid",
		table.border.top.width = px(2),
		table.border.top.color = "black",
		table.border.bottom.style = "solid",
		table.border.bottom.width = px(2),
		table.border.bottom.color = "black",
		heading.border.bottom.style = "solid",
		heading.border.bottom.width = px(2),
		heading.border.bottom.color = "black",
		column_labels.border.bottom.style = "none",
		table_body.hlines.style = "none",
		table_body.vlines.style = "none") %>%
	tab_style(
		style = cell_borders(
			sides = "top", color = "black", weight = px(2)),
		locations = cells_column_labels(columns = everything())) %>%
	tab_style(
		style = cell_borders(
			sides = "bottom", color = "black", weight = px(2)),
		locations = cells_column_labels(columns = everything())) %>%
	tab_style(
		style = cell_borders(
			sides = "bottom", color = "black", weight = px(2)),
		locations = cells_body(
			rows = c(nrow(summary_stats))))

# Save the table as an HTML file
html_path <- "/Users/serpent/Documents/MSc/Thesis/Code/analysis/descriptives/trait_descriptives.html"
gtsave(table, html_path)

## >>>>>>>>>>>>>>> Feature descriptives <<<<<<<<<<<<<<<

features <- data %>% 
	select(standage, temp_pc, soil_pc, rain_pc, elevation, soil_ph)

# Function to create histogram for a variable
create_histogram <- function(features, var, path) {
	p <- ggplot(features, aes_string(x = var)) +
		geom_histogram(binwidth = 0.1, fill = "navyblue", color = "navyblue") +
		theme_void() +
		theme(axis.title = element_blank(),
					axis.text = element_blank(),
					axis.ticks = element_blank())
	ggsave(path, plot = p, width = 200, height = 130, units = "mm", dpi = 600)
}

# Directory to save histograms
histogram_dir <- "/Users/serpent/Documents/MSc/Thesis/Code/analysis/descriptives/histograms/"
dir.create(histogram_dir, showWarnings = FALSE)

# Create a list of histograms and save them in the specified directory
histogram_paths <- sapply(names(features), function(var) {
	temp_path <- file.path(histogram_dir, paste0(var, "_hist.png"))
	create_histogram(features, var, temp_path)
	temp_path
})

names(histogram_paths) <- names(features)

# Summarize the features manually
summary_stats <- data.frame(
	Variable = names(features),
	Mean = round(sapply(features, mean), 2),
	Min = round(sapply(features, min), 2),
	Max = round(sapply(features, max), 2),
	Unit = c("years", "principal\ncomponent", "principal\ncomponent", "principal\ncomponent", "meters", "pH"))

# Function to format variable names
format_variable_name <- function(name) {
	name <- gsub("_", " ", name)
	name <- tolower(name)
	name <- sub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", name, perl = TRUE)
	return(name)
}

summary_stats <- summary_stats %>%
	mutate(Feature = sapply(Variable, format_variable_name)) %>%
	mutate(Feature = case_when(
		Feature == "Temp pc" ~ "Temperature",
		Feature == "Soil pc" ~ "Soil Water Retention",
		Feature == "Rain pc" ~ "Precipitation",
		
		Feature == "Elevation" ~ "Elevation",
		Feature == "Soil ph" ~ "Soil pH",
		Feature == "Standage" ~ "Standage",
		TRUE ~ NA_character_))

# Create the gt table
table <- gt(summary_stats) %>%
	tab_header(
		title = "Environmental Features Used") %>%
	text_transform(
		locations = cells_body(columns = "Variable"),
		fn = function(x) {
			lapply(x, function(var) {
				path <- histogram_paths[[var]]
				paste0('<img src="', path, '" style="width: 200px; height: 60px;" />')
			})
		}) %>%
	cols_label(
		Variable = "Distribution",
		Mean = "Mean",
		Min = "Min",
		Max = "Max",
		Unit = "Unit",
		Feature = "Feature") %>%
	cols_move_to_start(columns = "Feature") %>%
	tab_options(
		table.border.top.style = "solid",
		table.border.top.width = px(2),
		table.border.top.color = "black",
		table.border.bottom.style = "solid",
		table.border.bottom.width = px(2),
		table.border.bottom.color = "black",
		heading.border.bottom.style = "solid",
		heading.border.bottom.width = px(2),
		heading.border.bottom.color = "black",
		column_labels.border.bottom.style = "none",
		table_body.hlines.style = "none",
		table_body.vlines.style = "none") %>%
	tab_style(
		style = cell_borders(
			sides = "top", color = "black", weight = px(2)),
		locations = cells_column_labels(columns = everything())) %>%
	tab_style(
		style = cell_borders(
			sides = "bottom", color = "black", weight = px(2)),
		locations = cells_column_labels(columns = everything())) %>%
	tab_style(
		style = cell_borders(
			sides = "bottom", color = "black", weight = px(2)),
		locations = cells_body(
			rows = c(nrow(summary_stats))))

# Save the table as an HTML file
html_path <- "/Users/serpent/Documents/MSc/Thesis/Code/analysis/descriptives/feature_descriptives.html"
gtsave(table, html_path)

## >>>>>>>>>>>>>>> Biome descriptives <<<<<<<<<<<<<<<

summary_stats <- data %>%
	select(biome_boreal_forests_or_taiga, biome_flooded_grasslands, biome_mediterranean_woodlands, 
				 biome_temperate_broadleaf_forests, biome_temperate_conifer_forests, biome_temperate_grasslands,  
				 biome_tundra, biome_xeric_shrublands) %>%
	rename_with(~ str_replace_all(., "biome_", ""), everything()) %>%
	rename_with(~ str_replace_all(., "_", " "), everything()) %>%
	rename_with(~ str_to_title(.), everything()) %>%
	summarise(across(everything(), list(
		Found = ~sum(. == 1),
		'Not found' = ~sum(. == 0)))) %>%
	pivot_longer(cols = everything(), names_to = "Variable", values_to = "Count") %>%
	separate(Variable, into = c("Variable", "Value"), sep = "_") %>%
	pivot_wider(names_from = "Value", values_from = "Count")

table <- gt(summary_stats) %>%
	tab_header(
		title = "Biome Features Used"
	) %>%
	cols_label(
		Variable = "Variable",
		Found = "Count of ocurrences",
		'Not found' = "Count of non-occurrences"
	) %>%
	tab_options(
		table.border.top.style = "solid",
		table.border.top.width = px(2),
		table.border.top.color = "black",
		table.border.bottom.style = "solid",
		table.border.bottom.width = px(2),
		table.border.bottom.color = "black",
		heading.border.bottom.style = "solid",
		heading.border.bottom.width = px(2),
		heading.border.bottom.color = "black",
		column_labels.border.bottom.style = "none",
		table_body.hlines.style = "none",
		table_body.vlines.style = "none") %>%
	tab_style(
		style = cell_borders(
			sides = "top", color = "black", weight = px(2)),
		locations = cells_column_labels(columns = everything())) %>%
	tab_style(
		style = cell_borders(
			sides = "bottom", color = "black", weight = px(2)),
		locations = cells_column_labels(columns = everything())) %>%
	tab_style(
		style = cell_borders(
			sides = "bottom", color = "black", weight = px(2)),
		locations = cells_body(
			rows = c(nrow(summary_stats))))

# Save the table as an HTML file
html_path <- "/Users/serpent/Documents/MSc/Thesis/Code/analysis/descriptives/biome_descriptives.html"
gtsave(table, html_path)

