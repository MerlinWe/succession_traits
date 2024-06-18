## Read most recent fia_clean data 
fia_clean <- list.files("/Volumes/ritd-ag-project-rd01pr-dmayn10/merlin/data/fia_traits", full.names = TRUE) %>%
	file.info() %>%
	as_tibble(rownames = "file") %>%
	arrange(desc(mtime)) %>%
	slice(1) %>%
	pull(file) %>%
	read_csv()

########## Calculate Resource Use Score using Euclidean distance ##########

# Function to normalize a vector to a range [0, 1]
normalize <- function(x) {
	return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
}

# Function to invert normalized value for acquisitive traits
invert_normalized <- function(x) {
	return(1 - x)
}

# Normalize trait values appropriately
fia_clean <- og %>%
	filter(complete.cases(wmean_wood_density,wmean_bark_thickness, wmean_conduit_diam, wmean_leaf_n,
												wmean_shade_tolerance, wmean_specific_leaf_area, wmean_seed_dry_mass)) %>%
	mutate(
		norm_wood_density = normalize(wmean_wood_density),
		norm_bark_thickness = normalize(wmean_bark_thickness),
		norm_conduit_diam = invert_normalized(normalize(wmean_conduit_diam)),
		norm_leaf_n = invert_normalized(normalize(wmean_leaf_n)),
		norm_specific_leaf_area = invert_normalized(normalize(wmean_specific_leaf_area)),
		norm_seed_dry_mass = normalize(wmean_seed_dry_mass),
		norm_shade_tolerance = normalize(wmean_shade_tolerance)) %>%
	# Scale and set format
	mutate(across(starts_with("norm_"), scale)) %>%
	mutate(across(starts_with("norm_"), as.numeric)) %>%
	# Calculate score 
	rowwise() %>%
	mutate(
		composite_index = norm_wood_density + norm_bark_thickness + norm_conduit_diam + norm_leaf_n +
			norm_specific_leaf_area + norm_seed_dry_mass + norm_shade_tolerance) %>%
	ungroup() %>%
	mutate(resource_use_score = normalize(composite_index)) %>%
	dplyr::select(-norm_wood_density, -norm_bark_thickness, -norm_conduit_diam, -norm_leaf_n,
									-norm_specific_leaf_area, -norm_seed_dry_mass, -norm_shade_tolerance)

fia_clean %>%
	filter(standage < 500) %>%
	filter(resource_use_score > 0 & resource_use_score < 1) %>%
	select(standage, resource_use_score, biome) %>%
	ggplot(aes(x=standage, y=resource_use_score)) +
	geom_hex() +
	scale_fill_viridis_c() +
	geom_smooth(method = "lm", colour = "red") +
	theme_bw() +
	theme(legend.position = "none") +
	xlab("Stand age (years)") +
	ylab("Index values")
