fia_clean <- fia_clean %>%
	# Drop missing trait values 
	filter(complete.cases(
		wmean_wood_density,wmean_bark_thickness, wmean_conduit_diam, wmean_leaf_n,
		wmean_shade_tolerance, wmean_specific_leaf_area, wmean_seed_dry_mass)) %>%
	# Normalize trait values and set directions appropriately 
	mutate(
		norm_wood_density = normalize(wmean_wood_density),
		norm_bark_thickness = normalize(wmean_bark_thickness),
		norm_conduit_diam = invert_normalized(normalize(wmean_conduit_diam)),
		norm_leaf_n = invert_normalized(normalize(wmean_leaf_n)),
		norm_specific_leaf_area = invert_normalized(normalize(wmean_specific_leaf_area)),
		norm_seed_dry_mass = normalize(wmean_seed_dry_mass),
		norm_shade_tolerance = normalize(wmean_shade_tolerance)) %>%
	# Scale around 0 and set format
	mutate(across(starts_with("norm_"), scale)) %>%
	mutate(across(starts_with("norm_"), as.numeric)) %>%
	# Calculate composite and resource use score
	rowwise() %>%
	mutate(
		composite_index = norm_wood_density + norm_bark_thickness + norm_conduit_diam + norm_leaf_n +
			norm_specific_leaf_area + norm_seed_dry_mass + norm_shade_tolerance) %>%
	ungroup() %>%
	mutate(resource_use_score = normalize(composite_index)) %>%
	# Drop redundant columns 
	dplyr::select(-norm_wood_density, -norm_bark_thickness, -norm_conduit_diam, -norm_leaf_n,
								-norm_specific_leaf_area, -norm_seed_dry_mass, -norm_shade_tolerance)





# Normalize a vector to a range [0, 1]
normalize <- function(x) {
	return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
}

# Invert normalized value for acquisitive traits
invert_normalized <- function(x) {
	return(1 - x)
}

# Normalize trait values appropriately
fia <- fia_clean %>%
	mutate(
		wmean_conduit_diam = -wmean_conduit_diam,
		wmean_leaf_n = -wmean_leaf_n,
		wmean_specific_leaf_area = -wmean_specific_leaf_area)


# Prepare trait matrix for PCA
trait_matrix <- fia %>%
	select(standage, wmean_wood_density, wmean_bark_thickness, wmean_conduit_diam, wmean_leaf_n, 
				 wmean_shade_tolerance, wmean_specific_leaf_area, wmean_seed_dry_mass) %>%
	filter(complete.cases(.)) %>%
	select(-standage) %>%
	as.matrix()

# Perform PCA
pca_result <- prcomp(trait_matrix, scale. = FALSE)
summary(pca_result)

# Extract scores for the first few principal components (e.g., PC1, PC2, PC3)
pc_scores <- pca_result$x[, 1:3]

# Calculate the proportion of variance explained by the first three PCs
var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
weights <- var_explained[1:3]

# Calculate the weighted sum of the first three principal components
composite_index <- rowSums(t(t(pc_scores) * weights))

# Add the composite index to the dataframe
fia <- fia %>%
	select(standage, wmean_wood_density, wmean_bark_thickness, wmean_conduit_diam, wmean_leaf_n, 
				 wmean_shade_tolerance, wmean_specific_leaf_area, wmean_seed_dry_mass)%>%
	filter(complete.cases(.)) %>%
	mutate(composite_index = composite_index)

# Normalize PC1 scores to [0, 1] range
fia <- fia %>%
	mutate(resource_use_score = normalize(composite_index))

plot_grid(
	
fia %>%
	# managed yes or no
	filter(standage < 500) %>%
	# Keep PID, stand age, and mean traits
	select(standage, resource_use_score) %>%
	# Plot indices against stand age 
	ggplot(aes(x=standage, y=resource_use_score)) +
	geom_hex() +
	scale_fill_viridis_c() +
	geom_smooth(method = "lm", colour = "red") +
	theme_bw() +
	theme(legend.position = "none") +
	ggtitle("PCA") +
	xlab("Stand age (years)") +
	ylab("Index values"),

fia_clean %>%
	# managed yes or no
	filter(standage < 500) %>%
	# Keep PID, stand age, and mean traits
	select(standage, resource_use_score) %>%
	# Plot indices against stand age 
	ggplot(aes(x=standage, y=resource_use_score)) +
	geom_hex() +
	scale_fill_viridis_c() +
	geom_smooth(method = "lm", colour = "red") +
	theme_bw() +
	theme(legend.position = "none") +
	ggtitle("Composite") +
	xlab("Stand age (years)") +
	ylab("Index values"),

nrow = 1, ncol = 2)

# Examine the PCA loadings
pca_loadings <- pca_result$rotation
print(pca_loadings)

# Check how each trait loads onto the first three principal components
