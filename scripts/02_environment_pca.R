################################################################################
## succession_traits: 02 — Environmental PCA

## Identifies underlying environmental axes from climate and soil variables
## using PCA with varimax rotation. Projects plots onto rotated axes and
## writes the analysis-ready dataset used by all downstream modelling scripts.

################################################################################

rm(list = ls())
set.seed(42)

# ── Libraries ─────────────────────────────────────────────────────────────────
library(psych)        
library(factoextra)   
library(cowplot)      
library(tidyverse)

# ── Paths ─────────────────────────────────────────────────────────────────────
PATH_IN     <- "data_processed/plotlevel_data.rds"
PATH_OUT    <- "data_processed/fia_traits_clean.rds"
PATH_PLOTS  <- "figures/supplementary"
PATH_PCA    <- "data_processed/pca_rotation.rds"

export <- TRUE

# ── Helper functions ──────────────────────────────────────────────────────────
source("scripts/functions.R")

# ══════════════════════════════════════════════════════════════════════════════
# 1. Load data
# ══════════════════════════════════════════════════════════════════════════════

data <- read_rds(PATH_IN)

message(sprintf("Input: %d plots, %d columns", nrow(data), ncol(data)))


# ══════════════════════════════════════════════════════════════════════════════
# 2. Prepare climate variable matrix for PCA
# ══════════════════════════════════════════════════════════════════════════════

# Elevation and soil pH are excluded from the PCA:
#   - Elevation: included directly as a predictor; no latent factor assumed
#   - Soil pH: orthogonal to water-retention/sand axis; including it would
#     conflate two distinct soil dimensions into one component
# All remaining temperature, precipitation, and soil texture variables are
# included to capture the major axes of climatic and edaphic variation.

CLIMATE_VARS <- c(
	# Temperature
	"annual_mean_temperature",
	"max_temperature_of_warmest_month",
	"min_temperature_of_coldest_month",
	"mean_temperature_of_coldest_quarter",
	"mean_temperature_of_warmest_quarter",
	# Precipitation
	"annual_precipitation",
	"precipitation_of_driest_quarter",
	"precipitation_of_wettest_quarter",
	"precipitation_of_coldest_quarter",
	"precipitation_of_warmest_quarter",
	# Soil texture
	"sand_content_015cm",
	"sand_content_060cm",
	"water_capacity_015cm",
	"water_capacity_060cm"
)

climate_mat <- data %>%
	dplyr::select(all_of(CLIMATE_VARS)) %>%
	# Drop any rows with missing values for PCA (should be none after script 01)
	{ stopifnot("Missing values in climate matrix" = complete.cases(.) %>% all()); . }


# ══════════════════════════════════════════════════════════════════════════════
# 3. PCA suitability checks
# ══════════════════════════════════════════════════════════════════════════════

message("\n── KMO test ───────────────────────────────────────────────────────")
kmo_result <- KMO(climate_mat)
print(kmo_result)

# KMO > 0.8 = "meritorious"; > 0.9 = "marvelous" (Kaiser 1974)
if (kmo_result$MSA < 0.8) {
	warning(sprintf(
		"KMO = %.3f is below 0.8",
		kmo_result$MSA
	))
}

message("\n── Bartlett's test of sphericity ──────────────────────────────────")
# n must be number of observations, not variables
bartlett_result <- cortest.bartlett(cor(climate_mat), n = nrow(climate_mat))
print(bartlett_result)

if (bartlett_result$p.value > 0.05) {
	warning("Bartlett's test non-significant — correlation matrix may be identity-like.")
}


# ══════════════════════════════════════════════════════════════════════════════
# 4. Fit PCA
# ══════════════════════════════════════════════════════════════════════════════

climate_pca <- prcomp(climate_mat, scale. = TRUE, center = TRUE)

message("\n── PCA variance explained ─────────────────────────────────────────")
pca_summary <- summary(climate_pca)
print(pca_summary)

# Confirm at least 3 components have eigenvalue > 1 (Kaiser criterion)
eigenvalues <- climate_pca$sdev^2
n_components_kaiser <- sum(eigenvalues > 1)
message(sprintf(
	"\nComponents with eigenvalue > 1: %d (retaining 3 for rotation)",
	n_components_kaiser
))
if (n_components_kaiser < 3) {
	warning("Fewer than 3 components meet Kaiser criterion — reconsider retaining 3.")
}


# ══════════════════════════════════════════════════════════════════════════════
# 5. Varimax rotation
# ══════════════════════════════════════════════════════════════════════════════

varimax_result   <- varimax(climate_pca$rotation[, 1:3])
varimax_loadings <- varimax_result$loadings   # rotated loading matrix (14 × 3)

message("\n── Varimax-rotated loadings ────────────────────────────────────────")
print(varimax_loadings, cutoff = 0.2)   # suppress near-zero loadings for clarity

climate_scaled   <- scale(climate_mat)   # same centering/scaling as prcomp
rotated_scores   <- climate_scaled %*% varimax_loadings   # n_plots × 3 matrix

# Assign interpretable names based on dominant loadings
colnames(rotated_scores) <- c("temp_pc", "soil_pc", "rain_pc")

# Sanity check: rotated scores should be uncorrelated
score_cors <- cor(rotated_scores)
max_offdiag <- max(abs(score_cors[upper.tri(score_cors)]))
message(sprintf(
	"\nMax off-diagonal correlation among rotated scores: %.4f (should be low)",
	max_offdiag
))

# Save rotation matrix so new data can be projected onto same axes later
saveRDS(
	list(
		loadings       = varimax_loadings,
		center         = attr(climate_scaled, "scaled:center"),
		scale          = attr(climate_scaled, "scaled:scale"),
		climate_vars   = CLIMATE_VARS,
		component_names = colnames(rotated_scores)
	),
	file = PATH_PCA
)
message(sprintf("PCA rotation saved to: %s", PATH_PCA))


# ══════════════════════════════════════════════════════════════════════════════
# 6. Attach rotated scores to data and finalise analysis dataset
# ══════════════════════════════════════════════════════════════════════════════

data <- data %>%
	bind_cols(as_tibble(rotated_scores))

# Final column selection for modelling
# Trait columns: wmean_ prefix stripped for cleaner model formulae
# Retaining LAT/LON for spatial diagnostics but not as model predictors
data_clean <- data %>%
	dplyr::select(
		# Traits (log-transformed, z-scored CWMs from script 01)
		starts_with("wmean_"),
		
		# Succession predictor
		standage,
		
		# Environmental predictors (rotated PCs + direct controls)
		temp_pc, soil_pc, rain_pc,
		elevation, soil_ph_015cm,
		
		# Biome dummies (for stratification and leaf-type derivation)
		starts_with("biome_"),
		
		# Spatial coordinates (diagnostics only)
		LAT, LON,
		
		# Plot metadata (kept for traceability)
		PID, PID_rep, rep_measure, PID_measure, state, INVYR,
		
		# Shade tolerance data quality flag (from script 01)
		shade_ba_coverage, shade_low_coverage
	) %>%
	
	# Strip wmean_ prefix: cleaner formulae in modelling scripts
	rename_with(~ str_remove(., "^wmean_"), starts_with("wmean_")) %>%
	
	# Strip _015cm suffix from soil pH (already selected single depth above)
	rename_with(~ str_remove(., "_015cm$"), ends_with("_015cm")) %>%
	
	# Drop any rows with missing values in modelling columns
	# (log only — should be minimal after script 01 complete.cases filtering)
	{ 
		n_before <- nrow(.)
		out <- filter(., complete.cases(dplyr::select(.,
																									bark_thickness, conduit_diam, leaf_density, leaf_k,
																									root_depth, seed_dry_mass, specific_leaf_area, shade_tolerance,
																									height, standage, temp_pc, soil_pc, rain_pc, elevation, soil_ph
		)))
		n_dropped <- n_before - nrow(out)
		if (n_dropped > 0) message(sprintf(
			"complete.cases filter dropped %d rows (%.2f%% of data)",
			n_dropped, 100 * n_dropped / n_before
		))
		out
	}

# ══════════════════════════════════════════════════════════════════════════════
# 7. Supplementary Figure S-1: PCA validation compound figure
# ══════════════════════════════════════════════════════════════════════════════
# Four-panel PCA block (scree + varimax contributions) combined with
# PC scores by biome and predictor correlation matrix.
# Exported directly to figures/supplementary/pca/S1_pca_validation.png.

library(patchwork)
source("scripts/plot_theme.R")

DIR_S1 <- "figures/supplementary"
dir.create(DIR_S1, showWarnings = FALSE, recursive = TRUE)

# ── Scree plot ────────────────────────────────────────────────────────────────
scree_plot <- tibble(
	component  = paste0("PC", seq_along(eigenvalues)),
	eigenvalue = eigenvalues,
	pct_var    = 100 * eigenvalues / sum(eigenvalues)
) %>%
	slice(1:10) %>%
	mutate(component = fct_inorder(component)) %>%
	ggplot(aes(x = component, y = pct_var, group = 1)) +
	geom_col(fill = "#1B7E74", colour = "black", alpha = 0.8) +
	geom_point() +
	geom_line() +
	geom_vline(xintercept = 3.5, linetype = "dashed", colour = "red") +
	labs(y = "Variance explained (%)", x = "Principal component") +
	theme_succession(base_size = 8)

# ── Varimax contribution plots ────────────────────────────────────────────────
abbreviations <- c(
	"annual_mean_temperature"             = "AnnMeanTemp",
	"max_temperature_of_warmest_month"    = "MaxTempWarmest",
	"min_temperature_of_coldest_month"    = "MinTempColdest",
	"mean_temperature_of_coldest_quarter" = "MeanTempColdest",
	"mean_temperature_of_warmest_quarter" = "MeanTempWarmest",
	"annual_precipitation"                = "AnnPrecip",
	"precipitation_of_driest_quarter"     = "PrecipDriest",
	"precipitation_of_wettest_quarter"    = "PrecipWettest",
	"precipitation_of_coldest_quarter"    = "PrecipColdest",
	"precipitation_of_warmest_quarter"    = "PrecipWarmest",
	"sand_content_015cm"                  = "Sand015",
	"sand_content_060cm"                  = "Sand060",
	"water_capacity_015cm"                = "WaterCap015",
	"water_capacity_060cm"                = "WaterCap060"
)

varimax_contributions <- varimax_contrib(varimax_loadings)

var_contrib_pc1 <- tibble(
	name    = rownames(varimax_contributions),
	contrib = varimax_contributions[, 1]
) %>%
	create_contrib_plot("PC1 (Temperature)", "firebrick4", abbreviations)

var_contrib_pc2 <- tibble(
	name    = rownames(varimax_contributions),
	contrib = varimax_contributions[, 2]
) %>%
	create_contrib_plot("PC2 (Soil water retention)", "tan4", abbreviations)

var_contrib_pc3 <- tibble(
	name    = rownames(varimax_contributions),
	contrib = varimax_contributions[, 3]
) %>%
	create_contrib_plot("PC3 (Precipitation)", "dodgerblue3", abbreviations)

# ── PC scores by biome ────────────────────────────────────────────────────────
# Recover biome label from dummy columns
biome_cols <- names(data_clean) %>% str_subset("^biome_")
data_clean <- data_clean %>%
	mutate(biome_label = biome_cols %>%
				 	map_dfc(~ if_else(data_clean[[.x]] == 1L,
				 										str_remove(.x, "^biome_") %>%
				 											str_replace_all("_", " ") %>%
				 											str_to_title(),
				 										NA_character_)) %>%
				 	reduce(coalesce))

p_biome <- data_clean %>%
	pivot_longer(c(temp_pc, soil_pc, rain_pc),
							 names_to = "pc", values_to = "score") %>%
	mutate(
		pc = recode(pc,
								temp_pc = "Temperature PC",
								soil_pc = "Soil water retention PC",
								rain_pc = "Precipitation PC"),
		pc          = factor(pc, levels = c("Temperature PC",
																				"Soil water retention PC",
																				"Precipitation PC")),
		biome_label = fct_reorder(biome_label, score, median, na.rm = TRUE)
	) %>%
	ggplot(aes(biome_label, score, fill = pc)) +
	geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
	geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3, linewidth = 0.3) +
	facet_wrap(~ pc, ncol = 1, scales = "free_x") +
	coord_flip() +
	scale_fill_manual(values = c(
		"Temperature PC"          = "#E69F00",
		"Soil water retention PC" = "#009E73",
		"Precipitation PC"        = "#56B4E9"
	)) +
	labs(x = NULL, y = "Rotated PC score") +
	theme_succession(base_size = 8) +
	theme(
		legend.position = "none",
		strip.text      = element_text(face = "bold", size = 7,
																	 margin = margin(t = 4, b = 4))
	)

# ── Predictor correlation matrix ──────────────────────────────────────────────
p_cors <- c("standage", "temp_pc", "soil_pc",
						"rain_pc",  "elevation", "soil_ph") %>%
	{ cor(data_clean[.], use = "pairwise.complete.obs") } %>%
	as_tibble(rownames = "v1") %>%
	pivot_longer(-v1, names_to = "v2", values_to = "r") %>%
	mutate(
		flag = abs(r) > 0.5 & v1 != v2,
		v1   = recode(v1, !!!ALL_PREDICTOR_LABELS),
		v2   = recode(v2, !!!ALL_PREDICTOR_LABELS)
	) %>%
	ggplot(aes(v1, v2, fill = r)) +
	geom_tile(colour = "white") +
	geom_text(aes(label  = round(r, 2),
								colour = if_else(flag, "red", "black")),
						size = 2.5) +
	scale_fill_distiller(palette = "RdBu", direction = 1, limits = c(-1, 1),
											 name = "r") +
	scale_colour_identity() +
	scale_x_discrete(guide = guide_axis(angle = 35)) +
	coord_fixed() +
	labs(x = NULL, y = NULL, caption = "Red = |r| > 0.5") +
	theme_succession(base_size = 8) +
	theme(
		panel.grid      = element_blank(),
		legend.position = "none"
	)

# ── Assemble compound figure ──────────────────────────────────────────────────
pca_block <- (scree_plot | var_contrib_pc1) /
	(var_contrib_pc2 | var_contrib_pc3)

validation_block <- p_biome | p_cors +
	plot_layout(heights = c(1.3,1))

s1 <- (pca_block / validation_block) +
	plot_layout(nrow = 3, heights = c(1, 1, 2)) +
	plot_annotation(
		tag_levels = "a",
		theme      = theme(plot.margin = margin(4, 4, 4, 4))
	)

ggsave(
	file.path(DIR_S1, "S1_pca_validation.png"),
	plot   = s1,
	bg     = "white",
	width  = 220,
	height = 300,
	units  = "mm",
	dpi    = 300
)
message("S-1 saved to figures/supplementary/pca/S1_pca_validation.png")