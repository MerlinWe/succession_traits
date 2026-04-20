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
PATH_PLOTS  <- "plots"
PATH_PCA    <- "data_processed/pca_rotation.rds"

export <- TRUE

# ── Helper functions ──────────────────────────────────────────────────────────
source("scripts/functions.R")

# ── Trait label mapping (used in downstream scripts; defined once here) ───────
trait_labels <- c(
	"bark_thickness"     = "Bark Thickness",
	"conduit_diam"       = "Conduit Diameter",
	"height"             = "Tree Height",
	"leaf_density"       = "Leaf Density",
	"leaf_k"             = "Leaf Potassium",
	"root_depth"         = "Root Depth",
	"seed_dry_mass"      = "Seed Dry Mass",
	"shade_tolerance"    = "Shade Tolerance",
	"specific_leaf_area" = "Specific Leaf Area"
)


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
# 6. Supplementary figure: scree plot + varimax contribution plots (S1)
# ══════════════════════════════════════════════════════════════════════════════

# Variable abbreviations for axis labels
abbreviations <- c(
	"annual_mean_temperature"            = "AnnMeanTemp",
	"max_temperature_of_warmest_month"   = "MaxTempWarmest",
	"min_temperature_of_coldest_month"   = "MinTempColdest",
	"mean_temperature_of_coldest_quarter"= "MeanTempColdest",
	"mean_temperature_of_warmest_quarter"= "MeanTempWarmest",
	"annual_precipitation"               = "AnnPrecip",
	"precipitation_of_driest_quarter"    = "PrecipDriest",
	"precipitation_of_wettest_quarter"   = "PrecipWettest",
	"precipitation_of_coldest_quarter"   = "PrecipColdest",
	"precipitation_of_warmest_quarter"   = "PrecipWarmest",
	"sand_content_015cm"                 = "Sand015",
	"sand_content_060cm"                 = "Sand060",
	"water_capacity_015cm"               = "WaterCap015",
	"water_capacity_060cm"               = "WaterCap060"
)

# Scree plot (using raw eigenvalues, not fviz wrapper, for ggplot control)
scree_plot <- tibble(
	component = paste0("PC", seq_along(eigenvalues)),
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
	labs(
		y = "Variance explained (%)",
		x = "Principal component",
		title = "Scree plot"
	) +
	theme_bw(base_size = 8)

# Varimax contribution plots (one per retained component)
varimax_contributions <- varimax_contrib(varimax_loadings)

var_contrib_pc1 <- tibble(
	name   = rownames(varimax_contributions),
	contrib = varimax_contributions[, 1]
) %>%
	create_contrib_plot(
		"PC1 (Temperature) — Varimax contributions",
		"firebrick4", abbreviations
	)

var_contrib_pc2 <- tibble(
	name   = rownames(varimax_contributions),
	contrib = varimax_contributions[, 2]
) %>%
	create_contrib_plot(
		"PC2 (Soil water retention) — Varimax contributions",
		"tan4", abbreviations
	)

var_contrib_pc3 <- tibble(
	name   = rownames(varimax_contributions),
	contrib = varimax_contributions[, 3]
) %>%
	create_contrib_plot(
		"PC3 (Precipitation) — Varimax contributions",
		"dodgerblue3", abbreviations
	)

pca_plot <- plot_grid(
	scree_plot, var_contrib_pc1, var_contrib_pc2, var_contrib_pc3,
	ncol = 2, nrow = 2,
	labels = "auto", label_fontfamily = "sans", label_size = 8
)

if (export) {
	ggsave(
		file.path(PATH_PLOTS, "/pca/s1_pca.png"),
		plot   = pca_plot,
		bg     = "white",
		width  = 200, height = 130, units = "mm", dpi = 600
	)
	message(sprintf("Supplementary figure saved to: %s/s1_pca.png", PATH_PLOTS))
}


# ══════════════════════════════════════════════════════════════════════════════
# 7. Attach rotated scores to data and finalise analysis dataset
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
# 8. Final diagnostics and output
# ══════════════════════════════════════════════════════════════════════════════

message("\n── Final dataset summary ──────────────────────────────────────────")
message(sprintf("  Plots:              %d", nrow(data_clean)))
message(sprintf("  Predictors:         standage + temp_pc + soil_pc + rain_pc + elevation + soil_ph"))
message(sprintf("  Traits:             %d", length(trait_labels)))
message(sprintf("  Biome dummies:      %d", sum(str_starts(names(data_clean), "biome_"))))
message(sprintf("  Stand age range:    %d – %d years",
								min(data_clean$standage), max(data_clean$standage)))

# Check rotated score distributions
message("\n── Rotated PC score summaries ─────────────────────────────────────")
data_clean %>%
	dplyr::select(temp_pc, soil_pc, rain_pc) %>%
	summary() %>%
	print()

write_rds(data_clean, PATH_OUT)
message(sprintf("\nOutput written to: %s", PATH_OUT))

# ══════════════════════════════════════════════════════════════════════════════
# 9. Sanity check plots
# ══════════════════════════════════════════════════════════════════════════════

sc_save <- function(p, name, w = 200, h = 140)
  ggsave(file.path("plots/sanity_checks", name), p,
         bg = "white", width = w, height = h, units = "mm", dpi = 300)

# Recover biome label from dummy columns (one per row by construction)
biome_cols <- names(data_clean) %>% str_subset("^biome_")
data_clean <- data_clean %>%
  mutate(biome_label = biome_cols %>%
    map_dfc(~ if_else(data_clean[[.x]] == 1L,
                      str_remove(.x, "^biome_") %>% str_replace_all("_", " ") %>% str_to_title(),
                      NA_character_)) %>%
    reduce(coalesce))

# PC scores by biome — checks rotation is ecologically interpretable
p1 <- data_clean %>%
  pivot_longer(c(temp_pc, soil_pc, rain_pc), names_to = "pc", values_to = "score") %>%
  mutate(pc = recode(pc, temp_pc = "Temperature PC",
                         soil_pc = "Soil PC", rain_pc = "Precipitation PC"),
         biome_label = fct_reorder(biome_label, score, median, na.rm = TRUE)) %>%
  ggplot(aes(biome_label, score, fill = pc)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3, linewidth = 0.3) +
  facet_wrap(~ pc, ncol = 1, scales = "free_x") +
  coord_flip() +
  scale_fill_manual(values = c("Temperature PC" = "#d73027",
                                "Soil PC" = "#8c510a",
                                "Precipitation PC" = "#4575b4")) +
  labs(x = NULL, y = "Rotated PC score",
       title = "PC scores by biome") +
  theme_bw(base_size = 9) + theme(legend.position = "none")

p1
sc_save(p1, "01_pc_by_biome.png", w = 200, h = 200)

# 2. Spatial maps — checks for projection errors
p2 <- data_clean %>%
  pivot_longer(c(temp_pc, soil_pc, rain_pc), names_to = "pc", values_to = "score") %>%
  mutate(pc = recode(pc, temp_pc = "Temperature PC",
                         soil_pc = "Soil PC", rain_pc = "Precipitation PC")) %>%
  ggplot(aes(LON, LAT, colour = score)) +
  geom_point(size = 0.3, alpha = 0.4) +
  facet_wrap(~ pc, ncol = 1) +
  scale_colour_distiller(palette = "RdBu", direction = 1, name = "Score",
                         limits = function(x) c(-max(abs(x)), max(abs(x)))) +
  coord_fixed(1.3) +
  labs(title = "PC scores spatially") +
  theme_bw(base_size = 9)

p2 
p2 <- sc_save(p2, "02_pc_spatial.png", w = 160, h = 220)

# 3. Predictor correlations — checks collinearity among model predictors
p3 <- c("standage","temp_pc","soil_pc","rain_pc","elevation","soil_ph") %>%
  { cor(data_clean[.], use = "pairwise.complete.obs") } %>%
  as_tibble(rownames = "v1") %>%
  pivot_longer(-v1, names_to = "v2", values_to = "r") %>%
  mutate(flag = abs(r) > 0.5 & v1 != v2) %>%
  ggplot(aes(v1, v2, fill = r)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = round(r, 2),
                colour = if_else(flag, "red", "black")), size = 2.8) +
  scale_fill_distiller(palette = "RdBu", direction = 1, limits = c(-1, 1)) +
  scale_colour_identity() +
  scale_x_discrete(guide = guide_axis(angle = 30)) +
  coord_fixed() +
  labs(x = NULL, y = NULL, title = "Predictor correlations — red > |0.5|") +
  theme_bw(base_size = 9) + theme(panel.grid = element_blank())

p3
p3 <- sc_save(p3, "03_predictor_cors.png", w = 150, h = 130)

message("Sanity check plots saved to plots/sanity_checks/")

