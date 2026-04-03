#################################################################
##### 01 FIA data preparation & trait plus clim integration #####
#################################################################

## Cleans and assembles plot-level community-weighted mean trait data from
## per-state FIA feather files, joins with trait and climate data, and writes
## a single analysis-ready RDS.

# ========================================
# SETUP 
# ========================================

rm(list = ls()) # clean environment 

library(stringi)
library(tidyverse)

# ---- Paths ----
PATH_FIA_DIR    <- "/Volumes/External/Data/fia_data/cleaned_state_data"
PATH_TRAITS     <- "data/Estimated_trait_table_with_monos.csv"
PATH_SHADE      <- "data/rueda_clean.csv"        
PATH_CLIM       <- "data/clim_composite.csv"
PATH_ECOREGIONS <- "data/resolve_ecoregions_legend.csv"
PATH_LATLON     <- "data/lat_lon_lookup.csv"

# ── Trait names used in the analysis ──────────────────────────────────────────

# Canonical short names as they appear in the raw trait table (before cleaning)
TRAIT_NAMES_RAW <- c(
	"Root depth", "Leaf density", "Specific leaf area",
	"Leaf K", "Conduit diam.", "Seed dry mass", "Bark thickness"
)

# ══════════════════════════════════════════════════════════════════════════════
# 1. Load trait data
# ══════════════════════════════════════════════════════════════════════════════

traits <- read_csv(PATH_TRAITS, show_col_types = FALSE) %>%
	
	# Keep only the seven focal traits, phylogenetically imputed values only
	filter(trait_short %in% TRAIT_NAMES_RAW, fit == "phy") %>%
	
	# Standardise trait column names: spaces → "_", remove ".", lowercase
	mutate(
		trait_short = trait_short %>%
			str_replace_all(" ", "_") %>%
			str_replace_all("\\.", "") %>%
			str_to_lower()
	) %>%
	
	dplyr::select(accepted_bin, trait_short, pred_value) %>%
	
	pivot_wider(names_from = trait_short, values_from = pred_value) %>%
	
	# Join shade tolerance index (Rueda et al. 2018)
	left_join(
		read_csv(PATH_SHADE, show_col_types = FALSE) %>%
			dplyr::select(accepted_bin, shade),
		by = "accepted_bin"
	) %>%
	rename(shade_tolerance = shade) %>%
	as_tibble()

# Diagnostic: how many species have shade tolerance data?
n_shade <- sum(!is.na(traits$shade_tolerance))
message(sprintf(
	"Trait table: %d species loaded; %d (%.1f%%) have shade tolerance values.",
	nrow(traits), n_shade, 100 * n_shade / nrow(traits)
))


# ══════════════════════════════════════════════════════════════════════════════
# 2. Climate composite
# ══════════════════════════════════════════════════════════════════════════════

clim <- read_csv(PATH_CLIM, show_col_types = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# 3. Per-state FIA preparation function
# ══════════════════════════════════════════════════════════════════════════════
# NOTE: trait log-transformation and scaling is applied once to the full assembled 
# dataset after all states are joined so that scaling parameters are derived from 
# the complete species pool and not the state-level subset.

prep_fia <- function(arrow_path, trait_data) {
	
	require(arrow)
	require(tidyverse)
	
	# ── 3.1 Read and filter raw tree-level data ────────────────────────────────
	fia <- arrow::read_feather(arrow_path) %>%
		
		# Living trees only
		mutate(alive = as.integer(alive)) %>%
		filter(alive == 1) %>%
		
		# Species-level taxonomy only (tax_level == 2)
		mutate(tax_level = as.integer(tax_level)) %>%
		filter(tax_level == 2) %>%
		
		# Single-condition plots only
		mutate(
			CONDID     = as.integer(CONDID),
			mult_conds = as.integer(mult_conds)
		) %>%
		filter(CONDID == 1, mult_conds == 0) %>%
		
		# Standard FIA plot designs (fixed-area; incl. CA/OR/WA macroplot codes)
		mutate(DESIGNCD = as.integer(DESIGNCD)) %>%
		filter(DESIGNCD %in% c(1L, 501L, 505L)) %>%
		
		# Public ownership only (exclude OWNGRPCD == 40 = private)
		mutate(OWNGRPCD = as.integer(OWNGRPCD)) %>%
		filter(OWNGRPCD != 40L) %>%
		
		# Diameter measured at breast height only
		mutate(DIAHTCD = as.integer(DIAHTCD)) %>%
		filter(DIAHTCD == 1L) %>%
		
		# Modern inventory period only
		filter(INVYR >= 1980L) %>%
		
		dplyr::select(
			accepted_bin, TREEID, PID, SUBP, PLOT_TYPE,
			TPA_UNADJ, INVYR, FORTYPCD, STDAGE,
			DIA, CARBON_AG, biomass, height,
			ownership, biome, managed, ll_id, min_year, max_year
		) %>%
		
		# ── 3.2 Tree-level derived variables ──────────────────────────────────────
		
		# Basal area (cm²) from diameter at breast height
		mutate(BA = pi * (DIA / 2)^2) %>%
		
		# Scale stem-level quantities by trees-per-acre expansion factor,
		# converting individual measurements to per-acre equivalents
		mutate(across(
			c(CARBON_AG, BA, biomass),
			~ . * TPA_UNADJ,
			.names = "{col}_scaled"
		)) %>%
		dplyr::select(-TPA_UNADJ, -CARBON_AG, -BA, -biomass) %>%
		
		# Log-transform height and diameter for tree-level aggregation.
		# NOTE: these are NOT z-scored here. Scaling is deferred to the
		# full-dataset trait standardisation step (Section 5).
		mutate(
			DIA    = log(DIA),
			height = log(height)
		) %>%
		
		# Broad forest type from FIA forest-type code
		mutate(foresttype = case_when(
			FORTYPCD >= 101L & FORTYPCD <= 199L ~ "Coniferous",
			FORTYPCD >= 201L & FORTYPCD <= 299L ~ "Deciduous",
			FORTYPCD >= 301L & FORTYPCD <= 399L ~ "Mixed",
			TRUE                                ~ "Non-stocked/Other"
		)) %>%
		
		# Unique identifier combining plot × inventory year for repeated measures
		mutate(PID_rep = paste0(PID, "-", INVYR))
	
	# ── 3.3 Aggregate to species × plot level ─────────────────────────────────
	
	fia <- fia %>%
		group_by(
			accepted_bin, PID_rep,
			STDAGE, INVYR, FORTYPCD, foresttype,
			biome, ownership, managed, ll_id
		) %>%
		summarise(
			# BA-weighted mean diameter and height at species × plot level
			wmean_DIA_species    = weighted.mean(DIA,    BA_scaled, na.rm = TRUE),
			wmean_height_species = weighted.mean(height, BA_scaled, na.rm = TRUE),
			# Per-acre totals
			sum_BA_species       = sum(BA_scaled,       na.rm = TRUE),
			sum_biomass_species  = sum(biomass_scaled,   na.rm = TRUE),
			sum_CARBON_AG_species = sum(CARBON_AG_scaled, na.rm = TRUE),
			.groups = "drop"
		) %>%
		as_tibble()
	
	# ── 3.4 Join with trait data ───────────────────────────────────────────────
	# Traits are joined here in raw (untransformed) form. Log-transform and
	# z-scoring happen once on the full assembled dataset (Section 5).
	
	dat <- trait_data %>%
		filter(accepted_bin %in% fia$accepted_bin) %>%
		right_join(fia, by = "accepted_bin")
	
	# ── 3.5 Aggregate to plot level (community-weighted means) ────────────────
	
	fia_plot <- dat %>%
		group_by(
			PID_rep,
			standage = STDAGE,   # rename here for clarity in output
			INVYR, FORTYPCD, foresttype,
			biome, ownership, managed, ll_id
		) %>%
		summarise(
			
			# Functional traits (BA-weighted CWMs; raw scale — transformed later)
			wmean_bark_thickness    = weighted.mean(bark_thickness,    sum_BA_species, na.rm = TRUE),
			wmean_conduit_diam      = weighted.mean(conduit_diam,      sum_BA_species, na.rm = TRUE),
			wmean_leaf_density      = weighted.mean(leaf_density,      sum_BA_species, na.rm = TRUE),
			wmean_leaf_k            = weighted.mean(leaf_k,            sum_BA_species, na.rm = TRUE),
			wmean_root_depth        = weighted.mean(root_depth,        sum_BA_species, na.rm = TRUE),
			wmean_seed_dry_mass     = weighted.mean(seed_dry_mass,     sum_BA_species, na.rm = TRUE),
			wmean_specific_leaf_area = weighted.mean(specific_leaf_area, sum_BA_species, na.rm = TRUE),
			wmean_shade_tolerance   = weighted.mean(shade_tolerance,   sum_BA_species, na.rm = TRUE),
			
			# Tree size (already log-transformed at tree level; BA-weighted)
			wmean_dia    = weighted.mean(wmean_DIA_species,    sum_BA_species, na.rm = TRUE),
			wmean_height = weighted.mean(wmean_height_species, sum_BA_species, na.rm = TRUE),
			
			# Shade tolerance data coverage: fraction of BA represented
			shade_ba_coverage = sum(sum_BA_species[!is.na(shade_tolerance)], na.rm = TRUE) /
				sum(sum_BA_species, na.rm = TRUE),
			
			# Plot-level totals
			total_carbon_ag  = sum(sum_CARBON_AG_species, na.rm = TRUE),
			total_biomass    = sum(sum_biomass_species,   na.rm = TRUE),
			total_basal_area = sum(sum_BA_species,        na.rm = TRUE),
			
			.groups = "drop"
		) %>%
		# Replace any NaN values (from weighted.mean on all-NA inputs) with NA
		mutate(across(where(is.numeric), ~ if_else(is.nan(.), NA_real_, .))) %>%
		as_tibble()
	
	# ── 3.6 Repeated-measures bookkeeping ────────────────────────────────────
	# Reconstruct the original PID, flag repeated plots, and keep only the
	# three most recent measurements per plot to limit temporal autocorrelation.
	
	fia_plot <- fia_plot %>%
		mutate(PID = str_remove(PID_rep, "-\\d{4}$")) %>%
		group_by(PID) %>%
		mutate(rep_measure = n() > 1L) %>%
		# slice_max is explicit and unambiguous vs arrange + row_number
		slice_max(order_by = INVYR, n = 3L, with_ties = FALSE) %>%
		mutate(PID_measure = row_number()) %>%
		ungroup()
	
	return(fia_plot)
}

# ══════════════════════════════════════════════════════════════════════════════
# 4. Run pipeline across all state files
# ══════════════════════════════════════════════════════════════════════════════

prep_fia_safe <- safely(prep_fia, otherwise = NULL)

fia_data <- tibble(
	pathway = list.files(PATH_FIA_DIR, full.names = TRUE)
) %>%
	mutate(state = str_extract(pathway, "(?<=/FIA_)[^\\.]+"))

message(sprintf("Processing %d state file(s)...", nrow(fia_data)))

fia_data <- fia_data %>%
	mutate(
		result       = map(pathway, ~ prep_fia_safe(.x, traits)),
		data         = map(result, "result"),
		error        = map(result, "error"),
		has_problems = map_lgl(error, ~ !is.null(.x))
	)

# Report any failures before discarding them
failed_states <- fia_data %>% filter(has_problems)
if (nrow(failed_states) > 0) {
	warning(sprintf("%d state(s) failed and will be excluded:", nrow(failed_states)))
	walk2(failed_states$state, failed_states$error, ~ {
		message(sprintf("  ✗ %s: %s", .x, conditionMessage(.y)))
	})
} else {
	message("  ✓ All states processed successfully.")
}

# Assemble successful states only
fia_clean <- fia_data %>%
	filter(!has_problems) %>%
	dplyr::select(state, data) %>%
	unnest(data)

message(sprintf(
	"Assembled dataset: %d plots across %d state(s).",
	nrow(fia_clean),
	n_distinct(fia_clean$state)
))

# ══════════════════════════════════════════════════════════════════════════════
# 5. Trait transformation and scaling (full dataset, applied once)
# ══════════════════════════════════════════════════════════════════════════════

# All traits are log-transformed then z-scored using the grand mean and SD across the complete assembled dataset. 
# Tree size variables (wmean_dia, wmean_height) were log-transformed at the tree level before aggregation 
# and are z-scored here. shade_ba_coverage is a proportion and is not transformed.

TRAIT_COLS <- c(
	"wmean_bark_thickness", "wmean_conduit_diam", "wmean_leaf_density",
	"wmean_leaf_k",         "wmean_root_depth",   "wmean_seed_dry_mass",
	"wmean_specific_leaf_area", "wmean_shade_tolerance"
)

SIZE_COLS <- c("wmean_dia", "wmean_height")

fia_clean <- fia_clean %>%
	# Log-transform raw trait CWMs (all positive by construction after imputation)
	mutate(across(all_of(TRAIT_COLS), log)) %>%
	# Z-score all trait and size columns using full-dataset parameters
	mutate(across(all_of(c(TRAIT_COLS, SIZE_COLS)), ~ as.numeric(scale(.)))) %>%
	# Diagnostic: flag any plots where shade tolerance is poorly represented
	mutate(shade_low_coverage = shade_ba_coverage < 0.5)

# Report shade tolerance coverage
n_low <- sum(fia_clean$shade_low_coverage, na.rm = TRUE)
message(sprintf(
	"Shade tolerance coverage: %d plots (%.1f%%) have <50%% BA covered.",
	n_low, 100 * n_low / nrow(fia_clean)
))

# ══════════════════════════════════════════════════════════════════════════════
# 6. Join climate, ecoregion, and coordinate data
# ══════════════════════════════════════════════════════════════════════════════

ecoregion_labels <- read_csv(
	PATH_ECOREGIONS,
	locale = locale(encoding = "UTF-8"),
	show_col_types = FALSE
) %>%
	mutate(
		Name = iconv(Name, from = "UTF-8", to = "ASCII//TRANSLIT"),
		Name = stri_replace_all_regex(str_to_lower(Name), "\\s+", "_"),
		Code = as.integer(Code)
	)

latlon <- read_csv(PATH_LATLON, show_col_types = FALSE) %>%
	distinct() %>%
	mutate(ll_id = as.double(ll_id))

fia_clean <- fia_clean %>%
	mutate(ll_id = as.double(ll_id)) %>%
	left_join(clim, by = "ll_id") %>%
	left_join(ecoregion_labels, by = c("ecoregion" = "Code")) %>%
	mutate(ecoregion = Name) %>%
	dplyr::select(-Name) %>%
	left_join(latlon, by = "ll_id")


# ══════════════════════════════════════════════════════════════════════════════
# 7. Final column selection, type coercion, and integrity checks
# ══════════════════════════════════════════════════════════════════════════════

fia_clean <- fia_clean %>%
	dplyr::select(
		# Plot identifiers and metadata
		PID, PID_rep, rep_measure, PID_measure,
		state, standage, INVYR, FORTYPCD, foresttype,
		biome, ownership, managed, ll_id, LAT, LON, ecoregion,
		
		# Community-weighted mean traits (log-transformed, z-scored)
		wmean_bark_thickness, wmean_conduit_diam,   wmean_leaf_density,
		wmean_leaf_k,         wmean_root_depth,     wmean_seed_dry_mass,
		wmean_specific_leaf_area, wmean_shade_tolerance,
		wmean_dia,            wmean_height,
		
		# Shade tolerance data quality flag
		shade_ba_coverage, shade_low_coverage,
		
		# Plot-level ecosystem variables
		total_carbon_ag, total_biomass, total_basal_area,
		
		# Climate variables (raw; PCA is run in 02_environment_pca.R)
		annual_mean_temperature, annual_precipitation,
		isothermality, mean_diurnal_range, temperature_seasonality,
		temperature_annual_range,
		max_temperature_of_warmest_month, min_temperature_of_coldest_month,
		mean_temperature_of_warmest_quarter, mean_temperature_of_coldest_quarter,
		mean_temperature_of_driest_quarter, mean_temperature_of_wettest_quarter,
		precipitation_seasonality,
		precipitation_of_wettest_month, precipitation_of_driest_month,
		precipitation_of_wettest_quarter, precipitation_of_driest_quarter,
		precipitation_of_warmest_quarter, precipitation_of_coldest_quarter,
		elevation, pop_density,
		starts_with("sand_content_"),
		starts_with("soil_ph_"),
		starts_with("water_capacity_")
	) %>%
	
	# Type coercion
	mutate(
		PID           = as.character(PID),
		PID_rep       = as.character(PID_rep),
		rep_measure   = as.logical(rep_measure),
		PID_measure   = as.integer(PID_measure),
		state         = as.factor(state),
		standage      = as.integer(standage),
		INVYR         = as.integer(INVYR),
		FORTYPCD      = as.integer(FORTYPCD),
		foresttype    = as.factor(foresttype),
		biome         = as.factor(biome),
		ownership     = as.factor(ownership),
		managed       = as.factor(managed),
		ll_id         = as.character(ll_id),
		ecoregion     = as.character(ecoregion),
		across(where(is.numeric), as.numeric)
	) %>%
	
	# Drop rows with missing stand age — essential predictor
	filter(!is.na(standage))

# ══════════════════════════════════════════════════════════════════════════════
# 8. Analysis-ready filtering and biome encoding
# ══════════════════════════════════════════════════════════════════════════════
# These filters define the analytical sample and are applied here so that
# the output file exactly represents the data used in all downstream scripts.

# Exclude upper 10% of stand age (focus on early-to-late succession, 
# avoids influence of extremely old-growth plots)
standage_cap <- quantile(fia_clean$standage, 0.9, na.rm = TRUE)
message(sprintf("Standage cap (90th percentile): %d years", standage_cap))

fia_clean <- fia_clean %>%
	filter(standage < standage_cap) %>%
	# Exclude actively managed plots
	filter(managed == "0" | managed == 0)

# Keep only biomes with sufficient representation (>100 plots)
biome_counts <- count(fia_clean, biome)
biomes_keep  <- biome_counts %>% filter(n > 100) %>% pull(biome)

message(sprintf(
	"Biomes retained (n > 100 plots): %d of %d",
	length(biomes_keep), nrow(biome_counts)
))
message("Dropped biomes:")
biome_counts %>% filter(!biome %in% biomes_keep) %>%
	pwalk(~ message(sprintf("  ✗ %s: %d plots", ..1, ..2)))

fia_clean <- fia_clean %>%
	filter(biome %in% biomes_keep)

# Create biome dummy variables (one per retained biome, used in RF models)
biome_dummies <- model.matrix(~ biome - 1, data = fia_clean) %>%
	as_tibble() %>%
	rename_with(~ .x %>%
								str_replace_all("biome(?!_)", "biome_") %>%
								str_replace_all(" ", "_") %>%
								str_to_lower()
	)

fia_clean <- bind_cols(fia_clean, biome_dummies) %>%
	dplyr::select(-biome)

message(sprintf("Final analytical sample: %d plots", nrow(fia_clean)))


# ── Integrity checks ──────────────────────────────────────────────────────────

stopifnot(
	"Duplicate PID_rep rows detected" = !any(duplicated(fia_clean$PID_rep)),
	"PID_measure exceeds 3"           = max(fia_clean$PID_measure, na.rm = TRUE) <= 3L,
	"Negative basal area values"      = all(fia_clean$total_basal_area >= 0, na.rm = TRUE)
)

# Summary report
message("\n── Final dataset summary ──────────────────────────────────────────")
message(sprintf("  Plots (PID_rep):    %d", nrow(fia_clean)))
message(sprintf("  Unique plots (PID): %d", n_distinct(fia_clean$PID)))
message(sprintf("  States:             %d", n_distinct(fia_clean$state)))
message(sprintf("  Biomes:             %d", n_distinct(fia_clean$biome)))
message(sprintf("  Stand age range:    %d – %d years",
								min(fia_clean$standage), max(fia_clean$standage)))
message(sprintf("  Repeated measures:  %d plots (%.1f%%)",
								sum(fia_clean$rep_measure),
								100 * mean(fia_clean$rep_measure)))

# Trait completeness table
trait_completeness <- fia_clean %>%
	dplyr::select(all_of(TRAIT_COLS)) %>%
	summarise(across(everything(), ~ mean(!is.na(.)))) %>%
	pivot_longer(everything(), names_to = "trait", values_to = "completeness") %>%
	mutate(completeness = scales::percent(completeness, accuracy = 0.1))

message("\n── Trait completeness ─────────────────────────────────────────────")
print(trait_completeness, n = Inf)
message("───────────────────────────────────────────────────────────────────\n")

# ══════════════════════════════════════════════════════════════════════════════
# 8. Write output
# ══════════════════════════════════════════════════════════════════════════════

out_path <- paste0("data_processed/plotlevel_data.rds")
write_rds(fia_clean, file = out_path)
message(sprintf("Output written to: %s", out_path))

