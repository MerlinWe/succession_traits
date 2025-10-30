##### MSc Diss. FIA data prep w/ trait and climate composites #####
## Author: M. Weiss @ Maynard Lab UCL 


rm(list = ls()) # clean environment 

library(DT)
library(arrow)
library(stringi)
library(cowplot)
library(tidyverse)

## Get trait database
traits <- read_csv("data/Estimated_trait_table_with_monos.csv") %>%
	
	# Keep only relevant traits to the analysis 
	filter(trait_short %in% c("Root depth", "Leaf density", "Specific leaf area", "Leaf K", 
														"Conduit diam.", "Seed dry mass", "Bark thickness")) %>%
	
	# Filter by imputation: only phylogenetically imputed values 
	filter(fit == "phy") %>% 
	
	# Clean trait column (replace white spaces)
	mutate(trait_short = gsub(" ", "_", trait_short),
				 trait_short = gsub("\\.", "", trait_short),
				 trait_short = tolower(trait_short)) %>%
	
	# Keep only species and (predicted) trait
	dplyr::select(accepted_bin, trait_short, pred_value) %>%
	
	# Pivot wider
	pivot_wider(names_from = trait_short, 
							values_from = pred_value,
							names_sep = "_") %>%
	# Integrate symptoms data (from Rueda et al.)
	left_join(
		read_csv("/Users/merlin/Documents/MSc/Thesis/Data/traits/Rueda_etal/rueda_clean.csv") %>%
			dplyr::select(accepted_bin, shade), by = "accepted_bin") %>%
	dplyr::rename(shade_tolerance = shade) %>%
	
	# Set format 
	as_tibble() 

## Get climate database
clim <- read_csv("data/clim_composite.csv")

prep_fia <- function(arrow_path, trait_data){
	
	# Make sure libraries are open
	require(arrow)
	require(tidyverse)
	
	########## Run preparation pipeline on state FIA data ##########
	
	# Read state FIA file using pathname 
	fia <- arrow::read_feather(arrow_path) %>%
		
		## Filtering (while making sure formats are correct) 
		
		# Dead: No dead trees 
		mutate(alive = as.integer(alive)) %>%
		filter(alive == 1) %>%
		
		# Taxonomy: Only species level trees
		mutate(tax_level = as.integer(tax_level)) %>%
		filter(tax_level == 2) %>%
		
		# Condition: Only plots with condition 1
		mutate(CONDID = as.integer(CONDID)) %>%
		filter(CONDID == 1) %>%
		mutate(mult_conds = as.integer(mult_conds)) %>%
		filter(mult_conds == 0) %>%
		
		# Design: Only fixed design (plus with macros to add CA, OR, WA)
		mutate(DESIGNCD = as.integer(DESIGNCD)) %>%
		filter(DESIGNCD %in% c(505, 501, 1)) %>%
		
		# Ownership: Disregard private ownerships
		mutate(OWNGRPCD = as.integer(OWNGRPCD)) %>%
		filter(OWNGRPCD != 40 ) %>%
		
		# DIA: Only diameter values measured at breast height
		mutate(DIAHTCD = as.integer(DIAHTCD)) %>%
		filter(DIAHTCD == 1) %>%
		
		# Keep only inventory years after 1980
		filter(INVYR >= 1980) %>%
		
		# Remove columns with no variation and set format 
		dplyr::select(accepted_bin, TREEID, PID, SUBP, PLOT_TYPE, TPA_UNADJ, INVYR, FORTYPCD, 
									STDAGE, DIA, CARBON_AG, biomass, height,
									ownership, biome, managed, ll_id, min_year, max_year) %>%
		
		## Calculate Basal Area using the DIA 
		mutate(BA = pi * (DIA/2)^2) %>%
		
		## Scale numeric tree level variables by TPA
		mutate(across(
			c(CARBON_AG, BA, biomass), ~ . * TPA_UNADJ, .names = "{col}_scaled")) %>%
		dplyr::select(-TPA_UNADJ, -CARBON_AG, -BA, -biomass) %>%	
		
		# Log scale height and diameter
		mutate(DIA = as.numeric(scale(log(DIA)))) %>%
		mutate(height = as.numeric(scale(log(height)))) %>%
		
		# Classify forest type based on the forest type codes 
		mutate(foresttype = case_when(
			FORTYPCD >= 101 & FORTYPCD <= 199 ~ "Coniferous",
			FORTYPCD >= 201 & FORTYPCD <= 299 ~ "Deciduous",
			FORTYPCD >= 301 & FORTYPCD <= 399 ~ "Mixed",
			TRUE ~ "Non-stocked/Other"
		)) %>%
		
		# Add a unique ID for repeated plot measures
		mutate(PID_rep = paste0(PID, "-", INVYR))
	
	## Summarize to species per plot obs level
	
	# First do a safety check to make sure there are no more than 1 levels to PID per 
	# grouping variable (excluding accepted_bin): 
	
	grouping_vars <- c("STDAGE", "INVYR", "FORTYPCD", "foresttype", "biome","ownership", "managed", "ll_id")
	
	for (var in grouping_vars) {
		level_check <- fia %>%
			group_by(PID_rep) %>%
			summarise(levels_count = n_distinct(.data[[var]])) %>%
			filter(levels_count > 1)
		
		if (nrow(level_check) > 0) {
			stop(paste("Error: More than one level per PID in grouping variable", var, "detected!"))
		}
	}
	
	# Now summarize data to species per plot level
	fia <- fia %>%
		
		# Group by species and plot level variables
		group_by(accepted_bin, PID_rep, STDAGE, INVYR, FORTYPCD, foresttype, biome, ownership, managed, ll_id) %>%
		
		# Summarise to species per plot level 
		summarize(
			
			# Weighted means for DIA and height weighted by BA_scaled 
			wmean_DIA_species = weighted.mean(DIA, BA_scaled, na.rm = TRUE),
			wmean_height_species = weighted.mean(height, BA_scaled, na.rm = TRUE),
			
			# Sums for the rest
			sum_BA_species = sum(BA_scaled, na.rm = TRUE),
			sum_biomass_species = sum(biomass_scaled, na.rm = TRUE),
			sum_CARBON_AG_species = sum(CARBON_AG_scaled, na.rm = TRUE),
			
			.groups = 'drop') %>%
		
		## Finish and set format
		ungroup() %>%
		as_tibble()
	
	########## Join FIA data with trait data ##########
	
	## Either use raw data: 
	
	dat <- trait_data %>%
		
		# Keep only species found in FIA data 
		filter(accepted_bin %in% fia$accepted_bin) %>%
		
		## LOG-TRANSFORM and SCALE trait data if needed 
		mutate(across(where(is.double), log)) %>%
		mutate(across(where(is.double), scale)) %>%
		mutate(across(where(is.double), as.numeric)) %>%
		
		# Join traits with FIA data 
		right_join(fia, by = "accepted_bin")
	
	########## Summarise functional diversity and build database ##########
	
	# Summarise and scale trait data and format as matrix to calculate indices
	fun_traits <- dat %>%
		
		# Keep relevant columns 
		dplyr::select(accepted_bin, root_depth, bark_thickness, seed_dry_mass, leaf_k, 
									leaf_density, conduit_diam, specific_leaf_area, shade_tolerance,
									wmean_DIA_species, wmean_height_species, sum_BA_species) %>%
		
		# Group by species
		group_by(accepted_bin) %>%
		
		# Summarise traits to plot level 
		summarise(
			
			# Tree Means weighted by BA
			dia = mean(wmean_DIA_species, na.rm = TRUE),
			height = mean(wmean_height_species, na.rm = TRUE),
			
			# Species Means; equal to trait data
			root_depth = mean(root_depth, na.rm = TRUE),
			bark_thickness = mean(bark_thickness, na.rm = TRUE),
			seed_dry_mass = mean(seed_dry_mass, na.rm = TRUE),
			leaf_k = mean(leaf_k, na.rm = TRUE),
			leaf_density = mean(leaf_density, na.rm = TRUE),
			conduit_diam = mean(conduit_diam, na.rm = TRUE),
			specific_leaf_area = mean(specific_leaf_area, na.rm = TRUE),
			shade_tolerance = mean(shade_tolerance, na.rm = TRUE),
			.groups = 'drop') %>%
		
		# Drop missing values 
		filter(complete.cases(.)) %>%
		
		# Format as matrix 
		column_to_rownames(var = "accepted_bin") %>%
		as.matrix() %>%
		
		# Scale values
		scale()
	
	# Summarise plot data and format as presence/absence matrix 
	fun_plots <- dat %>%
		# Keep only relevant columns 
		dplyr::select(accepted_bin, PID_rep) %>%
		# Set presence to 1
		mutate(present = 1) %>%
		# Keep only species found in trait data 
		filter(accepted_bin %in% rownames(fun_traits)) %>%
		# Pivot to wider 
		pivot_wider(names_from = accepted_bin, values_from = present, values_fill = list(present = 0)) %>%
		# Format as matrix 
		column_to_rownames(var = "PID_rep") %>%
		as.matrix()
	
	## Summarise to plot level, calculate resource use score, and functional indices
	
	# Summarise data to plot level 
	fia <- dat %>%
		
		# Define grouping variables that are constant for every plot 
		group_by(PID_rep, STDAGE, INVYR, FORTYPCD, foresttype, biome, ownership, managed, ll_id) %>%
		rename(standage = STDAGE) %>%
		
		# Summarise to plot levels by weighted means (using BA)
		summarise(
			
			# Root depth
			wmean_root_depth = weighted.mean(root_depth, sum_BA_species, na.rm = TRUE),
			
			# Wood density 
			wmean_bark_thickness = weighted.mean(bark_thickness, sum_BA_species, na.rm = TRUE),
			
			# Seed dry mass
			wmean_seed_dry_mass = weighted.mean(seed_dry_mass, sum_BA_species, na.rm = TRUE),
			
			# Leaf potassium 
			wmean_leaf_k = weighted.mean(leaf_k, sum_BA_species, na.rm = TRUE),
			
			# Leaf density 
			wmean_leaf_density = weighted.mean(leaf_density, sum_BA_species, na.rm = TRUE),

			# Conduit diameter
			wmean_conduit_diam = weighted.mean(conduit_diam, sum_BA_species, na.rm = TRUE),
			
			# Specific leaf area
			wmean_specific_leaf_area = weighted.mean(specific_leaf_area, sum_BA_species, na.rm = TRUE),
			
			# Shade tolerance 
			wmean_shade_tolerance = weighted.mean(shade_tolerance, sum_BA_species, na.rm = TRUE),
			
			# Diameter
			wmean_dia = weighted.mean(wmean_DIA_species, sum_BA_species, na.rm = TRUE),
			
			# Height 
			wmean_height = weighted.mean(wmean_height_species, sum_BA_species, na.rm = TRUE),
			
			# Above the ground carbon (SUM)
			total_carbon_ag = sum(sum_CARBON_AG_species),
			
			# Biomass (SUM)
			total_biomass = sum(sum_biomass_species),
			
			# Basal area (SUM)
			total_basal_area = sum(sum_BA_species),
			.groups = 'drop') %>%
		
		ungroup() %>%
		
		# Set format
		mutate_all(~ ifelse(is.nan(.), NA, .)) %>%
		as_tibble() %>%
		
		## Assess repeated measures: Check if TRUE and rank measurements, keep only most recent 3  
		
		# Reconstruct original PID variable to use for grouping
		mutate(PID = str_remove(PID_rep, "-\\d{4}$")) %>%
		group_by(PID) %>%
		
		# Check if PID is measured more than once
		mutate(rep_measure = n() > 1 ) %>%
		# Sort by PID and descending INVYR
		arrange(PID, desc(INVYR)) %>%  
		# Rank measurements based on sorted order
		mutate(PID_measure = if_else(rep_measure, row_number(), 1L)) %>%  
		# Keep only the most recent 3 measurements
		filter(PID_measure <= 3) %>%  
		ungroup()
	
	## Return cleaned FIA data 
	return(fia)
}


# Wrap the FIA prep pipeline with safely
prep_fia_safe <- safely(prep_fia, otherwise = NA)

# Access FIA pathways 
fia_data <- tibble(pathway = list.files(path = "/Volumes/External/MSc Thesis/fia_data/cleaned_state_data",
																				full.names = TRUE)) %>%
	# Extract state name from file name 
	mutate(state = str_extract(pathway, "(?<=/FIA_)[^\\.]+"))

# Pass data through preparation pipeline while skipping problematic files 
fia_data <- fia_data %>%
	mutate(result = map(pathway, ~ prep_fia_safe(.x, traits)),
				 data = map(result, "result"),
				 has_problems = map_lgl(result, ~ !is.null(.x$error))) %>%
	select(state, data, has_problems)

# Get available data, merge with climate data, and write file 
fia_clean <- fia_data %>%
	unnest(data) %>%
	
	# Final variable selection and order 
	dplyr::select(PID, PID_rep, rep_measure, PID_measure, state, standage, INVYR, FORTYPCD, foresttype, biome, ownership, managed, ll_id,
								
								wmean_root_depth, wmean_bark_thickness, wmean_conduit_diam, wmean_leaf_k, wmean_specific_leaf_area,
								wmean_leaf_density, wmean_seed_dry_mass, wmean_shade_tolerance, wmean_dia, wmean_height, 
								
								total_carbon_ag, total_biomass, total_basal_area)

fia_clean <- fia_clean %>%
	
	# Integrate climate data 
	left_join(clim, by = "ll_id") %>%
	
	# Add Ecoregion labels
	left_join(
		read_csv("data/resolve_ecoregions_legend.csv", locale = locale(encoding = "UTF-8")) %>%
			mutate(
				Name = iconv(Name, from = "UTF-8", to = "ASCII//TRANSLIT"), 
				Name = stri_replace_all_regex(tolower(Name), "\\s", "_")) %>%
			mutate(Code = as.integer(Code)),
		by = c("ecoregion" = "Code")) %>%
	mutate(ecoregion = Name) %>%
	select(-Name) %>%
	# Add full coordinates
	mutate(ll_id = as.double(ll_id)) %>%
	left_join(
		distinct(
			read_csv("data/lat_lon_lookup.csv")
			),
		by = "ll_id") %>%
	as_tibble()

# Set formats 
fia_clean <- fia_clean %>%
	mutate(
		PID = as.character(PID),
		PID_rep = as.character(PID_rep),
		rep_measure = as.logical(rep_measure),
		PID_measure = as.integer(PID_measure),
		state = as.factor(state),
		standage = as.integer(standage),
		INVYR = as.integer(INVYR),
		FORTYPCD = as.integer(FORTYPCD),
		foresttype = as.factor(foresttype),
		biome = as.factor(biome),
		ownership = as.factor(ownership),
		managed = as.factor(managed),
		ll_id = as.character(ll_id),
		across(starts_with("wmean_"), as.numeric),
		across(starts_with("total_"), as.numeric),
		annual_mean_temperature = as.numeric(annual_mean_temperature),
		annual_precipitation = as.numeric(annual_precipitation),
		isothermality = as.numeric(isothermality),
		max_temperature_of_warmest_month = as.numeric(max_temperature_of_warmest_month),
		mean_diurnal_range = as.numeric(mean_diurnal_range),
		mean_temperature_of_coldest_quarter = as.numeric(mean_temperature_of_coldest_quarter),
		mean_temperature_of_driest_quarter = as.numeric(mean_temperature_of_driest_quarter),
		mean_temperature_of_warmest_quarter = as.numeric(mean_temperature_of_warmest_quarter),
		mean_temperature_of_wettest_quarter = as.numeric(mean_temperature_of_wettest_quarter),
		min_temperature_of_coldest_month = as.numeric(min_temperature_of_coldest_month),
		precipitation_seasonality = as.numeric(precipitation_seasonality),
		precipitation_of_coldest_quarter = as.numeric(precipitation_of_coldest_quarter),
		precipitation_of_driest_month = as.numeric(precipitation_of_driest_month),
		precipitation_of_driest_quarter = as.numeric(precipitation_of_driest_quarter),
		precipitation_of_warmest_quarter = as.numeric(precipitation_of_warmest_quarter),
		precipitation_of_wettest_month = as.numeric(precipitation_of_wettest_month),
		precipitation_of_wettest_quarter = as.numeric(precipitation_of_wettest_quarter),
		temperature_annual_range = as.numeric(temperature_annual_range),
		temperature_seasonality = as.numeric(temperature_seasonality),
		elevation = as.numeric(elevation),
		pop_density = as.numeric(pop_density),
		ecoregion = as.character(ecoregion),
		across(starts_with("sand_content_"), as.numeric),
		across(starts_with("soil_ph_"), as.numeric),
		across(starts_with("water_capacity_"), as.numeric),
		LAT = as.numeric(LAT),
		LON = as.numeric(LON)) %>%
	
	# Drop missing standage rows
	filter(complete.cases(standage))

fia_clean %>% write_csv(file = paste0("data/plotlevel_data_", Sys.Date(),".csv"))
