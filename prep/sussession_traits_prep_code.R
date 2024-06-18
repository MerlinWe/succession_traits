##### MSc Diss. FIA data prep w/ trait and climate composites #####
## Author: M. Weiss @ Maynard Lab UCL 

rm(list = ls()) # clean environment 

library(DT)
library(arrow)
library(fundiversity)
library(cowplot)
library(tidyverse)

rd_pathway <- "/Volumes/ritd-ag-project-rd01pr-dmayn10" # path to UCL research drive

## Get trait database
traits <- read_csv("/Volumes/ritd-ag-project-rd01pr-dmayn10/trait_data/Estimated_trait_table_with_monos.csv") %>%
	
	# Keep only relevant traits to the analysis 
	filter(trait_short %in% c("Wood density", "Bark thickness", "Seed dry mass", 
														"Leaf N", "Specific leaf area", "Conduit diam.")) %>%
	
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
		read_csv("/Volumes/ritd-ag-project-rd01pr-dmayn10/merlin/data/traits/rueda_clean.csv") %>%
			dplyr::select(accepted_bin, shade), by = "accepted_bin") %>%
	dplyr::rename(shade_tolerance = shade) %>%
	
	# Set format 
	as_tibble() 

## Get climate database
clim <- read_csv("/Volumes/ritd-ag-project-rd01pr-dmayn10/fia_data/composite/composite_FIA_alastair.csv") %>%
	
	# Keep only CHELSA BIO climate variables and filter to FIA
	dplyr::select(ll_id, starts_with("CHELSA_BIO")) %>%
	rename_with(~ str_replace_all(., "CHELSA_BIO_", "") %>% tolower()) %>%
	
	# Set formats 
	mutate(annual_mean_temperature = annual_mean_temperature/10) %>%
	mutate(isothermality = isothermality/10) %>%
	mutate(max_temperature_of_warmest_month = max_temperature_of_warmest_month/10) %>%
	mutate(mean_temperature_of_coldest_quarter = mean_temperature_of_coldest_quarter/10) %>%
	mutate(mean_temperature_of_driest_quarter = mean_temperature_of_driest_quarter/10) %>%
	mutate(mean_temperature_of_warmest_quarter = mean_temperature_of_warmest_quarter/10) %>%
	mutate(mean_temperature_of_wettest_quarter = mean_temperature_of_wettest_quarter/10) %>%
	mutate(min_temperature_of_coldest_month = min_temperature_of_coldest_month/10) %>%
	mutate(temperature_annual_range = temperature_annual_range/10) %>%
	mutate(temperature_seasonality = temperature_seasonality/10) %>%
	
	# Write file 
	write_csv(file = "/Volumes/ritd-ag-project-rd01pr-dmayn10/merlin/data/climate/eth_clim_composite.csv")


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
		
		# Keep only the most recent observation per PID
		group_by(PID) %>%
		filter(INVYR == max(INVYR)) %>%
		ungroup() %>%
		
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
		))
	
	## Summarize to species per plot level
	
	# First do a safety check to make sure there are no more than 1 levels to PID per 
	# grouping variable (excluding accepted_bin): 
	
	grouping_vars <- c("STDAGE", "INVYR", "FORTYPCD", "foresttype", "biome","ownership", "managed", "ll_id")
	
	for (var in grouping_vars) {
		level_check <- fia %>%
			group_by(PID) %>%
			summarise(levels_count = n_distinct(.data[[var]])) %>%
			filter(levels_count > 1)
		
		if (nrow(level_check) > 0) {
			stop(paste("Error: More than one level per PID in grouping variable", var, "detected!"))
		}
	}
	
	# Now summarize data to species per plot level
	fia <- fia %>%
		
		# Group by species and plot level variables
		group_by(accepted_bin, PID, STDAGE, INVYR, FORTYPCD, foresttype, biome, ownership, managed, ll_id) %>%
		
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
	
	########## Calculate functional diversity indices and build database ##########
	
	# Summarise and scale trait data and format as matrix to calculate indices
	fun_traits <- dat %>%
		
		# Keep relevant columns 
		dplyr::select(accepted_bin, wood_density, leaf_n, bark_thickness, seed_dry_mass, 
									conduit_diam, shade_tolerance, specific_leaf_area,
									wmean_DIA_species, wmean_height_species, sum_BA_species) %>%
		
		# Group by species
		group_by(accepted_bin) %>%
		
		# Summarise traits to plot level 
		summarise(
			
			# Tree Means weighted by BA
			dia = mean(wmean_DIA_species, na.rm = TRUE),
			height = mean(wmean_height_species, na.rm = TRUE),
			
			# Species Means; equal to trait data
			wood_density = mean(wood_density, na.rm = TRUE),
			bark_thickness = mean(bark_thickness, na.rm = TRUE),
			conduit_diam = mean(conduit_diam, na.rm = TRUE),
			leaf_n = mean(leaf_n, na.rm = TRUE),
			specific_leaf_area = mean(specific_leaf_area, na.rm = TRUE),
			shade_tolerance = mean(shade_tolerance, na.rm = TRUE),
			seed_dry_mass = mean(seed_dry_mass, na.rm = TRUE),
			
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
		dplyr::select(accepted_bin, PID) %>%
		# Set presence to 1
		mutate(present = 1) %>%
		# Keep only species found in trait data 
		filter(accepted_bin %in% rownames(fun_traits)) %>%
		# Pivot to wider 
		pivot_wider(names_from = accepted_bin, values_from = present, values_fill = list(present = 0)) %>%
		# Format as matrix 
		column_to_rownames(var = "PID") %>%
		as.matrix()
	
	## Summarise to plot level, calculate resource use score, and functional indices
	
	# Summarise data to plot level 
	fia <- dat %>%
		
		# Define grouping variables that are constant for every plot 
		group_by(PID, STDAGE, INVYR, FORTYPCD, foresttype, biome, ownership, managed, ll_id) %>%
		rename(standage = STDAGE) %>%
		
		# Summarise to plot levels by weighted means (using BA)
		summarise(
			
			# Wood density 
			wmean_wood_density = weighted.mean(wood_density, sum_BA_species, na.rm = TRUE),
			
			# Bark thickness
			wmean_bark_thickness = weighted.mean(bark_thickness, sum_BA_species, na.rm = TRUE),
			
			# Conduit diameter
			wmean_conduit_diam = weighted.mean(conduit_diam, sum_BA_species, na.rm = TRUE),
			
			# Leaf nitrogen 
			wmean_leaf_n = weighted.mean(leaf_n, sum_BA_species, na.rm = TRUE),
			
			# Specific leaf area
			wmean_specific_leaf_area = weighted.mean(specific_leaf_area, sum_BA_species, na.rm = TRUE),
			
			# Seed dry mass
			wmean_seed_dry_mass = weighted.mean(seed_dry_mass, sum_BA_species, na.rm = TRUE),
			
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
		
		## Calculate functional diversity indices per plot
		
		# Functional divergence 
		left_join(fd_fdiv(fun_traits, fun_plots) %>% rename(PID = site), by = "PID") %>%
		rename(fun_div = FDiv) %>%
		
		# Functional dispersion
		left_join(fd_fdis(fun_traits, fun_plots) %>% rename(PID = site), by = "PID") %>%
		rename(fun_disp = FDis) %>%
		
		# Functional evenness
		left_join(fd_feve(fun_traits, fun_plots) %>% rename(PID = site), by = "PID") %>%
		rename(fun_even = FEve) %>%
		
		# Set format
		mutate_all(~ ifelse(is.nan(.), NA, .)) %>%
		as_tibble()
	
	########## Calculate Resource Use Score ##########
	
	# Function to normalize a vector to a range [0, 1]
	normalize <- function(x) {
		return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
	}
	
	# Function to invert normalized value for acquisitive traits
	invert_normalized <- function(x) {
		return(1 - x)
	}
	
	## Normalize trait values appropriately and calculate composite score 
	fia <- fia %>%
		
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
	
	## Return cleaned FIA data 
	return(fia)
}

# Wrap the FIA prep pipeline with safely
prep_fia_safe <- safely(prep_fia, otherwise = NA)

# Access FIA pathways 
fia_data <- tibble(pathway = list.files(path = "/Volumes/ritd-ag-project-rd01pr-dmayn10/fia_data/cleaned_state_data",
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
	dplyr::select(PID, state, standage, INVYR, FORTYPCD, foresttype, biome, ownership, managed, ll_id,
								wmean_wood_density, wmean_bark_thickness, wmean_conduit_diam, wmean_leaf_n, wmean_specific_leaf_area,
								wmean_seed_dry_mass, wmean_shade_tolerance, wmean_dia, wmean_height, 
								total_carbon_ag, total_biomass, total_basal_area, fun_div, fun_disp, fun_even, resource_use_score)

fia_clean <- fia_clean %>%
	
	# Integrate climate data 
	left_join(clim, by = "ll_id") %>%
	filter(complete.cases(standage)) %>%
	
	# Write file 
	write_csv(file = paste0("/Volumes/ritd-ag-project-rd01pr-dmayn10/merlin/data/fia_traits/plotlvl_data_", Sys.Date(),".csv"))

## Plot resource use score

plot_grid(
	
	plot_grid(
		
		fia_clean %>%
			filter(standage < 500) %>%
			ggplot(aes(x=standage, y=resource_use_score)) +
			geom_hex() +
			scale_fill_viridis_c() +
			geom_smooth(method = "lm", colour = "red") +
			theme_bw() +
			theme(legend.position = "none",
						text = element_text(family = "Palatino")) +
			ggtitle("All plots") +
			xlab("Stand age (years)") +
			ylab("Resource use score"),
		
		fia_clean %>%
			filter(standage < 500) %>%
			filter(complete.cases(standage, resource_use_score, foresttype)) %>%
			ggplot(aes(x=standage, y=resource_use_score)) +
			geom_hex() +
			scale_fill_viridis_c() +
			geom_smooth(method = "lm", colour = "red") +
			facet_wrap(~foresttype, scales = "free") +
			theme_bw() +
			theme(legend.position = "none",
						text = element_text(family = "Palatino")) +
			ggtitle("Forest types") +
			xlab("Stand age (years)") +
			ylab(NULL),
		
		ncol = 2, nrow = 1, rel_widths = c(.4, .6)),
	
	
	fia_clean %>%
		filter(standage < 500) %>%
		filter(complete.cases(standage, resource_use_score, biome)) %>%
		filter(biome != "Mangroves") %>%
		ggplot(aes(x=standage, y=resource_use_score, group = biome)) +
		geom_hex() +
		scale_fill_viridis_c() +
		geom_smooth(method = "lm", colour = "red") +
		facet_wrap(~biome, scales = "free") +
		theme_bw() +
		theme(legend.position = "none",
					text = element_text(family = "Palatino")) +
		ggtitle("Biomes") +
		xlab("Stand age (years)") +
		ylab(NULL), 
	
	nrow = 2, ncol = 1)

## plot functinal diversity indices

fia_clean %>%
	# managed yes or no
	filter(standage < 500) %>%
	# Keep PID, stand age, and mean traits
	select(PID, standage, fun_div, fun_disp, fun_even) %>%
	# Pivot longer
	pivot_longer(
		cols = c(fun_div, fun_disp, fun_even),       
		names_to = "Index",               
		values_to = "value") %>%
	# Cut off extreme values?
	filter(value > 0) %>%
	filter(value < 1) %>%
	# Plot indices against stand age 
	ggplot(aes(x=standage, y=value, group = Index)) +
	geom_hex() +
	scale_fill_viridis_c() +
	geom_smooth(method = "lm", colour = "red") +
	facet_wrap(~Index, scales = "free") +
	theme_bw() +
	theme(legend.position = "none",
				text = element_text(family = "Palatino")) +
	ggtitle("Functional diversity indices") +
	xlab("Stand age (years)") +
	ylab("Index values")

# plot weighted trait means 

fia_clean %>%
	filter(standage < 500) %>%
	# Keep PID, stand age, and mean traits
	select(PID, standage, wmean_wood_density, wmean_bark_thickness, wmean_leaf_n, wmean_conduit_diam,
				 wmean_shade_tolerance, wmean_specific_leaf_area, wmean_seed_dry_mass, wmean_height, wmean_dia) %>%
	# Pivot longer
	pivot_longer(
		cols = c(wmean_wood_density, wmean_bark_thickness, wmean_conduit_diam, wmean_leaf_n, 
						 wmean_shade_tolerance, wmean_specific_leaf_area, wmean_seed_dry_mass, wmean_height, wmean_dia),       
		names_to = "Trait",               
		values_to = "value") %>%
	
	# Plot mean traits against stand age 
	ggplot(aes(x=standage, y=value, group = Trait)) +
	geom_hex() +
	scale_fill_viridis_c() +
	geom_smooth(method = "lm", colour = "red") +
	facet_wrap(~Trait, scales = "free") +
	theme_bw() +
	theme(legend.position = "none",
				text = element_text(family = "Palatino")) +
	ggtitle("Weighted trait means") +
	xlab("Stand age (years)") +
	ylab("Values")

## plot total trait sums

fia_clean %>%
	# Stand age filter
	filter(standage < 500) %>%
	# Keep PID, stand age, and mean traits
	select(PID, standage, total_basal_area, total_biomass, total_carbon_ag, wmean_dia) %>%
	# Pivot longer
	pivot_longer(
		cols = c(total_basal_area, total_biomass, total_carbon_ag, wmean_dia),       
		names_to = "Trait",               
		values_to = "value") %>%
	# Plot mean traits against stand age 
	ggplot(aes(x=standage, y=value, group = Trait)) +
	geom_hex() +
	scale_fill_viridis_c() +
	geom_smooth(method = "lm", colour = "red") +
	facet_wrap(~Trait, scales = "free") +
	theme_bw() +
	theme(legend.position = "none",
				text = element_text(family = "Palatino")) +
	ggtitle("Total trait sums") +
	xlab("Stand age (years)") +
	ylab("Values")