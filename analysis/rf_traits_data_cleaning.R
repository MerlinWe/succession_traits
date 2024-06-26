## Prep a clean data file for random forest trait modelling 

fia_clean <- list.files("/Volumes/ritd-ag-project-rd01pr-dmayn10/merlin/data/fia_traits", full.names = TRUE) %>%
	file.info() %>%
	as_tibble(rownames = "file") %>%
	arrange(desc(mtime)) %>%
	slice(1) %>%
	pull(file) %>%
	read_csv() %>% 
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
		across(starts_with("fun_"), as.numeric),
		resource_use_score = as.numeric(resource_use_score),
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
		ecoregion = as.integer(ecoregion),
		across(starts_with("sand_content_"), as.numeric),
		across(starts_with("soil_ph_"), as.numeric),
		across(starts_with("water_capacity_"), as.numeric))

# Add ecoregions 
ecoregions <- read_csv("/Volumes/ritd-ag-project-rd01pr-dmayn10/merlin/data/composite/resolve_ecoregions_legend.csv", locale = locale(encoding = "UTF-8")) %>%
	mutate(
		Name = iconv(Name, from = "UTF-8", to = "ASCII//TRANSLIT"), 
		Name = stri_replace_all_regex(tolower(Name), "\\s", "_")) %>%
	mutate(Code = as.integer(Code))

## Final adjustments and filtering 
fia_clean <- fia_clean %>%
	left_join(ecoregions, by = c("ecoregion" = "Code")) %>%
	mutate(ecoregion = Name) %>%
	select(-Name) %>%
	select(starts_with("wmean_"), standage, ecoregion, managed, 
				 annual_mean_temperature, annual_precipitation, temperature_seasonality, mean_diurnal_range, 
				 min_temperature_of_coldest_month, max_temperature_of_warmest_month,
				 elevation, pop_density, sand_content_015cm, soil_ph_015cm, water_capacity_015cm) %>%
	filter(complete.cases(.)) %>%
	as_tibble() %>%
	write_csv("traits_rf_clean.csv", file = "/Volumes/ritd-ag-project-rd01pr-dmayn10/merlin/data/fia_traits/traits_rf_clean.csv")
