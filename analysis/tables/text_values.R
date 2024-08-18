library(dplyr)
library(tidyr)

pdp_data <- read_csv("/Volumes/ritd-ag-project-rd01pr-dmayn10/merlin/traits_output/including_standage/tables/pdp_data.csv")

# Calculate the percentage difference and fold change for each group and trait
pdp_data <- pdp_data %>%
	group_by(trait, group) %>%
	summarise(
		yhat_start = first(yhat),
		yhat_end = last(yhat),
		perc_diff = ((yhat_end - yhat_start) / abs(yhat_start)) * 100,  # Calculate percentage difference
		fold_change = yhat_end / yhat_start  # Calculate fold change
	)

# Calculate the overall mean of the percentage difference and fold change across all groups
overall_mean_perc_diff <- mean(pdp_data$perc_diff)
overall_sd_perc_diff <- sd(pdp_data$perc_diff)
overall_mean_fold_change <- mean(pdp_data$fold_change)
overall_sd_fold_change <- sd(pdp_data$fold_change)

cat("Overall mean percentage difference between first and last trait values across all groups: ", 
		round(overall_mean_perc_diff, 2), "±", round(overall_sd_perc_diff, 2), " (", round(overall_mean_fold_change, 2), "±", round(overall_sd_fold_change, 2), " fold)\n")

# Identify the traits with the greatest and lowest mean percentage difference across standage
trait_diff_means <- pdp_data %>%
	group_by(trait) %>%
	summarise(mean_perc_diff = mean(perc_diff),
						mean_fold_change = mean(fold_change)) %>%
	arrange(desc(mean_perc_diff))

max_trait_perc_diff <- max(trait_diff_means$mean_perc_diff)
min_trait_perc_diff <- min(trait_diff_means$mean_perc_diff)
max_trait_fold_change <- max(trait_diff_means$mean_fold_change)
min_trait_fold_change <- min(trait_diff_means$mean_fold_change)
max_trait_name <- trait_diff_means$trait[which.max(trait_diff_means$mean_perc_diff)]
min_trait_name <- trait_diff_means$trait[which.min(trait_diff_means$mean_perc_diff)]

cat("Trait with the greatest mean percentage difference across standage: ", 
		max_trait_name, "(", round(max_trait_perc_diff, 2), "%, ", round(max_trait_fold_change, 2), " fold)\n")
cat("Trait with the lowest mean percentage difference across standage: ", 
		min_trait_name, "(", round(min_trait_perc_diff, 2), "%, ", round(min_trait_fold_change, 2), " fold)\n")



pdp_data <- read_csv("/Volumes/ritd-ag-project-rd01pr-dmayn10/merlin/traits_output/including_standage/tables/pdp_data.csv")

# Define the stability threshold (e.g., 10%)
stability_threshold <- 25

# Calculate the relative change
trait_stability <- pdp_data %>%
	group_by(trait, group, variable) %>%
	summarize(relative_change = (last(yhat) - first(yhat)) / abs(first(yhat)) * 100) %>%
	mutate(trend = case_when(
		relative_change > stability_threshold ~ "increased",
		relative_change < -stability_threshold ~ "decreased",
		TRUE ~ "stable"
	))

# Summarize the results
trait_stability_summary <- trait_stability %>%
	group_by(trait, variable, trend) %>%
	summarize(n = n()) %>%
	arrange(variable, trait, desc(n))

print(trait_stability_summary)




pdp_data <- read_csv("/Volumes/ritd-ag-project-rd01pr-dmayn10/merlin/traits_output/including_standage/tables/pdp_data.csv")

# Calculate the percentage difference and fold change for each group and trait
pdp_data <- pdp_data %>%
	filter(variable == "Temperature (PC)") %>%
	group_by(trait, group) %>%
	summarise(
		yhat_start = first(yhat),
		yhat_end = last(yhat),
		perc_diff = ((yhat_end - yhat_start) / abs(yhat_start)) * 100,  # Calculate percentage difference
		fold_change = yhat_end / yhat_start  # Calculate fold change
	)

# Calculate the overall mean of the percentage difference and fold change across all groups
overall_mean_perc_diff <- mean(pdp_data$perc_diff)
overall_sd_perc_diff <- sd(pdp_data$perc_diff)
overall_mean_fold_change <- mean(pdp_data$fold_change)
overall_sd_fold_change <- sd(pdp_data$fold_change)

cat("Overall mean percentage difference between first and last trait values across all groups: ", 
		round(overall_mean_perc_diff, 2), "±", round(overall_sd_perc_diff, 2), " (", round(overall_mean_fold_change, 2), "±", round(overall_sd_fold_change, 2), " fold)\n")

pdp_data <- read_csv("/Volumes/ritd-ag-project-rd01pr-dmayn10/merlin/traits_output/including_standage/tables/pdp_data.csv")

# Define the stability threshold (e.g., 10%)
stability_threshold <- 25

# Calculate the relative change
trait_stability <- pdp_data %>%
	filter(variable == "Temperature (PC)") %>%
	group_by(trait, group, variable) %>%
	summarize(relative_change = (last(yhat) - first(yhat)) / abs(first(yhat)) * 100) %>%
	mutate(trend = case_when(
		relative_change > stability_threshold ~ "increased",
		relative_change < -stability_threshold ~ "decreased",
		TRUE ~ "stable"
	))

# Summarize the results
trait_stability_summary <- trait_stability %>%
	group_by(trait, variable, trend) %>%
	summarize(n = n()) %>%
	arrange(variable, trait, desc(n))

print(trait_stability_summary)


pdp_data <- read_csv("/Volumes/ritd-ag-project-rd01pr-dmayn10/merlin/traits_output/including_standage/tables/pdp_data.csv")

# Calculate the percentage difference and fold change for each group and trait
pdp_data <- pdp_data %>%
	filter(variable == "Precipitation (PC)") %>%
	group_by(trait, group) %>%
	summarise(
		yhat_start = first(yhat),
		yhat_end = last(yhat),
		perc_diff = ((yhat_end - yhat_start) / abs(yhat_start)) * 100,  # Calculate percentage difference
		fold_change = yhat_end / yhat_start  # Calculate fold change
	)

# Calculate the overall mean of the percentage difference and fold change across all groups
overall_mean_perc_diff <- mean(pdp_data$perc_diff)
overall_sd_perc_diff <- sd(pdp_data$perc_diff)
overall_mean_fold_change <- mean(pdp_data$fold_change)
overall_sd_fold_change <- sd(pdp_data$fold_change)

cat("Overall mean percentage difference between first and last trait values across all groups: ", 
		round(overall_mean_perc_diff, 2), "±", round(overall_sd_perc_diff, 2), " (", round(overall_mean_fold_change, 2), "±", round(overall_sd_fold_change, 2), " fold)\n")



pdp_data <- read_csv("/Volumes/ritd-ag-project-rd01pr-dmayn10/merlin/traits_output/including_standage/tables/pdp_data.csv")

# Define the stability threshold (e.g., 10%)
stability_threshold <- 25

# Calculate the relative change
trait_stability <- pdp_data %>%
	filter(variable == "Precipitation (PC)") %>%
	group_by(trait, group, variable) %>%
	summarize(relative_change = (last(yhat) - first(yhat)) / abs(first(yhat)) * 100) %>%
	mutate(trend = case_when(
		relative_change > stability_threshold ~ "increased",
		relative_change < -stability_threshold ~ "decreased",
		TRUE ~ "stable"
	))

# Summarize the results
trait_stability_summary <- trait_stability %>%
	group_by(trait, variable, trend) %>%
	summarize(n = n()) %>%
	arrange(variable, trait, desc(n))

print(trait_stability_summary)
arrange(trait_stability_summary, n)

