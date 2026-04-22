################################################################################
## succession_traits: 06 — Predictability analysis (VEcv)
## Assesses whether trait expression becomes more predictable through
## succession and whether predictability diverges across environmental
## gradients.
##
## Method: repeated k-fold cross-validation within each leaf type, stratified
## post-hoc by environmental variable quantiles and stand age bins. VEcv
## (Variance Explained by Cross-validation) is computed per stratum × bin
## as a scale-free, cross-validated skill metric comparable across traits.
##
## We use VEcv over global model residuals (MSE) because:
##   - VEcv is scale-free and directly comparable across traits on different
##     scales, even after z-scoring
##   - Out-of-fold predictions avoid optimistic in-sample error estimates
##   - The CV framework matches standard predictive modelling practice
##
## Input:  data_processed/fia_traits_clean.rds
##         tables/perf_broadleaf.csv / perf_coniferous.csv  (hyperparameters)
##
## Output: tables/vecv_raw.rds        (all repeat × bin × stratum VEcv values)
##         tables/vecv_summary.rds    (medians + CIs)
##         plots/vecv/fig3_vecv.png   (RQ3 main figure)
##         plots/vecv/supp_divergence.png
##
## Author: M. Weiss @ Maynard Lab UCL / ETH Zürich
################################################################################

rm(list = ls())
set.seed(42)

# ── Libraries ─────────────────────────────────────────────────────────────────
library(ranger)
library(doParallel)
library(foreach)
library(tidyverse)

# ── Configuration ─────────────────────────────────────────────────────────────
# repeats = 30 gives stable VEcv distributions with manageable compute time.
# Set higher (e.g. 100) for final publication run if server time permits.
N_REPEATS       <- 30L
N_FOLDS         <- 10L
PROBS           <- c(0.25, 0.75)
STANDAGE_BREAKS <- seq(0, 150, by = 10)
MIN_BIN_N       <- 30L    # minimum plots per bin × stratum to report VEcv
PARALLEL        <- TRUE
N_CORES         <- 16L

PATH_DATA   <- "data_processed/fia_traits_clean.rds"
PATH_TABLES <- "tables"
PATH_PLOTS  <- "plots/vecv"

dir.create(PATH_PLOTS,  recursive = TRUE, showWarnings = FALSE)
dir.create(PATH_TABLES, showWarnings = FALSE)

source("scripts/functions.R")

# ── Vocabulary ────────────────────────────────────────────────────────────────
TRAIT_LABELS <- c(
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
TRAITS     <- names(TRAIT_LABELS)
COVARIATES <- c("standage", "temp_pc", "soil_pc", "rain_pc", "elevation", "soil_ph")
LEAF_TYPES <- c("broadleaf", "coniferous")
ENV_VARS   <- c("temp_pc", "soil_pc", "rain_pc", "elevation", "soil_ph")

ENV_LABELS <- c(
	"temp_pc"   = "Temperature PC",
	"soil_pc"   = "Soil water retention PC",
	"rain_pc"   = "Precipitation PC",
	"elevation" = "Elevation",
	"soil_ph"   = "Soil pH"
)


# ══════════════════════════════════════════════════════════════════════════════
# 1. Load data and hyperparameters
# ══════════════════════════════════════════════════════════════════════════════

data <- read_rds(PATH_DATA)

# Derive leaf type (consistent with all upstream scripts)
data <- data %>%
	mutate(
		leaf_type = case_when(
			biome_boreal_forests_or_taiga    == 1 |
				biome_temperate_conifer_forests  == 1 ~ "coniferous",
			biome_temperate_broadleaf_forests == 1 |
				biome_mediterranean_woodlands    == 1 ~ "broadleaf",
			TRUE ~ NA_character_
		)
	) %>%
	filter(!is.na(leaf_type))

message(sprintf(
	"Data: %d plots (%d broadleaf / %d coniferous)",
	nrow(data),
	sum(data$leaf_type == "broadleaf"),
	sum(data$leaf_type == "coniferous")
))

# Load pre-tuned hyperparameters
hyper_params <- map(LEAF_TYPES, function(lt) {
	read_csv(
		file.path(PATH_TABLES, sprintf("perf_%s.csv", lt)),
		show_col_types = FALSE
	) %>%
		dplyr::select(trait, num_trees, mtry, min_node_size)
}) %>% set_names(LEAF_TYPES)


# ══════════════════════════════════════════════════════════════════════════════
# 2. Run repeated CV skill computation
# ══════════════════════════════════════════════════════════════════════════════
# Parallelise across traits × leaf types. Each combination is independent.
# oof_skill_by_bins() handles the inner CV loop sequentially (num_threads = 1
# inside ranger) to avoid nested parallelism.

if (PARALLEL) {
	cl <- makeCluster(N_CORES)
	registerDoParallel(cl)
	message(sprintf("Parallel backend: %d cores", N_CORES))
} else {
	registerDoSEQ()
}

message(sprintf(
	"\nRunning repeated CV (v = %d, repeats = %d) for %d traits × %d leaf types...",
	N_FOLDS, N_REPEATS, length(TRAITS), length(LEAF_TYPES)
))

# Build all trait × leaf type combinations
jobs <- expand.grid(
	trait     = TRAITS,
	leaf_type = LEAF_TYPES,
	stringsAsFactors = FALSE
)

vecv_raw <- foreach(
	i         = seq_len(nrow(jobs)),
	.combine  = bind_rows,
	.packages = c("ranger", "dplyr", "purrr", "tidyr", "tibble"),
	.export   = c("oof_skill_by_bins", "fit_rf_model",
								"VEcv", "E1", "COVARIATES", "ENV_VARS",
								"STANDAGE_BREAKS", "PROBS", "N_FOLDS", "N_REPEATS")
) %dopar% {
	
	tr <- jobs$trait[i]
	lt <- jobs$leaf_type[i]
	
	df_lt <- data %>% dplyr::filter(leaf_type == lt)
	
	tryCatch(
		oof_skill_by_bins(
			trait           = tr,
			data            = df_lt,
			covariates      = COVARIATES,
			env_vars        = ENV_VARS,
			hyper_grid      = hyper_params[[lt]],
			standage_breaks = STANDAGE_BREAKS,
			probs           = PROBS,
			v               = N_FOLDS,
			repeats         = N_REPEATS
		) %>%
			dplyr::mutate(leaf_type = lt),
		error = function(e) {
			warning(sprintf("Failed: trait=%s, lt=%s\n  %s", tr, lt, conditionMessage(e)))
			NULL
		}
	)
}

if (PARALLEL) stopCluster(cl)

message(sprintf(
	"CV complete: %d rows across %d trait × leaf type combinations",
	nrow(vecv_raw), nrow(jobs)
))

# Filter bins with too few observations — VEcv is unreliable with small n
n_dropped <- sum(vecv_raw$n < MIN_BIN_N, na.rm = TRUE)
if (n_dropped > 0)
	message(sprintf(
		"  Dropping %d bin × stratum records with n < %d",
		n_dropped, MIN_BIN_N
	))

vecv_raw <- vecv_raw %>%
	filter(n >= MIN_BIN_N) %>%
	mutate(
		trait_label    = recode(trait,    !!!TRAIT_LABELS),
		variable_label = recode(variable, !!!ENV_LABELS),
		# Extract numeric stand age midpoint from bin label for plotting
		standage_mid   = as.numeric(str_extract(
			as.character(standage_bin), "(?<=\\[)\\d+"
		)) + 5
	)

write_rds(vecv_raw, file.path(PATH_TABLES, "vecv_raw.rds"))
message(sprintf("Raw VEcv saved to %s/vecv_raw.rds", PATH_TABLES))


# ══════════════════════════════════════════════════════════════════════════════
# 3. Summarise across repeats
# ══════════════════════════════════════════════════════════════════════════════
# Collapse the repeat dimension, retaining median and 95% CIs per
# trait × leaf type × environmental variable × stratum × stand age bin.

vecv_summary <- vecv_raw %>%
	group_by(trait, trait_label, leaf_type, variable, variable_label,
					 env_group, standage_bin, standage_mid) %>%
	summarise(
		n_med    = median(n,    na.rm = TRUE),
		VEcv_med = median(VEcv, na.rm = TRUE),
		VEcv_lwr = quantile(VEcv, 0.025, na.rm = TRUE),
		VEcv_upr = quantile(VEcv, 0.975, na.rm = TRUE),
		E1_med   = median(E1,   na.rm = TRUE),
		E1_lwr   = quantile(E1,  0.025, na.rm = TRUE),
		E1_upr   = quantile(E1,  0.975, na.rm = TRUE),
		.groups  = "drop"
	)

write_rds(vecv_summary, file.path(PATH_TABLES, "vecv_summary.rds"))


# ══════════════════════════════════════════════════════════════════════════════
# 4. Divergence index
# ══════════════════════════════════════════════════════════════════════════════
# ΔVECV = VEcv_high - VEcv_low per trait × variable × stand age bin × repeat.
# Positive: high-environment plots are more predictable at that age.
# Negative: low-environment plots are more predictable.
# The trajectory of ΔVECV through succession is the key RQ3 result:
#   - ΔVECV converging toward 0: environmental context matters less over time
#   - ΔVECV growing away from 0: assembly becomes more context-dependent

vecv_divergence <- vecv_raw %>%
	dplyr::select(trait, trait_label, leaf_type, variable, variable_label,
								env_group, standage_bin, standage_mid, repeat_id, VEcv) %>%
	pivot_wider(names_from = env_group, values_from = VEcv,
							names_prefix = "VEcv_") %>%
	mutate(delta_VEcv = VEcv_high - VEcv_low)

# Summarise divergence across repeats
divergence_summary <- vecv_divergence %>%
	group_by(trait, trait_label, leaf_type, variable, variable_label,
					 standage_bin, standage_mid) %>%
	summarise(
		delta_med  = median(delta_VEcv, na.rm = TRUE),
		delta_lwr  = quantile(delta_VEcv, 0.025, na.rm = TRUE),
		delta_upr  = quantile(delta_VEcv, 0.975, na.rm = TRUE),
		# Flag bins where CI excludes zero (robust divergence)
		sig_divergence = (delta_lwr > 0 | delta_upr < 0),
		.groups    = "drop"
	)

write_rds(divergence_summary, file.path(PATH_TABLES, "vecv_divergence.rds"))


# ══════════════════════════════════════════════════════════════════════════════
# 5. Figures
# ══════════════════════════════════════════════════════════════════════════════

# ── 5a. Main figure: VEcv trajectories (Figure 3 candidate) ──────────────────
# Shows VEcv vs stand age for high and low environmental groups.
# Averaged across environmental variables within each group for the main panel.
# Ribbons = 95% CI across repeats.
# Split by leaf type (facet row) and trait (facet column).
#
# This replaces the MSE figure in the current manuscript. Key improvements:
#   - VEcv is scale-free so traits are comparable on the same axis
#   - The y-axis has a natural interpretation (0 = no predictive skill,
#     1 = perfect prediction)
#   - CI ribbons reflect genuine out-of-fold uncertainty

# Average VEcv across environmental variables first (one line per group)
vecv_avg <- vecv_raw %>%
	group_by(trait, trait_label, leaf_type, env_group,
					 standage_bin, standage_mid, repeat_id) %>%
	summarise(VEcv = mean(VEcv, na.rm = TRUE), .groups = "drop") %>%
	group_by(trait, trait_label, leaf_type, env_group,
					 standage_bin, standage_mid) %>%
	summarise(
		VEcv_med = median(VEcv, na.rm = TRUE),
		VEcv_lwr = quantile(VEcv, 0.025, na.rm = TRUE),
		VEcv_upr = quantile(VEcv, 0.975, na.rm = TRUE),
		.groups  = "drop"
	) %>%
	mutate(
		env_group  = factor(env_group,
												levels = c("low", "high"),
												labels = c("Lower environmental quantile",
																	 "Upper environmental quantile")),
		leaf_type  = str_to_title(leaf_type)
	)

p_vecv <- ggplot(vecv_avg,
								 aes(x = standage_mid, colour = env_group, fill = env_group)) +
	geom_ribbon(aes(ymin = VEcv_lwr, ymax = VEcv_upr),
							alpha = 0.2, colour = NA) +
	geom_line(aes(y = VEcv_med), linewidth = 0.8) +
	geom_hline(yintercept = 0, linetype = "dashed",
						 colour = "grey50", linewidth = 0.3) +
	facet_grid(leaf_type ~ trait_label, scales = "free_y") +
	scale_colour_manual(
		values = c(
			"Lower environmental quantile" = "#1B9E77",
			"Upper environmental quantile" = "#D95F02"
		)
	) +
	scale_fill_manual(
		values = c(
			"Lower environmental quantile" = "#1B9E77",
			"Upper environmental quantile" = "#D95F02"
		)
	) +
	scale_x_continuous(breaks = c(0, 50, 100, 150)) +
	labs(
		x      = "Stand age (years)",
		y      = "VEcv (cross-validated variance explained)",
		colour = NULL, fill = NULL,
		title  = "Trait predictability through succession by environmental context"
	) +
	theme_bw(base_size = 9) +
	theme(
		legend.position  = "top",
		strip.text       = element_text(face = "bold", size = 7),
		strip.background = element_rect(fill = "white", colour = "black",
																		linewidth = 0.4),
		panel.grid.minor = element_blank(),
		axis.text.x      = element_text(size = 7)
	)

ggsave(
	file.path(PATH_PLOTS, "fig3_vecv.png"),
	plot = p_vecv,
	width = 260, height = 160, units = "mm", dpi = 400
)
message("Figure 3 VEcv saved")

# ── 5b. Divergence trajectories ───────────────────────────────────────────────
# Shows ΔVECV (high - low) vs stand age per trait × environmental variable.
# Dashed line at 0 = no divergence.
# CI ribbon excludes 0 where divergence is robust.
# Split by leaf type.
# This directly answers: "does environmental heterogeneity promote divergence
# or convergence in trait predictability over succession?"

p_divergence <- divergence_summary %>%
	mutate(leaf_type = str_to_title(leaf_type)) %>%
	ggplot(aes(x = standage_mid, y = delta_med,
						 colour = variable_label, fill = variable_label)) +
	geom_ribbon(aes(ymin = delta_lwr, ymax = delta_upr),
							alpha = 0.15, colour = NA) +
	geom_line(linewidth = 0.7) +
	geom_hline(yintercept = 0, linetype = "dashed",
						 colour = "grey40", linewidth = 0.4) +
	facet_grid(leaf_type ~ trait_label, scales = "free_y") +
	scale_colour_brewer(palette = "Set1", name = "Environmental variable") +
	scale_fill_brewer(  palette = "Set1", name = "Environmental variable") +
	scale_x_continuous(breaks = c(0, 50, 100, 150)) +
	labs(
		x     = "Stand age (years)",
		y     = expression(Delta * "VEcv (high \u2212 low environmental quantile)"),
		title = "Divergence in trait predictability across environmental gradients"
	) +
	theme_bw(base_size = 9) +
	theme(
		legend.position  = "bottom",
		legend.key.size  = unit(3, "mm"),
		strip.text       = element_text(face = "bold", size = 7),
		strip.background = element_rect(fill = "white", colour = "black",
																		linewidth = 0.4),
		panel.grid.minor = element_blank(),
		axis.text.x      = element_text(size = 7)
	)

ggsave(
	file.path(PATH_PLOTS, "supp_divergence.png"),
	plot = p_divergence,
	width = 260, height = 160, units = "mm", dpi = 350
)
message("Supplementary divergence figure saved")

# ── 5c. Environmental variable ranking: which drives most divergence? ─────────
# Mean |ΔVECV| averaged across traits and stand age bins per env variable.
# Bar chart, split by leaf type.

p_env_rank <- divergence_summary %>%
	group_by(leaf_type, variable_label) %>%
	summarise(
		mean_abs_delta = mean(abs(delta_med), na.rm = TRUE),
		.groups        = "drop"
	) %>%
	mutate(
		leaf_type      = str_to_title(leaf_type),
		variable_label = fct_reorder(variable_label, mean_abs_delta)
	) %>%
	ggplot(aes(x = mean_abs_delta, y = variable_label, fill = leaf_type)) +
	geom_col(position = "dodge", colour = "black",
					 linewidth = 0.3, alpha = 0.8) +
	scale_fill_manual(
		values = c("Broadleaf" = "#228B22", "Coniferous" = "#d95f02")
	) +
	labs(
		x    = "Mean |ΔVEcv| across traits and stand age bins",
		y    = NULL,
		fill = NULL,
		title = "Environmental drivers of predictability divergence"
	) +
	theme_bw(base_size = 10) +
	theme(
		legend.position  = "bottom",
		panel.grid.minor = element_blank(),
		panel.grid.major.y = element_blank()
	)

ggsave(
	file.path(PATH_PLOTS, "env_ranking.png"),
	plot = p_env_rank,
	width = 160, height = 100, units = "mm", dpi = 350
)
message("Environmental ranking figure saved")


# ══════════════════════════════════════════════════════════════════════════════
# 6. Console diagnostics
# ══════════════════════════════════════════════════════════════════════════════

message("\n── Mean VEcv by trait and leaf type ────────────────────────────────")
vecv_summary %>%
	group_by(trait_label, leaf_type) %>%
	summarise(mean_VEcv = mean(VEcv_med, na.rm = TRUE), .groups = "drop") %>%
	pivot_wider(names_from = leaf_type, values_from = mean_VEcv) %>%
	mutate(across(where(is.numeric), ~ round(., 3))) %>%
	print(n = Inf)

message("\n── Traits with strongest divergence (mean |ΔVEcv|) ─────────────────")
divergence_summary %>%
	group_by(trait_label, leaf_type) %>%
	summarise(
		mean_abs_delta = mean(abs(delta_med), na.rm = TRUE),
		prop_sig       = mean(sig_divergence, na.rm = TRUE),
		.groups        = "drop"
	) %>%
	arrange(desc(mean_abs_delta)) %>%
	mutate(across(where(is.numeric), ~ round(., 3))) %>%
	print(n = Inf)

message("\n── Environmental variables ranked by divergence ─────────────────────")
divergence_summary %>%
	group_by(variable_label, leaf_type) %>%
	summarise(
		mean_abs_delta = mean(abs(delta_med), na.rm = TRUE),
		.groups        = "drop"
	) %>%
	arrange(leaf_type, desc(mean_abs_delta)) %>%
	mutate(across(where(is.numeric), ~ round(., 3))) %>%
	print(n = Inf)

message("\n── 06_vecv.R complete ──────────────────────────────────────────────")
message(sprintf("Tables: %s/", PATH_TABLES))
message(sprintf("Plots:  %s/", PATH_PLOTS))