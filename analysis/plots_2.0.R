library(ggthemes)
library(ggalt)    
library(forcats)
library(ggrepel)
library(scales)
library(tidyverse)
path_out <- "/Users/merlin/Documents/MSc/Thesis/Code/output/plots/fig2"

## FIG 2 PANEL A
sum_cat <- shap_summary_boot %>%
  group_by(trait, category) %>%
  summarise(med = median(weighted_shap),
            lwr = quantile(weighted_shap, 0.025),
            upr = quantile(weighted_shap, 0.975),
            .groups = "drop")
wide <- sum_cat %>% select(trait, category, med) %>%
  pivot_wider(names_from = category, values_from = med)
trait_order2 <- wide %>%
  mutate(delta_med = `Successional Filtering` - `Environmental Filtering`) %>%
  arrange(desc(delta_med)) %>% pull(trait)
sum_cat <- sum_cat %>% mutate(trait = factor(trait, levels = trait_order2))
wide    <- wide    %>% mutate(trait = factor(trait, levels = trait_order2))

# Dumbbell of medians
p_db <- ggplot(wide,
               aes(y = trait, x = `Environmental Filtering`, xend = `Successional Filtering`)) +
  geom_dumbbell(size = 1.2, color = "grey70",
                colour_x = "black", colour_xend = "darkgreen", size_x = 2.5, size_xend = 2.5) +
  labs(x = "Weighted SHAP (median across bootstraps)", y = NULL) +
  theme_few() +
  theme(
    strip.text = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),  
    legend.position = "bottom",
    plot.margin = margin(t = 5, r = 20, b = 5, l = 10),
    panel.background = element_blank(),
    plot.background = element_blank(),
    strip.background = element_blank()
  )

p_db
ggsave(file.path(path_out, "fig2_panel_a.svg"), p_db,
       width = 150, height = 120, units = "mm")


## FIGURE 2 PANEL B
pdp_stats <- bootstrap_pdp %>%
  mutate(group = recode_group(group)) %>%
  # Compute bootstrapped slopes and intercepts
  group_by(trait, variable, group, iteration) %>%
  summarise(slope = coef(lm(yhat ~ standage, data = cur_data()))[2],
            intercept = coef(lm(yhat ~ standage, data = cur_data()))[1],
            .groups = "drop") %>%
  # Pivot to compare high vs. low environmental groups
  pivot_wider(names_from = group, values_from = c(slope, intercept), names_prefix = "") %>%
  # Compute absolute differences in slopes & intercepts
  mutate(slope_diff = abs(slope_high - slope_low),
         intercept_diff = abs(intercept_high - intercept_low)) %>%
  mutate(trait = factor(trait, levels = names(trait_labels), labels = trait_labels))

## Plot results along a shap weighted mean 
pdp_stats_ellipses <- bootstrap_shap %>%
  # Compute total SHAP per trait-variable-iteration
  dplyr::group_by(trait, variable, bootstrap_id) %>%
  dplyr::summarise(total_shap = sum(abs(shap_value)), .groups = "drop") %>%
  # Filter to environmental variables and recode names
  filter(variable %in% c("elevation", "rain_pc", "soil_pc", "soil_ph", "temp_pc")) %>%
  mutate(variable = recode(variable,
                           elevation = "Elevation",
                           rain_pc = "Precipitation (PC)",
                           soil_pc = "Soil - Water Retention (PC)",
                           soil_ph = "Soil pH",
                           temp_pc = "Temperature (PC)"),
         trait = recode(trait,
                        bark_thickness = "Bark Thickness",
                        conduit_diam = "Conduit Diameter",
                        height = "Tree Height",
                        leaf_n = "Leaf Nitrogen",
                        seed_dry_mass = "Seed Dry Mass",
                        shade_tolerance = "Shade Tolerance",
                        specific_leaf_area = "Specific Leaf Area",
                        wood_density = "Wood Density")) %>%
  
  # Compute per-bootstrap SHAP weights
  group_by(trait, bootstrap_id) %>%
  mutate(weight = total_shap / sum(total_shap)) %>%
  rename(iteration = bootstrap_id) %>%
  select(trait, iteration, variable, weight) %>%
  # Join to PDP slope/intercept diffs and compute weighted means
  right_join(pdp_stats, by = c("trait", "variable", "iteration")) %>%
  group_by(trait, iteration) %>%
  summarise(
    slope_diff = sum(slope_diff * weight),
    intercept_diff = sum(intercept_diff * weight),
    variable = "SHAP-Weighted Mean",
    .groups = "drop"
  ) %>%
  
  # Combine with unweighted PDP stats and set facet levels
  bind_rows(pdp_stats, .) %>%
  mutate(variable = factor(variable, levels = c(
    "SHAP-Weighted Mean", "Elevation",
    "Temperature (PC)", "Precipitation (PC)",
    "Soil - Water Retention (PC)", "Soil pH"
  )))

pdp_stats_summary <- pdp_stats %>%
  group_by(trait, variable) %>%
  summarise(
    
    # Wilcoxon paired test
    slope_p_value = wilcox.test(slope_high, slope_low, paired = TRUE, exact = FALSE)$p.value,
    intercept_p_value = wilcox.test(intercept_high, intercept_low, paired = TRUE, exact = FALSE)$p.value,
    
    # Compute Cohen’s d for effect size
    slope_cohen_d = effsize::cohen.d(slope_high, slope_low, paired = TRUE, hedges.correction = TRUE)$estimate,
    intercept_cohen_d = effsize::cohen.d(intercept_high, intercept_low, paired = TRUE, hedges.correction = TRUE)$estimate,
    
    # Significance flags
    significant_slope = slope_p_value < 0.05,
    significant_intercept = intercept_p_value < 0.05,
    
    slope_median = median(slope_diff),
    slope_lower = quantile(slope_diff, 0.025),
    slope_upper = quantile(slope_diff, 0.975),
    intercept_median = median(intercept_diff),
    intercept_lower = quantile(intercept_diff, 0.025),
    intercept_upper = quantile(intercept_diff, 0.975),
    .groups = "drop") 

# Summarise bootstrap distributions per trait × variable
pdp_stats_summary_extended <- pdp_stats_ellipses %>%
  dplyr::group_by(trait, variable) %>%
  dplyr::summarise(
    intercept_median = median(intercept_diff, na.rm = TRUE),
    intercept_lwr    = quantile(intercept_diff, 0.025, na.rm = TRUE),
    intercept_upr    = quantile(intercept_diff, 0.975, na.rm = TRUE),
    slope_median     = median(slope_diff, na.rm = TRUE),
    slope_lwr        = quantile(slope_diff, 0.025, na.rm = TRUE),
    slope_upr        = quantile(slope_diff, 0.975, na.rm = TRUE),
    r_med            = sqrt(intercept_median^2 + slope_median^2),  # optional magnitude
    .groups = "drop"
  )

# Build the plot (Fig 3)
pdp_plot <- pdp_stats_ellipses %>% ggplot(aes(x = intercept_diff, y = slope_diff, fill = trait, shape = trait)) +
  # reference lines at zero (if you switch to signed deltas this becomes very informative)
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey55") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey55") +
  
  stat_ellipse(aes(group = trait), geom = "polygon", alpha = 0.2, color = NA) +
  geom_point(data = pdp_stats_summary_extended,
             aes(x = intercept_median, y = slope_median, fill = trait, shape = trait),
             size = 3, color = "black", stroke = 0.5, inherit.aes = FALSE) +
  scale_fill_viridis_d(name = NULL) +
  scale_shape_manual(values = c(
    "Bark Thickness" = 21, "Conduit Diameter" = 22, "Tree Height" = 23,
    "Leaf Nitrogen" = 24, "Seed Dry Mass" = 25, "Shade Tolerance" = 21,
    "Specific Leaf Area" = 22, "Wood Density" = 23), name = NULL) +
  guides(fill = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) + 
  facet_wrap(~variable, scales = "fixed", ncol = 2) +
  labs(
    x = expression("Initial Environmental Filtering (" * Delta * " intercept)"),
    y = expression("Spatio-temporal Interaction (" * Delta * " slope)")) +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "white", colour = "black", linewidth = .75),
    legend.position = "top",
    legend.box = "horizontal",
    legend.key = element_rect(fill = "white", colour = "black"))

pdp_plot
ggsave(file.path(path_out, "fig2_panel_b.svg"), pdp_plot + theme(legend.position = "none"),
       width = 150, height = 120, units = "mm")

ggsave(file.path(path_out, "fig2_panel_b_legend.svg"), pdp_plot,
       width = 150, height = 120, units = "mm")
