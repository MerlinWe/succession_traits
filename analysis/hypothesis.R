library(ggplot2)

quadrants <- data.frame(
	x = c(1, 1, 3, 3),
	y = c(1, 3, 1, 3),
	label = c("Strong environmental filtering,\nhigh stochasitcity", 
						"Strong environmental filtering,\nhigh determinism", 
						"No environmental filtering\n strongly affected by biotic drivers",
						"Partial environmental filtering,\nstrong biotic residual effect"))

hypo <- ggplot() +
	geom_rect(data = quadrants, aes(xmin = ifelse(x == 1, -1.95, 0.05), xmax = ifelse(x == 3, 1.95, -0.05), 
																	ymin = ifelse(y == 1, -1.95, 0.05), ymax = ifelse(y == 3, 1.95, -0.05), 
																	fill = label), alpha = 0.3, colour = "black", linewidth = .3) +
	geom_text(data = quadrants, aes(x = ifelse(x == 1, -1, 1), y = ifelse(y == 1, -1, 1), 
																	label = label), hjust = "center", vjust = "center", size = 2.5, family = "sans") +
	scale_fill_manual(values = c("palegreen","darkgreen", "midnightblue", "lightskyblue")) +
	xlim(-2, 2) +
	ylim(-2, 2) +
	xlab("Signal across successional time") +
	ylab("Determinism") +
	theme_classic() +
	theme(legend.position = "none",
				text = element_text(size = 9),
				axis.text = element_blank(),
				axis.ticks = element_blank())

# hypo <- ggplot() +
# 	geom_rect(data = quadrants, aes(xmin = ifelse(x == 1, -1.99, 0.05), xmax = ifelse(x == 3, 1.99, -0.05), 
# 																	ymin = ifelse(y == 1, -1.99, 0.05), ymax = ifelse(y == 3, 1.99, -0.05), 
# 																	fill = label), alpha = 0.2, colour = "black", linewidth = .3) +
# 	geom_text(data = quadrants, aes(x = ifelse(x == 1, -1, 1), y = ifelse(y == 1, -1, 1), 
# 																	label = label), hjust = "center", vjust = "center", size = 3.5, family = "sans") +
# 	scale_fill_manual(values = c("darkgreen","midnightblue", "palegreen", "lightskyblue")) +
# 	xlim(-2, 2) +
# 	ylim(-2, 2) +
# 	xlab("Trait change over succession") +
# 	ylab("Determinism") +
# 	theme_bw() +
# 	theme(legend.position = "none",
# 				text = element_text(size = 9),
# 				axis.text = element_blank(),
# 				axis.ticks = element_blank())

ggsave(filename = "/Users/serpent/Documents/MSc/Thesis/Code/analysis/plots/hypothesis.png",
			 plot = hypo, 
			 bg = "white",
			 width = 100, 
			 height = 80, 
			 units = "mm", 
			 dpi = 600)

