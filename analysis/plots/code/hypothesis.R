library(tidyverse)

# ----- Main Figure ------

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
																	fill = label), alpha = .1, colour = "transparent", linetype = "solid", linewidth = .3) +
	scale_fill_viridis_d(direction = -1) +
	xlim(-2, 2) +
	ylim(-2, 2) +
	xlab("Signal across successional time") +
	ylab("Predictability through the environment") +
	theme_classic() +
	theme(legend.position = "none",
				text = element_text(size = 10, family = "sans"),
				axis.text = element_blank(),
				axis.ticks = element_blank())

hypo
ggsave(filename = "/Users/serpent/Documents/MSc/Thesis/Code/analysis/plots/hypthesis/hypothesis.png",
			 plot = hypo, 
			 bg = "white",
			 width = 100, 
			 height = 90, 
			 units = "mm", 
			 dpi = 600)

## ----- Sub Plots ------

# Upper left
set.seed(43) 
points <- data.frame(x = rnorm(30, mean = 0), y = rnorm(30, mean = 5, sd = 0.1))  
p1 <- ggplot() +
	geom_point(data = points, aes(x, y), size = 3, shape = 20) + 
	geom_line(aes(x=c(-2.5, 2), y=c(5,5)), linetype = "solid", color = "black", size = 1) +
	labs(x = "Successional Time", y = "Trait Expression") +
	ylim(c(4,6))+
	theme_classic() +
	theme(legend.position = "none",
				text = element_text(size = 12, family = "sans"),
				axis.ticks = element_blank(),
				axis.text = element_blank(),
				panel.background = element_blank(), 
				plot.background = element_blank())
p1
ggsave(filename = "/Users/serpent/Documents/MSc/Thesis/Code/analysis/plots/hypthesis/upper_left.png",
			 plot = p1, 
			 bg = "transparent",
			 width = 100, 
			 height = 70, 
			 units = "mm", 
			 dpi = 600)


# Upper right
set.seed(123) 
data <- data.frame(x = rnorm(30, mean = 0), y = rnorm(30, mean = 5, sd = 0.05))  
data$y <- data$y + 0.5 * data$x  
p2 <- ggplot(data, aes(x, y)) +
	geom_point(size = 3, shape = 20) + 
	geom_smooth(method = "lm", se = FALSE, color = "black", size = 1) +  
	labs(x = "Successional Time", y = "Trait Expression") +
	theme_classic() +
	theme(legend.position = "none",
				text = element_text(size = 12, family = "sans"),
				axis.ticks = element_blank(),
				axis.text = element_blank(),
				panel.background = element_blank(), 
				plot.background = element_blank())
p2
ggsave(filename = "/Users/serpent/Documents/MSc/Thesis/Code/analysis/plots/hypthesis/upper_right.png",
			 plot = p2, 
			 bg = "transparent",
			 width = 100, 
			 height = 70, 
			 units = "mm", 
			 dpi = 600)


# Lower left
set.seed(42)
points <- data.frame(x = runif(30, min = -2.5, max = 2), y = runif(30, min = 4.75, max = 5.25))
p3 <- ggplot() +
	geom_point(data = points, aes(x, y), size = 3, shape = 20) +
	geom_line(aes(x=c(-2.5, 2), y=c(5,5)), linetype = "solid", color = "black", size = 1) +
	labs(x = "Successional Time", y = "Trait Expression") +
	ylim(c(4.75,5.25)) +
	theme_classic() +
	theme(
		legend.position = "none",
		text = element_text(size = 12, family = "sans"),
		axis.ticks = element_blank(),
		axis.text = element_blank(),
		panel.background = element_blank(),
		plot.background = element_blank()
	)

p3
ggsave(filename = "/Users/serpent/Documents/MSc/Thesis/Code/analysis/plots/hypthesis/lower_left.png",
			 plot = p3, 
			 bg = "transparent",
			 width = 100, 
			 height = 70, 
			 units = "mm", 
			 dpi = 600)


# Lower right
set.seed(42) 
data <- data.frame(x = rnorm(30, mean = 0), y = rnorm(30, mean = 5, sd = .7))  
data$y <- data$y + 0.5 * data$x  
p4 <- ggplot(data, aes(x, y)) +
	geom_point(size = 3, shape = 20) + 
	geom_smooth(method = "lm", se = FALSE, color = "black", size = 1) +  
	labs(x = "Successional Time", y = "Trait Expression") +
	theme_classic() +
	theme(legend.position = "none",
				text = element_text(size = 12, family = "sans"),
				axis.ticks = element_blank(),
				axis.text = element_blank(),
				panel.background = element_blank(), 
				plot.background = element_blank())
p4
ggsave(filename = "/Users/serpent/Documents/MSc/Thesis/Code/analysis/plots/hypthesis/lower_right.png",
			 plot = p4, 
			 bg = "transparent",
			 width = 100, 
			 height = 70, 
			 units = "mm", 
			 dpi = 600)
