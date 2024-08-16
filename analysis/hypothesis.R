library(tidyverse) 

regions <- rbind(
	
	c("A", 0, 1),
	c("A", 0, 0.6),
	c("A", 0.225, 0.375),
	c("A", 0.6, 0.4),
	c("A", 0.65, .75),
	c("A", 0.4, 1),
	
	c("B", 0.225, 0.375),
	c("B", 0.6, 0.4),
	c("B", 0.65, .75),
	c("B", 1, 0.4),
	c("B", 1, 0),
	c("B", 0.6, 0),
	
	c("C", 0, 0),
	c("C", 0, 0.575),
	c("C", 0.575, 0),
	
	c("D", 1, 0.425),
	c("D", 1, 1),
	c("D", 0.425, 1)
) %>% 
	as.data.frame(stringsAsFactors = FALSE) %>% 
	setNames(c("region", "x", "y")) %>% 
	mutate(
		x = as.numeric(x),
		y = as.numeric(y)
	)

# Create the plot
hypo <- ggplot() +
	geom_polygon(data = regions, aes(x = x, y = y, group = region, fill = region, linetype = region), 
							 color = "black", alpha = 0.25) +
	scale_fill_manual(values = c("A" = "midnightblue",
															 "B" = "forestgreen", 
															 "C" = "gray",
															 "D" = "gray")) +
	scale_linetype_manual(values = c("A" = "solid",
																	 "B" = "solid", 
																	 "C" = "dashed",
																	 "D" = "dashed")) +

	
	labs(x = "Successional Signal", y = "Determinism") +
	xlab("Signal across successional time") +
	ylab("Predictability through the environment") +
	theme_classic() +
	theme(legend.position = "none",
				text = element_text(size = 9),
				axis.text = element_blank(),
				axis.ticks = element_blank())

ggsave(filename = "/Users/serpent/Documents/MSc/Thesis/Code/analysis/plots/hypothesis.png",
			 plot = hypo, 
			 bg = "white",
			 width = 100, 
			 height = 100, 
			 units = "mm", 
			 dpi = 600)

