#!/usr/bin/Rscript

library(ggplot2)
library(scales)

starbeast_nodes = read.csv("species.trees.nodes.csv")[100009:1000008,]
simulated_nodes = read.csv("simulated.species.trees.nodes.csv")

starbeast_nodes$method = "StarBEAST2"
simulated_nodes$method = "Simulated"

combined_nodes = rbind(starbeast_nodes, simulated_nodes)
combined_nodes = subset(combined_nodes, node_rank > 4)
combined_nodes$method = factor(combined_nodes$method, levels = c("Simulated", "StarBEAST2"))
combined_nodes$node_rank = as.factor(combined_nodes$node_rank)

ggplot(combined_nodes, aes(x = node_height, linetype = method, color = method, group = paste(method, node_rank))) +
	geom_density(size = 1.0, adjust = 1.0) +
	scale_x_log10(limits = c(10.0^-6, 10.0^-1), breaks = 10^seq(-10, 10), labels = comma) +
	labs(x = "Node height", y = "Density") +
	theme(legend.position = "bottom", legend.title = element_blank())
ggsave("species_node_heights.pdf", units = "mm", height = 150, width = 150)
ggsave("species_node_heights.png", units = "mm", height = 150, width = 150)
