#!/usr/bin/Rscript

library(ggplot2)
library(scales)

starbeast_nodes = read.csv("species.trees.nodes.csv")[100009:1000008,]
simulated_nodes = read.csv("simulated.species.trees.nodes.csv")

starbeast_nodes$method = "StarBEAST2"
simulated_nodes$method = "Simulated"

combined_nodes = rbind(starbeast_nodes, simulated_nodes)
combined_nodes$method = factor(combined_nodes$method, levels = c("Simulated", "StarBEAST2"))

ggplot(combined_nodes, aes(x = branch_rate, linetype = method, color = method)) +
	geom_density(size = 1.0, adjust = 0.2) +
	scale_x_log10(limits = c(0.25, 4.0), breaks = 2^seq(-10, 10), labels = comma) +
	theme(legend.position = "bottom", legend.title = element_blank())
ggsave("species_branch_rates.pdf", units = "mm", height = 150, width = 150)
ggsave("species_branch_rates.png", units = "mm", height = 150, width = 150)
