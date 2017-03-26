#!/usr/bin/Rscript

library(ggplot2)
library(scales)

starbeast_gene0_nodes = read.csv("gene0.trees.nodes.csv")[100009:1000008,]
starbeast_gene1_nodes = read.csv("gene1.trees.nodes.csv")[100009:1000008,]
simulated_gene0_nodes = read.csv("simulated.gene0.trees.nodes.csv")
simulated_gene1_nodes = read.csv("simulated.gene1.trees.nodes.csv")

starbeast_gene0_nodes$overall_rate = "Half"
simulated_gene0_nodes$overall_rate = "Half"
starbeast_gene1_nodes$overall_rate = "Double"
simulated_gene1_nodes$overall_rate = "Double"

starbeast_nodes = rbind(starbeast_gene0_nodes, starbeast_gene1_nodes)
simulated_nodes = rbind(simulated_gene0_nodes, simulated_gene1_nodes)

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
ggsave("gene_node_heights.pdf", units = "mm", height = 150, width = 150)
ggsave("gene_node_heights.png", units = "mm", height = 150, width = 150)
