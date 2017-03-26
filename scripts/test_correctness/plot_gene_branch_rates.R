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
combined_nodes = subset(combined_nodes, node_rank < 8)
combined_nodes$method = factor(combined_nodes$method, levels = c("Simulated", "StarBEAST2"))
combined_nodes$overall_rate = factor(combined_nodes$overall_rate, levels = c("Half", "Double"))

ggplot(combined_nodes, aes(x = branch_rate, linetype = method, color = method, group = paste(method, overall_rate))) +
	geom_density(size = 1.0, adjust = 0.2) +
	scale_x_log10(limits = c(0.25, 4.0), breaks = 2^seq(-10, 10), labels = comma) +
	labs(x = "Branch substitution rate", y = "Density") +
	theme(legend.position = "bottom", legend.title = element_blank())
ggsave("gene_branch_rates.pdf", units = "mm", height = 150, width = 150)
ggsave("gene_branch_rates.png", units = "mm", height = 150, width = 150)
