#!/usr/bin/Rscript

library(ggplot2)

starbeast_gene0_topologies = read.csv("gene0.trees.topologies.csv")[11113:111112,]
starbeast_gene1_topologies = read.csv("gene1.trees.topologies.csv")[11113:111112,]
simulated_gene0_topologies = read.csv("simulated.gene0.trees.topologies.csv")
simulated_gene1_topologies = read.csv("simulated.gene1.trees.topologies.csv")

colnames(starbeast_gene0_topologies) = c("gene_tree_id", "gene_topology", "gene_topology_code")
colnames(starbeast_gene1_topologies) = c("gene_tree_id", "gene_topology", "gene_topology_code")
colnames(simulated_gene0_topologies) = c("gene_tree_id", "gene_topology", "gene_topology_code")
colnames(simulated_gene1_topologies) = c("gene_tree_id", "gene_topology", "gene_topology_code")

starbeast_gene0_topologies$method = "StarBEAST2"
starbeast_gene1_topologies$method = "StarBEAST2"
simulated_gene0_topologies$method = "Simulated"
simulated_gene1_topologies$method = "Simulated"

starbeast_gene0_topologies$rep = 0
starbeast_gene1_topologies$rep = 1
simulated_gene0_topologies$rep = 0
simulated_gene1_topologies$rep = 1

combined_topologies = rbind(starbeast_gene0_topologies, starbeast_gene1_topologies, simulated_gene0_topologies, simulated_gene1_topologies)
combined_topologies$method = factor(combined_topologies$method, levels = c("Simulated", "StarBEAST2"))
combined_topologies$rep = factor(combined_topologies$rep, levels = c(0, 1))

gene_topology_table <- table(combined_topologies$gene_topology_code)
gene_topology_order <- names(gene_topology_table[order(gene_topology_table, decreasing = TRUE)])
combined_topologies$gene_topology_code <- factor(combined_topologies$gene_topology_code, levels = gene_topology_order)

ggplot(combined_topologies, aes(x = gene_topology_code, color = method, shape = rep)) + geom_point(size = 1.0, stat = "count") +
	scale_y_log10(breaks = c(500, 1000, 2000), limits = c(400, 2500)) +
	labs(x = "Unique tree topology (ordered by frequency)", y = "Tree topology frequency") +
	theme(legend.position = "bottom", legend.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave("gene_topology_frequencies.pdf", units = "mm", height = 150, width = 150)
ggsave("gene_topology_frequencies.png", units = "mm", height = 150, width = 150)
