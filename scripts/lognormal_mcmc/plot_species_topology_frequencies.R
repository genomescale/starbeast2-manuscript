#!/usr/bin/Rscript

library(ggplot2)

starbeast_topologies = read.csv("species.trees.topologies.csv")[11113:111112,]
simulated_topologies = read.csv("simulated.species.trees.topologies.csv")

colnames(starbeast_topologies) = c("species_tree_id", "species_topology", "species_topology_code")
colnames(simulated_topologies) = c("species_tree_id", "species_topology", "species_topology_code")

starbeast_topologies$method  = "StarBEAST2"
simulated_topologies$method  = "Simulated"

combined_topologies = rbind(starbeast_topologies, simulated_topologies)
combined_topologies$method = factor(combined_topologies$method, levels = c("Simulated", "StarBEAST2"))

species_topology_table <- table(combined_topologies$species_topology_code)
species_topology_order <- names(species_topology_table[order(species_topology_table, decreasing = TRUE)])
combined_topologies$species_topology_code <- factor(combined_topologies$species_topology_code, levels = species_topology_order)

ggplot(combined_topologies, aes(x = species_topology_code, color = method)) + geom_point(size = 1.0, stat = "count") +
	scale_y_log10(breaks = c(500, 1000, 2000), limits = c(400, 2500)) +
	labs(x = "Unique tree topology (ordered by frequency)", y = "Tree topology frequency") +
	theme(legend.position = "bottom", legend.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave("species_topology_frequencies.pdf", units = "mm", height = 150, width = 150)
ggsave("species_topology_frequencies.png", units = "mm", height = 150, width = 150)
