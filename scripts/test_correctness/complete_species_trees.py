#!/usr/bin/python2.7

import dendropy
import scipy.stats
import random
import numpy

biopy_trees_filename = "biopy.species.trees"
complete_trees_filename = "simulated.species.trees"

species_trees = dendropy.TreeList.get(path = biopy_trees_filename, schema = "nexus", rooting = "force-rooted")

random.seed()
numpy.random.seed()

n_bins = 100
bin_rates = numpy.zeros(n_bins)
for bin_i in range(n_bins):
	bin_percentile = (bin_i + 0.5) / n_bins
	log_bin_rate = scipy.stats.norm.ppf(bin_percentile, loc = -0.0128, scale = 0.16)
	bin_rates[bin_i] = numpy.exp(log_bin_rate)

for species_tree in species_trees:
	n_rates = len(species_tree.nodes()) - 1
	species_tree_categories = numpy.random.randint(100, size = 9)

	node_i = 0
	for node in species_tree.nodes():
		node_category = species_tree_categories[node_i]
		node.annotations["rate"] = bin_rates[node_category]
		node_i += 1

species_trees.write(path = complete_trees_filename, schema = "nexus", suppress_rooting = True, suppress_annotations = False)
