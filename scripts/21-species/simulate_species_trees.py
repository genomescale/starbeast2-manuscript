#!/usr/bin/python2.7

import dendropy
import numpy
import os
import random
import scipy.stats
import subprocess
import sys

branch_rates_stdev = 0.3

n_species = int(sys.argv[1])
n_species_trees = int(sys.argv[2])

biopy_tree_filename = "biopy.species.tree"
complete_tree_filename = "simulated.species.tree"

biopy_cmd = ["sample_species_tree", "-o", biopy_tree_filename, "-p", "g,2,0.002", "-b", "100.0", "-d", "30.0", "-n", "1", str(n_species)]

random.seed()
numpy.random.seed()

for rep_i in range(n_species_trees):
	rep_folder = "rep-%02d" % (rep_i)
	if not os.path.exists(rep_folder):
		os.mkdir(rep_folder)

	biopy_tree_path = os.path.join(rep_folder, biopy_tree_filename)
	complete_tree_path = os.path.join(rep_folder, complete_tree_filename)

	subprocess.check_call(biopy_cmd, cwd = rep_folder)
	species_tree = dendropy.Tree.get(path = biopy_tree_path, schema = "nexus", rooting = "force-rooted")

	n_branches = len(species_tree.nodes()) - 1
	branch_rates = scipy.stats.lognorm.rvs(s = branch_rates_stdev, scale = numpy.exp(-branch_rates_stdev * branch_rates_stdev * 0.5), size = n_branches)
	norm_scale = 1.0 / numpy.mean(branch_rates)

	total_rate = 0.0
	node_i = 0
	for node in species_tree.seed_node.postorder_iter():
		if node == species_tree.seed_node:
			node.annotations["rate"] = 1.0
		else:
			node.annotations["rate"] = branch_rates[node_i] * norm_scale
			total_rate += branch_rates[node_i] * norm_scale
			node_i += 1

	species_tree.write(path = complete_tree_path, schema = "nexus", suppress_rooting = True, suppress_annotations = False)
