#!/usr/bin/python2.7

import collections
import dendropy
import numpy
import os
import random
import scipy.stats
import subprocess
import sys

clock_rates_stdev = 0.6

n_species_trees = int(sys.argv[1])
n_gene_trees = int(sys.argv[2])

def recurse_coalescence(gene_node, species_node, species_tree, substitution_lengths, assigned_nodes):
	species_branch_rate = float(species_node.annotations["rate"].value)
	floor = max(gene_node.age, species_node.age)
	stay_in_current_branch = False
	if species_node == species_tree.seed_node:
		stay_in_current_branch = True
	elif gene_node.parent_node.age < species_node.parent_node.age:
		stay_in_current_branch = True

	if stay_in_current_branch:
		ceiling = gene_node.parent_node.age
		if gene_node.parent_node not in assigned_nodes:
			assigned_nodes.add(gene_node.parent_node)
			recurse_coalescence(gene_node.parent_node, species_node, species_tree, substitution_lengths, assigned_nodes)
	else:
		ceiling = species_node.parent_node.age
		recurse_coalescence(gene_node, species_node.parent_node, species_tree, substitution_lengths, assigned_nodes)

	substitution_lengths[gene_node] += species_branch_rate * (ceiling - floor)

species_tree_filename = "simulated.species.tree"
biopy_trees_filename = "biopy.gene.trees"
complete_trees_filename = "simulated.gene.trees"
clock_rates_filename = "clock_rates.txt"

biopy_cmd = ["genetree_in_sptree_sim", "-n", str(n_gene_trees), "-t", "2", "-o", biopy_trees_filename, species_tree_filename]

random.seed()
numpy.random.seed()

for rep_i in range(n_species_trees):
	rep_folder = "rep-%02d" % (rep_i)
	species_tree_path     = os.path.join(rep_folder, species_tree_filename)
	biopy_trees_path      = os.path.join(rep_folder, biopy_trees_filename)
	complete_trees_path   = os.path.join(rep_folder, complete_trees_filename)
	clock_rates_path      = os.path.join(rep_folder, clock_rates_filename)

	species_tree = dendropy.Tree.get(path = species_tree_path, schema = "nexus", rooting = "force-rooted", preserve_underscores = True)
	species_tree.calc_node_ages()

	subprocess.check_call(biopy_cmd, cwd = rep_folder)
	gene_trees = dendropy.TreeList.get(path = biopy_trees_path, schema = "nexus", rooting = "force-rooted", preserve_underscores = True)

	clock_rates = None
	clock_rates_file = open(clock_rates_path, "w")
	for gene_tree_i, gene_tree in enumerate(gene_trees):
		gene_tree.calc_node_ages()

		clock_rates = scipy.stats.lognorm.rvs(s = clock_rates_stdev, scale = numpy.exp(-clock_rates_stdev * clock_rates_stdev * 0.5), size = n_gene_trees)

		overall_clock_rate = clock_rates[gene_tree_i]
		clock_rates_file.write("%s,%f\n" % (gene_tree.label, overall_clock_rate))

		substitution_lengths = collections.defaultdict(float)
		assigned_nodes = set([gene_tree.seed_node])
		for gene_leaf in gene_tree.leaf_nodes():
			species_name = gene_leaf.taxon.label.split("_")[0]
			for species_leaf in species_tree.leaf_nodes():
				if species_leaf.taxon.label == species_name:
					recurse_coalescence(gene_leaf, species_leaf, species_tree, substitution_lengths, assigned_nodes)

		for gene_node_i, gene_node in enumerate(gene_tree):
			if gene_node != gene_tree.seed_node:
				gene_node.edge_length = substitution_lengths[gene_node] * overall_clock_rate

	clock_rates_file.close()
	gene_trees.write(path = complete_trees_path, schema = "nexus", suppress_rooting = True, suppress_annotations = False)
