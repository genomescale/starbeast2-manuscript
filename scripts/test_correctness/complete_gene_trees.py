#!/usr/bin/python2.7

import collections
import dendropy
import sys

gene_tree_number = int(sys.argv[1])
overall_clock_rate = float(sys.argv[2])

def recurse_coalescence(gene_node, species_node, species_tree, unnormalized_rates, assigned_nodes):
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
			recurse_coalescence(gene_node.parent_node, species_node, species_tree, unnormalized_rates, assigned_nodes)
	else:
		ceiling = species_node.parent_node.age
		recurse_coalescence(gene_node, species_node.parent_node, species_tree, unnormalized_rates, assigned_nodes)

	unnormalized_rates[gene_node] += species_branch_rate * (ceiling - floor)

species_trees_filename = "simulated.species.trees"
gene_trees_filename = "biopy.gene%d.trees" % (gene_tree_number)
complete_trees_filename = "simulated.gene%d.trees" % (gene_tree_number)

gene_trees = dendropy.TreeList.get(path = gene_trees_filename, schema = "nexus", rooting = "force-rooted")
species_trees = dendropy.TreeList.get(path = species_trees_filename, schema = "nexus", rooting = "force-rooted")

for tree_i, gene_tree in enumerate(gene_trees):
	species_tree = species_trees[tree_i]

	gene_tree.calc_node_ages()
	species_tree.calc_node_ages()

	unnormalized_rates = collections.defaultdict(float)
	assigned_nodes = set([gene_tree.seed_node])
	for gene_leaf in gene_tree.leaf_nodes():
		for species_leaf in species_tree.leaf_nodes():
			species_name = gene_leaf.taxon.label[:8]
			if species_name == species_leaf.taxon.label:
				recurse_coalescence(gene_leaf, species_leaf, species_tree, unnormalized_rates, assigned_nodes)

	for gene_node_i, gene_node in enumerate(gene_tree):
		if gene_node == gene_tree.seed_node:
			gene_branch_rate = overall_clock_rate
		else:
			try:
				gene_branch_rate = overall_clock_rate * unnormalized_rates[gene_node] / gene_node.edge_length
			except Exception:
				print(tree_i, gene_node_i, overall_clock_rate, unnormalized_rates[gene_node], gene_node.edge_length)

		gene_node.annotations.add_new("rate", gene_branch_rate)

gene_trees.write(path = complete_trees_filename, schema = "nexus", suppress_rooting = True, suppress_annotations = False)
