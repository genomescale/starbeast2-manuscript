#!/usr/bin/python2.7

import dendropy
import numpy
import sys

mean_popsizes = []
mean_lengths = []
mean_coalescent_lengths = []
mean_heights = []

n_burnin = 0

if sys.argv[-1].isdigit():
	n_burnin = int(sys.argv[-1])
	tree_paths = sys.argv[1:-1]
else:
	tree_paths = sys.argv[1:]

for tree_path in tree_paths:
	all_popsizes = []
	all_lengths = []
	all_coalescent_lengths = []

	if tree_path.endswith("tree"):
		tree = dendropy.Tree.get(path = tree_path, schema = "nexus")
		mean_heights.append(tree.seed_node.distance_from_tip())
		for tree_node in tree.nodes():
			if tree_node != tree.seed_node:
				dmv = tree_node.annotations["dmv"].value
				diploid_popsize = float(dmv) * 0.5
				branch_length = tree_node.edge_length
				coalescent_length = (branch_length / (2.0 * diploid_popsize))

				all_popsizes.append(diploid_popsize)
				all_lengths.append(branch_length)
				all_coalescent_lengths.append(coalescent_length)
	else:
		all_trees = dendropy.TreeList.get(path = tree_path, schema = "nexus")
		for tree in all_trees[n_burnin:]:
			mean_heights.append(tree.seed_node.distance_from_tip())
			for tree_node in tree.nodes():
				if tree_node != tree.seed_node:
					branch_length = tree_node.edge_length
					diploid_popsize = 1.0
					try:
						dmv = tree_node.annotations["dmv"].value[0]
						diploid_popsize = float(dmv)
					except Exception:
						try:
							theta = tree_node.annotations["theta"].value
							diploid_popsize = float(theta) * 0.25
						except Exception:
							print("No DMV found...")

					coalescent_length = (branch_length / (2.0 * diploid_popsize))

					all_popsizes.append(diploid_popsize)
					all_lengths.append(branch_length)
					all_coalescent_lengths.append(coalescent_length)

		print(numpy.mean(all_popsizes), numpy.mean(all_lengths), numpy.mean(all_coalescent_lengths))

	mean_popsizes.append(numpy.mean(all_popsizes))
	mean_lengths.append(numpy.mean(all_lengths))
	mean_coalescent_lengths.append(numpy.mean(all_coalescent_lengths))

mean_popsize = numpy.mean(mean_popsizes)
mean_length = numpy.mean(mean_lengths)
mean_coalescent_length = numpy.mean(mean_coalescent_lengths)
mean_height = numpy.mean(mean_heights)

print(mean_popsize, mean_length, mean_coalescent_length, mean_height)
