#!/usr/bin/python2.7

import dendropy
from dendropy.calculate import treecompare
import numpy
import os
import random
import scipy.stats
import subprocess
import sys

n_species_trees = int(sys.argv[1])

catapillar_newick = "((((s0,s1),s2),s3),s4);"
catapillar_tree = dendropy.Tree.get(data = catapillar_newick, schema = "newick", rooting = "force-rooted")

biopy_tree_filename = "biopy.species.tree"
complete_tree_filename = "simulated.species.tree"

biopy_cmd = ["sample_species_tree", "-o", biopy_tree_filename, "-p", "i,3,0.2", "-b", "10.0", "-d", "0.0", "-n", "1000", "5"]

random.seed()
numpy.random.seed()

for rep_i in range(n_species_trees):
	rep_folder = "rep-%02d" % (rep_i)
	if not os.path.exists(rep_folder):
		os.mkdir(rep_folder)

	biopy_tree_path = os.path.join(rep_folder, biopy_tree_filename)

	subprocess.check_call(biopy_cmd, cwd = rep_folder)
	biopy_trees = dendropy.TreeList.get(path = biopy_tree_path, schema = "nexus", rooting = "force-rooted", taxon_namespace = catapillar_tree.taxon_namespace)

	is_catapillar = False
	for bt in biopy_trees:
		print(bt)
		print(catapillar_tree)
		print(treecompare.symmetric_difference(catapillar_tree, bt))
		if treecompare.symmetric_difference(catapillar_tree, bt) == 0:
			species_tree = bt

	complete_tree_path = os.path.join(rep_folder, complete_tree_filename)
	species_tree.write(path = complete_tree_path, schema = "nexus", suppress_rooting = True, suppress_annotations = False)
