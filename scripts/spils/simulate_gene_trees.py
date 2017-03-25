#!/usr/bin/python2.7

import collections
import dendropy
import numpy
import os
import random
import scipy.stats
import subprocess
import sys

n_species_trees = int(sys.argv[1])
n_gene_trees = int(sys.argv[2])

species_tree_filename = "simulated.species.tree"
biopy_trees_filename = "biopy.gene.trees"
complete_trees_filename = "simulated.gene.trees"

biopy_cmd = ["genetree_in_sptree_sim", "-n", str(n_gene_trees), "-t", "1", "-o", biopy_trees_filename, species_tree_filename]

random.seed()
numpy.random.seed()

for rep_i in range(n_species_trees):
	rep_folder = "rep-%02d" % (rep_i)
	species_tree_path     = os.path.join(rep_folder, species_tree_filename)
	biopy_trees_path      = os.path.join(rep_folder, biopy_trees_filename)
	complete_trees_path   = os.path.join(rep_folder, complete_trees_filename)

	species_tree = dendropy.Tree.get(path = species_tree_path, schema = "nexus", rooting = "force-rooted", preserve_underscores = True)
	species_tree.calc_node_ages()

	subprocess.check_call(biopy_cmd, cwd = rep_folder)
	gene_trees = dendropy.TreeList.get(path = biopy_trees_path, schema = "nexus", rooting = "force-rooted", preserve_underscores = True)
	gene_trees.write(path = complete_trees_path, schema = "nexus", suppress_rooting = True, suppress_annotations = False)
