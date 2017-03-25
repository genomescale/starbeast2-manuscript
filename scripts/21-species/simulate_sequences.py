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
n_bases = int(sys.argv[2])

gene_trees_filename = "simulated.gene.trees"
sequences_filename = "simulated.sequences"
seqgen_cmd = ["seq-gen", "-l" + str(n_bases), "-mHKY", "-a0.2", "-g4", "-t1.5", gene_trees_filename]

for rep_i in range(n_species_trees):
	rep_folder = "rep-%02d" % (rep_i)
	sequences_path = os.path.join(rep_folder, sequences_filename)

	sequences_file = open(sequences_path, "w")
	simulated_sequences = subprocess.check_output(seqgen_cmd, cwd = rep_folder)
	sequences_file.write(simulated_sequences)
	sequences_file.close()
