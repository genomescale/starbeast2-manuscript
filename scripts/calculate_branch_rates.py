#!/usr/bin/python2.7

import csv
import numpy
import os
import dendropy

def tree_to_dict(tree):
	tree_dict = {}
	for n in tree.postorder_node_iter():
		tree_dict[n.leafset_bitmask] = n

	return tree_dict

n_replicates = 30

simulation_folder = "21-species"

branch_rates_filename = "branch_rates.csv"
branch_rates_path = os.path.join(simulation_folder, branch_rates_filename)
branch_rates_file = open(branch_rates_path, "w")
branch_rates_writer = csv.writer(branch_rates_file)

output_header = ["method", "n_loci", "clock", "sequence", "rate", "estimated_rate", "estimated_rate_low", "estimated_rate_high"]
branch_rates_writer.writerow(output_header)

truth_filename = "simulated.species.tree"

for rep_i in range(n_replicates):
	rep_subfolder = "rep-%02d" % (rep_i)
	rep_folder = os.path.join(simulation_folder, rep_subfolder)

	truth_path = os.path.join(rep_folder, truth_filename)
	truth_tree = dendropy.Tree.get(path = truth_path, schema = "nexus", rooting = "force-rooted")
	truth_tree.encode_bipartitions()
	truth_tree.calc_node_ages()
	truth_dict = tree_to_dict(truth_tree)

	rep_filenames = os.listdir(rep_folder)
	summary_filenames = [fn for fn in rep_filenames if fn.endswith("truth_heights.tree")]
	for summary_filename in summary_filenames:
		analysis_name, random_seed, extension = summary_filename.split(".", 2)
		method_name, n_loci, config = analysis_name.split("-")
		if config[0] == "0": clock = "strict"
		elif config[0] == "1": clock = "gt-ucln"
		else: clock = "st-ucln"

		if method_name == "beast" and config[1] == "1": sequence = "diplotype"
		else: sequence = "haplotype"

		summary_path = os.path.join(rep_folder, summary_filename)
		summary_tree = dendropy.Tree.get(path = summary_path, schema = "nexus", rooting = "force-rooted", taxon_namespace = truth_tree.taxon_namespace)
		summary_tree.encode_bipartitions()
		summary_tree.calc_node_ages()
		summary_dict = tree_to_dict(summary_tree)

		tip_sum_truth = 0.0
		tip_sum_estimate = 0.0
		tip_sum_error = 0.0

		internal_sum_truth = 0.0
		internal_sum_estimate = 0.0
		internal_sum_error = 0.0

		tip_truth_in_hpd = 0.0
		internal_truth_in_hpd = 0.0
		for bitmask in sorted(truth_dict):
			truth_node = truth_dict[bitmask]
			summary_node = summary_dict[bitmask]
			if truth_node != truth_tree.seed_node: # not the root
				truth_rate_annotation = truth_node.annotations["rate"]
				truth_rate = float(truth_rate_annotation.value)

				summary_rate_annotation = summary_node.annotations["rate"]
				if summary_rate_annotation.value == "":
					summary_rate = -1.0
				else:
					summary_rate = float(summary_rate_annotation.value)

				summary_hpd_annotation = summary_node.annotations["rate_95%_HPD"]
				if summary_hpd_annotation.value == "":
					summary_hpd_low = summary_rate
					summary_hpd_high = summary_rate
				else:
					summary_hpd_low = float(summary_hpd_annotation.value[0])
					summary_hpd_high = float(summary_hpd_annotation.value[1])

				output_row = [method_name, n_loci, clock, sequence, truth_rate, summary_rate, summary_hpd_low, summary_hpd_high]
				branch_rates_writer.writerow(output_row)

branch_rates_file.close()
