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
n_tip_nodes = 21
n_internal_nodes = n_tip_nodes - 2

cred_interval = 0.95

simulation_folder = "21-species"

branch_length_accuracy_filename = "branch_length_accuracy.csv"
branch_length_accuracy_path = os.path.join(simulation_folder, branch_length_accuracy_filename)
branch_length_accuracy_file = open(branch_length_accuracy_path, "w")
branch_length_accuracy_writer = csv.writer(branch_length_accuracy_file)

output_header = ["method", "n_loci", "clock", "sequence", "tip_length_error", "tip_length_bias", "tip_truth_in_hpd", "internal_length_error", "internal_length_bias", "internal_truth_in_hpd"]
branch_length_accuracy_writer.writerow(output_header)

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
				truth_length = truth_node.edge_length
				summary_length = summary_node.edge_length
				summary_hpd = summary_node.annotations["length_95%_HPD"]
				if summary_hpd.value == "":
					summary_hpd_low = summary_length
					summary_hpd_high = summary_length
				else:
					summary_hpd_low = float(summary_hpd.value[0])
					summary_hpd_high = float(summary_hpd.value[1])

				truth_in_hpd = (truth_length >= summary_hpd_low) and (truth_length <= summary_hpd_high)

				if truth_node.is_leaf():
					if truth_in_hpd: tip_truth_in_hpd += 1.0 / float(n_tip_nodes)
					tip_sum_truth += truth_length
					tip_sum_estimate += summary_length
					tip_sum_error += abs(summary_length - truth_length)
				else:
					if truth_in_hpd: internal_truth_in_hpd += 1.0 / float(n_internal_nodes)
					internal_sum_truth += truth_length
					internal_sum_estimate += summary_length
					internal_sum_error += abs(summary_length - truth_length)

		tip_length_error = tip_sum_error / tip_sum_truth
		tip_length_bias = (tip_sum_estimate - tip_sum_truth) / tip_sum_truth

		internal_length_error = internal_sum_error / internal_sum_truth
		internal_length_bias = (internal_sum_estimate - internal_sum_truth) / internal_sum_truth

		output_row = [method_name, n_loci, clock, sequence, tip_length_error, tip_length_bias, tip_truth_in_hpd, internal_length_error, internal_length_bias, internal_truth_in_hpd]
		branch_length_accuracy_writer.writerow(output_row)

branch_length_accuracy_file.close()
