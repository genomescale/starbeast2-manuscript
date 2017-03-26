#!/usr/bin/python2.7

import collections
import csv
import dendropy
import os
import subprocess
import sys

n_reps = 96

internal_node_names = {
	"s1": "s0,s1",
	"s2": "(s0,s1),s2",
	"s3": "((s0,s1),s2),s3"
}

false_slow_counts = {
	"u": collections.defaultdict(int),
	"v": collections.defaultdict(int)
}
false_fast_counts = {
	"u": collections.defaultdict(int),
	"v": collections.defaultdict(int)
}

false_short_counts = {
	"u": collections.defaultdict(int),
	"v": collections.defaultdict(int)
}
false_long_counts = {
	"u": collections.defaultdict(int),
	"v": collections.defaultdict(int)
}

branch_rates_filename = "branch_rates.csv"
branch_rates_file = open(branch_rates_filename, "w")
branch_rates_writer = csv.writer(branch_rates_file)

header_row = ["rep", "analysis", "branch", "rate_mean", "rate_min", "rate_max", "true_length", "estimated_length"]
branch_rates_writer.writerow(header_row)

simulated_tree_filename = "simulated.species.tree"

for rep_i in range(n_reps):
	rep_folder = "rep-%02d" % (rep_i)

	simulated_tree_path = os.path.join(rep_folder, simulated_tree_filename)
	simulated_tree = dendropy.Tree.get(path = simulated_tree_path, schema = "nexus")
	simulated_lengths = {}
	leaf_nodes = simulated_tree.leaf_nodes()
	for leaf_node in leaf_nodes:
		leaf_name = leaf_node.taxon.label
		simulated_lengths[leaf_name] = leaf_node.edge_length
		if leaf_name in internal_node_names:
			internal_node = leaf_node.parent_node
			internal_name = internal_node_names[leaf_name]
			simulated_lengths[internal_name] = internal_node.edge_length

	rep_filenames = os.listdir(rep_folder)
	rep_filenames.sort(reverse = True)

	archive_filenames = [fn for fn in rep_filenames if fn.endswith(".tar.xz")]
	incomplete_filenames = [fn for fn in rep_filenames if fn.endswith(".incomplete")]
	incomplete_run_names = set()
	for incomplete_filename in incomplete_filenames:
		run_name, extension = incomplete_filename.split(".", 1)
		incomplete_run_names.add(run_name)

	treeannotator_bin = sys.argv[1]

	for archive_filename in archive_filenames:
		run_name, group_number_str, extension = archive_filename.split(".", 2)
		alt_name, analysis_code, random_seed = run_name.split("-")

		group_number = int(group_number_str)
		if run_name not in incomplete_run_names:
			incomplete_run_names.add(run_name)
			summary_filename = run_name + ".species_tree"
			summary_path = os.path.join(rep_folder, summary_filename)
			if not os.path.isfile(summary_path):
				posterior_filename = "%s.%02d.trees" % (run_name, group_number)
				posterior_path = os.path.join(rep_folder, posterior_filename)

				extract_cmd = ["tar", "xJvf", archive_filename, posterior_filename]
				subprocess.check_call(extract_cmd, cwd = rep_folder)

				treeannotator_cmd = [treeannotator_bin, "-burnin", "10", posterior_filename, summary_filename]
				subprocess.check_call(treeannotator_cmd, cwd = rep_folder)
				os.remove(posterior_path)

			summary_tree = dendropy.Tree.get(path = summary_path, schema = "nexus")
			leaf_nodes = summary_tree.leaf_nodes()
			for leaf_node in leaf_nodes:
				leaf_name = leaf_node.taxon.label
				rate_mean = float(leaf_node.annotations.get_value("rate"))
				rate_range = leaf_node.annotations.get_value("rate_95%_HPD")
				if rate_range == None:
					rate_min = rate_mean
					rate_max = rate_mean
				else:
					rate_min = float(rate_range[0])
					rate_max = float(rate_range[1])

				true_length = simulated_lengths[leaf_name]
				length_range = leaf_node.annotations.get_value("length_95%_HPD")
				length_min = float(length_range[0])
				length_max = float(length_range[1])

				output_row = [rep_i, analysis_code, leaf_name, rate_mean, rate_min, rate_max, true_length, leaf_node.edge_length]
				branch_rates_writer.writerow(output_row)

				if rate_max < 1.0: false_slow_counts[analysis_code][leaf_name] += 1
				elif rate_min > 1.0: false_fast_counts[analysis_code][leaf_name] += 1
				if length_max < true_length: false_short_counts[analysis_code][leaf_name] += 1
				elif length_min > true_length: false_long_counts[analysis_code][leaf_name] += 1

				if leaf_name in internal_node_names:
					internal_node = leaf_node.parent_node
					internal_name = internal_node_names[leaf_name]
					rate_mean = float(internal_node.annotations.get_value("rate"))
					rate_range = internal_node.annotations.get_value("rate_95%_HPD")
					if rate_range == None:
						rate_min = rate_mean
						rate_max = rate_mean
					else:
						rate_min = float(rate_range[0])
						rate_max = float(rate_range[1])

					true_length = simulated_lengths[internal_name]
					length_range = internal_node.annotations.get_value("length_95%_HPD")
					length_min = float(length_range[0])
					length_max = float(length_range[1])

					output_row = [rep_i, analysis_code, internal_name, rate_mean, rate_min, rate_max, true_length, internal_node.edge_length]
					branch_rates_writer.writerow(output_row)

					if rate_max < 1.0: false_slow_counts[analysis_code][internal_name] += 1
					elif rate_min > 1.0: false_fast_counts[analysis_code][internal_name] += 1
					if length_max < true_length: false_short_counts[analysis_code][internal_name] += 1
					elif length_min > true_length: false_long_counts[analysis_code][internal_name] += 1

branch_rates_file.close()

print(false_slow_counts)
print(false_fast_counts)

print(false_short_counts)
print(false_long_counts)
