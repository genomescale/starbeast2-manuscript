#!/usr/bin/python2.7

import collections
import csv
import numpy
import os

def bootstrap(values, n_reps = 2000, ci = 0.95):
	ci_low = (100.0 - (ci * 100.0)) * 0.5
	ci_high = 100.0 - ci_low

	n_values = len(values)
	sample_mean = numpy.mean(values)

	bootstrap_means = numpy.zeros((n_reps))
	for rep_i in range(n_reps):
		bootstrap_values = numpy.random.choice(values, size = n_values, replace = True)
		bootstrap_means[rep_i] = numpy.mean(bootstrap_values)

	x = numpy.percentile(bootstrap_means, ci_low)
	y = numpy.percentile(bootstrap_means, ci_high)

	return (sample_mean, x, y)

n_replicates = 30

cred_interval = 0.95

simulation_folder = "21-species"

method_in_credible_set = collections.defaultdict(list)
method_mcc_rf_distance = collections.defaultdict(list)

for rep_i in range(n_replicates):
	rep_subfolder = "rep-%02d" % (rep_i)
	rep_folder = os.path.join(simulation_folder, rep_subfolder)
	rep_filenames = os.listdir(rep_folder)
	topfq_filenames = [fn for fn in rep_filenames if fn.endswith("topfq.csv")]
	for topfq_filename in topfq_filenames:
		analysis_name, random_seed, extension = topfq_filename.split(".", 2)
		topfq_path = os.path.join(rep_folder, topfq_filename)
		print(topfq_path)

		topfq_file = open(topfq_path)
		topfq_reader = csv.reader(topfq_file)
		mcc_rf_distance = 0
		cumulative_probability = 0.0
		in_credible_set = 0
		for row_i, row in enumerate(topfq_reader):
			if row_i >= 1:
				topology_probability = float(row[1])
				topfq_distance = int(row[2])
				mcc_distance = int(row[3])

				if mcc_distance == 0:
					mcc_rf_distance = topfq_distance

				if topfq_distance == 0 and cumulative_probability <= cred_interval:
					in_credible_set = 1

				cumulative_probability += topology_probability

		topfq_file.close()

		method_in_credible_set[analysis_name].append(in_credible_set)
		method_mcc_rf_distance[analysis_name].append(mcc_rf_distance)

topology_accuracy_filename = "topology_accuracy.csv"
topology_accuracy_path = os.path.join(simulation_folder, topology_accuracy_filename)
topology_accuracy_file = open(topology_accuracy_path, "w")
topology_accuracy_writer = csv.writer(topology_accuracy_file)

output_header = ["method", "n_loci", "clock", "sequence", "mcc_rf_mean", "mcc_rf_low", "mcc_rf_high", "in_cs_prop", "in_cs_low", "in_cs_high"]
topology_accuracy_writer.writerow(output_header)

for analysis_name in sorted(method_in_credible_set):
	method_name, n_loci, config = analysis_name.split("-")
	if config[0] == "0": clock = "strict"
	elif config[0] == "1": clock = "gt-ucln"
	else: clock = "st-ucln"

	if method_name == "beast" and config[1] == "1": sequence = "diplotype"
	else: sequence = "haplotype"

	in_credible_set_values = method_in_credible_set[analysis_name]
	mcc_rf_distance_values = method_mcc_rf_distance[analysis_name]
	output_row = [method_name, n_loci, clock, sequence]
	output_row.extend(bootstrap(mcc_rf_distance_values))
	output_row.extend(bootstrap(in_credible_set_values))

	topology_accuracy_writer.writerow(output_row)

topology_accuracy_file.close()
