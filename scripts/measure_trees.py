#!/usr/bin/python2.7

import csv
import dendropy
import dendropy.calculate
import os
import subprocess
import sys

n_burnin = 1025

rep_i = int(sys.argv[1])
simulation_folder = "21-species"
rep_subfolder = "rep-%02d" % (rep_i)
rep_folder = os.path.join(simulation_folder, rep_subfolder)

reference_filename = "simulated.species.tree"
reference_path =  os.path.join(rep_folder, reference_filename)

reference_tree = dendropy.Tree.get(path = reference_path, schema = "nexus", rooting = "force-rooted")
reference_tree.encode_bipartitions()
reference_tree.calc_node_ages()

rep_filenames = os.listdir(rep_folder)
rep_filenames.sort(reverse = True)

archive_filenames = [fn for fn in rep_filenames if fn.endswith(".tar.xz")]
incomplete_run_names = set()

# Commented out - just measure all trees (ESSes are good enough)
# incomplete_filenames = [fn for fn in rep_filenames if fn.endswith(".incomplete")]
# for incomplete_filename in incomplete_filenames:
# 	run_name, extension = incomplete_filename.split(".", 1)
# 	incomplete_run_names.add(run_name)

output_header = ["tree_index", "taxa", "subclade_a", "subclade_b", "dmv", "rate", "height", "length"]
topfq_header = ["newick", "count", "truth_distance", "mcc_distance"]

for archive_filename in archive_filenames:
	run_name, random_seed, group_number_str, extension = archive_filename.split(".", 3)
	alt_name, n_loci, analysis_code = run_name.split("-")

	if run_name not in incomplete_run_names and (alt_name != "starbeast2" or analysis_code.endswith("3")):
		incomplete_run_names.add(run_name)
		group_number = int(group_number_str)

		truth_filename = "%s.%s.truth.csv" % (run_name, random_seed)
		beast_filename = "%s.%s.beast.csv" % (run_name, random_seed)
		topfq_filename = "%s.%s.topfq.csv" % (run_name, random_seed)
		truth_path = os.path.join(rep_folder, truth_filename)
		beast_path = os.path.join(rep_folder, beast_filename)
		topfq_path = os.path.join(rep_folder, topfq_filename)
		if not (os.path.isfile(truth_path) and os.path.isfile(beast_path)):
			truth_file = open(truth_path, "w")
			truth_writer = csv.writer(truth_file)
			truth_writer.writerow(output_header)

			for node in reference_tree.nodes():
				taxa = node.bipartition.leafset_bitmask
				dmv = float(node.annotations["dmv"].value)
				rate = float(node.annotations["rate"].value)
				height = float(node.age)
				if node == reference_tree.seed_node:
					length = 0.0
				else:
					length = node.edge_length

				if node.is_leaf():
					subclade_a = 0
					subclade_b = 0
				else:
					subnode_a, subnode_b = node.child_nodes()
					subclade_a = subnode_a.bipartition.leafset_bitmask
					subclade_b = subnode_b.bipartition.leafset_bitmask

				output_row = [rep_i, taxa, subclade_a, subclade_b, dmv, rate, height, length]
				truth_writer.writerow(output_row)

			truth_file.close()

			posterior_filename = "%s.%s.%02d.trees" % (run_name, random_seed, group_number)
			extract_cmd = ["tar", "xJvf", archive_filename, posterior_filename]
			subprocess.check_call(extract_cmd, cwd = rep_folder)
			posterior_path = os.path.join(rep_folder, posterior_filename)
			posterior_trees = dendropy.TreeList.get(path = posterior_path, schema = "nexus", rooting = "force-rooted", taxon_namespace = reference_tree.taxon_namespace, tree_offset = n_burnin)
			os.remove(posterior_path)

			beast_file = open(beast_path, "w")
			beast_writer = csv.writer(beast_file)
			beast_writer.writerow(output_header)

			for tree_i, tree in enumerate(posterior_trees):
				tree.encode_bipartitions()
				tree.calc_node_ages()

				for node in tree.nodes():
					taxa = node.bipartition.leafset_bitmask

					dmv_str = node.annotations["dmv"].value
					if dmv_str == "":
						dmv = 0.0
					else:
						dmv = float(dmv_str[0])

					rate_str = node.annotations["rate"].value
					if rate_str == "":
						rate = 0.0
					else:
						rate = float(rate_str)

					height = float(node.age)

					if node == tree.seed_node:
						length = 0.0
					else:
						length = node.edge_length

					if node.is_leaf():
						subclade_a = 0
						subclade_b = 0
					else:
						subnode_a, subnode_b = node.child_nodes()
						subclade_a = subnode_a.bipartition.leafset_bitmask
						subclade_b = subnode_b.bipartition.leafset_bitmask

					output_row = [tree_i, taxa, subclade_a, subclade_b, dmv, rate, height, length]
					beast_writer.writerow(output_row)

			beast_file.close()

			mcc_tree = posterior_trees.maximum_product_of_split_support_tree()
			mcc_tree.encode_bipartitions()

			topfq_file = open(topfq_path, "w")
			topfq_writer = csv.writer(topfq_file)
			topfq_writer.writerow(topfq_header)

			posterior_treearray = dendropy.TreeArray(taxon_namespace = reference_tree.taxon_namespace, is_rooted_trees = True)
			posterior_treearray.add_trees(posterior_trees, is_bipartitions_updated = True)

			beast_topologies = posterior_treearray.topologies(sort_descending = True)
			is_mcc = 0
			for topology in beast_topologies:
				topology.encode_bipartitions()
				topology_newick = topology.as_string(schema = "newick", suppress_rooting = True).strip()

				truth_rf_distance = dendropy.calculate.treecompare.symmetric_difference(reference_tree, topology, is_bipartitions_updated = True)
				mcc_rf_distance = dendropy.calculate.treecompare.symmetric_difference(mcc_tree, topology, is_bipartitions_updated = True) 

				topfq_output = [topology_newick, topology.frequency, truth_rf_distance, mcc_rf_distance]
				topfq_writer.writerow(topfq_output)

			topfq_file.close()
