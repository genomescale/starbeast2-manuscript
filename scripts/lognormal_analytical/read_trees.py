#!/usr/bin/python2.7

import base64
import csv
import dendropy
import struct
import sys

trees_path = sys.argv[1]

if trees_path.find("test") == -1:
	offset = 0
else:
	offset = 1

species_namespace = dendropy.TaxonNamespace()
gene_namespace = dendropy.TaxonNamespace()
for species_i in range(5):
	species_taxon = dendropy.Taxon("species%d" % (species_i))
	species_namespace.add_taxon(species_taxon)
	gene_taxon = dendropy.Taxon("species%d_tip0" % (species_i))
	gene_namespace.add_taxon(gene_taxon)

if "species" in trees_path:
	trees = dendropy.TreeList.get(path = trees_path, schema = "nexus", rooting = "force-rooted", tree_offset = offset, taxon_namespace = species_namespace, preserve_underscores = True)
	n_internal_nodes = 3
else:
	trees = dendropy.TreeList.get(path = trees_path, schema = "nexus", rooting = "force-rooted", tree_offset = offset, taxon_namespace = gene_namespace, preserve_underscores = True)
	n_internal_nodes = 3

tree_properties_path = trees_path + ".topologies.csv"
node_properties_path = trees_path + ".nodes.csv"

tree_properties_file = open(tree_properties_path, "w")
node_properties_file = open(node_properties_path, "w")

tree_properties_writer = csv.writer(tree_properties_file)
node_properties_writer = csv.writer(node_properties_file)

tree_properties_header = ["tree_id", "tree_topology", "topology_code"]
node_properties_header = ["tree_id", "node_rank", "node_height", "branch_rate", "pop_size"]

tree_properties_writer.writerow(tree_properties_header)
node_properties_writer.writerow(node_properties_header)

code_to_newick = {}
# all_topology_codes = []
tree_i = 0
for tree in trees:
	internal_node_codes = []
	tree.encode_bipartitions()
	node_i = 0
	for node in tree.ageorder_node_iter():
		node_height = node.age

		if not (node.is_leaf() or node == tree.seed_node):
			internal_node_codes.append(node.bipartition.leafset_bitmask)

		rate = node.annotations.find(name = "rate")
		if rate == None:
			branch_rate = "NA"
		else:
			branch_rate = float(rate.value)

		dmv = node.annotations.find(name = "dmv")
		if dmv == None:
			pop_size = "NA"
		elif type(dmv.value) == list:
			pop_size = float(dmv.value[0])
		else:
			pop_size = float(dmv.value)

		node_properties_writer.writerow([tree_i, node_i, node_height, branch_rate, pop_size])
		node_i += 1

	internal_node_codes.sort()
	tree_topology_code = struct.pack("B" * n_internal_nodes, *internal_node_codes)
	# all_topology_codes.append(tree_topology_code)

	if tree_topology_code not in code_to_newick:
		tree_topology_newick = tree.as_string(schema = "newick", suppress_edge_lengths = True, suppress_rooting = True)[:-2]
		code_to_newick[tree_topology_code] = tree_topology_newick
	else:
		tree_topology_newick = code_to_newick[tree_topology_code]

	base64_code = base64.standard_b64encode(tree_topology_code)
	tree_properties_writer.writerow([tree_i, tree_topology_newick, base64_code])
	tree_i += 1

tree_properties_file.close()
node_properties_file.close()

# tree_i = 0
# for tree in trees:
# 	tree_topology_code = all_topology_codes[tree_i]

# 	for comparison_code, comparison_newick in code_to_newick.items():
# 		comparison_tree = dendropy.Tree.get(data = comparison_newick + ";\n", schema = "newick", rooting = "force-rooted", taxon_namespace = trees.taxon_namespace)
# 		symmetric_difference = dendropy.calculate.treecompare.symmetric_difference(tree, comparison_tree)
		
# 		if tree_topology_code == comparison_code:
# 			assert symmetric_difference == 0
# 		else:
# 			assert symmetric_difference > 0

# 	tree_topology = code_to_newick[tree_topology_code]
# 	tree_i += 1
