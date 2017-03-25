#!/usr/bin/python2.7

import os
import dendropy
import sys

seq_folder = sys.argv[1]

for input_filename in os.listdir(seq_folder):
	input_path = os.path.join(seq_folder, input_filename)
	base_name, file_ext = os.path.splitext(input_filename)
	if file_ext == ".nex" or file_ext == ".nexus":
		output_filename = base_name + ".fasta"
		output_path = os.path.join(seq_folder, output_filename)

		sequence_data = dendropy.DnaCharacterMatrix.get(path = input_path, schema = "nexus")
		sequence_data.write(path = output_path, schema = "fasta")
