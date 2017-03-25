#!/usr/bin/python2.7

import csv
import random
import numpy
import collections
import os
import string
import xml.etree.ElementTree
import sys

crocidura = [
	("beatus", ["FMNH166459", "KU167037", "KU167039", "KU165969"]),
	("mindorus", ["FMNH221890"]),
	("palawanensis", ["FMNH195992", "FMNH195991"]),
	("grayi", ["FMNH218425", "KU165178", "KU165176", "KU165912", "FMNH186719"]),
	("negrina", ["KU165049", "KU165103"]),
	("panayensis", ["KU164878", "KU164877"]),
	("ninoyi", ["FMNH145685"]),
	("sp", ["FMNH146788"]),
	("orientalis", ["FMNH212778"])
]

def read_fasta(fasta_path):
	fasta_data = {}
	fasta_file = open(fasta_path)
	label = ""
	sequence = ""
	l = fasta_file.readline()
	while l != "":
		if l[0] == ">":
			if label != "":
				species, individual = label.split()
				fasta_data[individual] = sequence
			label = l[1:].strip()
			sequence = ""
		else:
			sequence += l.strip().upper()

		l = fasta_file.readline()

	if label != "":
		species, individual = label.split()
		fasta_data[individual] = sequence

	fasta_file.close()

	return fasta_data

hky_models = {'F81', 'HKY', 'JC69', 'K2P'}

substitution_models_filename = "site_models.csv"
substitution_models_file = open(substitution_models_filename)
substitution_models_file.readline() # skip header
substitution_models_reader = csv.reader(substitution_models_file)

uce_folder = "new-incomplete-50percent-all-species-padded"

hky_loci = {}
for row in substitution_models_reader:
	uce_name = row[0].lower()
	substitution_model = row[5].split("+")[0]

	if substitution_model in hky_models:
		uce_filename = uce_name + ".fasta"
		uce_path = os.path.join(uce_folder, uce_filename)
		hky_loci[uce_name] = read_fasta(uce_path)

pbs_template = """#!/bin/bash
#PBS -P PROJECT_NAME
#PBS -q normal
#PBS -l walltime=%(walltime)s
#PBS -l mem=8GB
#PBS -l ncpus=1
#PBS -l wd

module load java

if [ -f %(incomplete_filename)s ]
	then
		python ../run_raijin.py %(analysis)s %(seed)d 0 %(seconds)s
fi

if [ -f %(incomplete_filename)s ]
	then
		qsub %(analysis)s.%(seed)010d.001.sh
fi
"""

def create_timestring(estimated_runtime):
	estimated_hours, remaining_minutes = divmod(estimated_runtime, 3600)
	estimated_minutes, estimated_seconds = divmod(remaining_minutes, 60)
	timestring = "%02d:%02d:%02d" % (estimated_hours, estimated_minutes, estimated_seconds)

	return timestring

walltime_seconds = 14400
walltime_string = create_timestring(walltime_seconds)

random.seed()
numpy.random.seed()

hky_locus_names = list(hky_loci.keys())
numpy.random.shuffle(hky_locus_names)

n_reps = 30

templates_folder = "templates"
all_templates = {}
for folder_item in os.listdir(templates_folder):
	if folder_item.endswith(".xml"):
		template_name = folder_item.rsplit(".", 1)[0]
		template_path = os.path.join(templates_folder, folder_item)
		template_xml = xml.etree.ElementTree.parse(template_path)
		all_templates[template_name] = template_xml

startup_filename = "start.sh"

print(len(hky_locus_names))

for rep_i in range(n_reps):
	numpy.random.shuffle(hky_locus_names)

	rep_folder = "rep-%02d" % (rep_i)
	if not os.path.exists(rep_folder):
		os.mkdir(rep_folder)

	startup_path = os.path.join(rep_folder, startup_filename)
	startup_file = open(startup_path, "w")

	for template_name in sorted(all_templates):
		template_xml = all_templates[template_name]
		root = template_xml.getroot()

		for data_element in root.findall("data"):
			locus_name = data_element.attrib["id"]
			locus_i = int(locus_name[3:])
			uce_name = hky_locus_names[locus_i]
			for sequence_element in data_element.findall("sequence"):
				xml_species, xml_individual = sequence_element.attrib["taxon"].split("_")
				species_i = int(xml_species[1:])
				individual_i = int(xml_individual[3:])
				croc_individual = crocidura[species_i][1][individual_i]
				sequence_element.attrib["value"] = hky_loci[uce_name][croc_individual]

		if rep_i == 0:
			for taxonset in root.iter("taxonset"):
				if taxonset.attrib["id"] == "taxonsuperset":
					for taxon in taxonset.findall("taxon"):
						species_i = int(taxon.attrib["id"][1:])
						species_name = crocidura[species_i][0]
						taxon.attrib["id"] = species_name

		random_seed = random.randint(0, 2**31 - 1)
		run_name = "%s.%010d" % (template_name, random_seed)

		for logger_element in root.iter("logger"):
			if "fileName" in logger_element.attrib:
				old_filename = logger_element.attrib["fileName"]
				logfile_ext = old_filename.rsplit(".", 1)[1]
				logger_element.attrib["fileName"] = run_name + "." + logfile_ext

		output_xml_filename = run_name + ".xml"
		output_xml_path = os.path.join(rep_folder, output_xml_filename)
		template_xml.write(output_xml_path)

		incomplete_filename = run_name + ".incomplete"
		incomplete_path = os.path.join(rep_folder, incomplete_filename)
		incomplete_file = open(incomplete_path, "w")
		incomplete_file.close()

		pbs_script = pbs_template % {"walltime": walltime_string, "analysis": template_name, "seed": random_seed, "seconds": walltime_seconds, "incomplete_filename": incomplete_filename}
		script_filename = run_name + ".000.sh"
		script_path = os.path.join(rep_folder, script_filename)
		script_file = open(script_path, "w")
		script_file.write(pbs_script)
		script_file.close()

		startup_file.write("qsub %s\n" % (script_filename))

	startup_file.close()
