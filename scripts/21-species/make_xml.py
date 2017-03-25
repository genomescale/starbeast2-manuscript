#!/usr/bin/python2.7

import csv
import random
import numpy
import collections
import os
import string
import xml.etree.ElementTree
import sys

ambig_code_map = {
	"AA": "A",
	"CC": "C",
	"GG": "G",
	"TT": "T",
	"AC": "M",
	"AG": "R",
	"AT": "W",
	"CG": "S",
	"CT": "Y",
	"GT": "K",
}

ambig_codes = {}
for haploid_codes, diploid_code in ambig_code_map.items():
	nt_a, nt_b = haploid_codes
	if nt_a in ambig_codes:
		ambig_codes[nt_a][nt_b] = diploid_code
	else:
		ambig_codes[nt_a] = {nt_b: diploid_code}

	if nt_b in ambig_codes:
		ambig_codes[nt_b][nt_a] = diploid_code
	else:
		ambig_codes[nt_b] = {nt_a: diploid_code}

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

walltime_seconds = 18000
walltime_string = create_timestring(walltime_seconds)

random.seed()
numpy.random.seed()

n_reps = 30

templates_folder = "templates"
all_templates = {}
for folder_item in os.listdir(templates_folder):
	if folder_item.endswith(".xml"):
		template_name = folder_item.rsplit(".", 1)[0]
		template_path = os.path.join(templates_folder, folder_item)
		template_xml = xml.etree.ElementTree.parse(template_path)
		all_templates[template_name] = template_xml

sequences_filename = "simulated.sequences"
startup_filename = "start.sh"

for rep_i in range(n_reps):
	rep_folder = "rep-%02d" % (rep_i)

	startup_path = os.path.join(rep_folder, startup_filename)
	startup_file = open(startup_path, "w")

	sequences_path = os.path.join(rep_folder, sequences_filename)
	sequences_file = open(sequences_path)
	haploid_loci = []
	haploid_locus = None
	for l in sequences_file.readlines():
		if l.startswith(" ") or l.strip() == "":
			if haploid_locus != None:
				haploid_loci.append(haploid_locus)
			haploid_locus = {}
		else:
			taxon, sequence = l.strip().lstrip("'").split("'")
			species, tip = taxon.split("_")
			if species in haploid_locus:
				haploid_locus[species][tip] = sequence
			else:
				haploid_locus[species] = {tip: sequence}

	haploid_loci.append(haploid_locus)
	sequences_file.close()

	diploid_loci = []
	for haploid_locus in haploid_loci:
		diploid_locus = {}
		for species in sorted(haploid_locus):
			sequence_a = haploid_locus[species]["tip0"]
			sequence_b = haploid_locus[species]["tip1"]

			locus_length = len(sequence_a)
			diploid_sequence = ""
			for i in range(locus_length):
				nt_a = sequence_a[i]
				nt_b = sequence_b[i]
				diploid_sequence += ambig_codes[nt_a][nt_b]

			diploid_locus[species] = diploid_sequence

		diploid_loci.append(diploid_locus)

	for template_name in sorted(all_templates):
		template_xml = all_templates[template_name]
		root = template_xml.getroot()

		for data_element in root.findall("data"):
			locus_name = data_element.attrib["id"]
			locus_i = int(locus_name.lstrip("loc"))
			for sequence_element in data_element.findall("sequence"):
				taxon_name = sequence_element.attrib["taxon"]
				if template_name.startswith("beast"):
					if template_name.endswith("0"):
						sequence_element.attrib["value"] = haploid_loci[locus_i][taxon_name]["tip0"]
					else:
						sequence_element.attrib["value"] = diploid_loci[locus_i][taxon_name]
				else:
					species, tip = taxon_name.split("_")
					sequence_element.attrib["value"] = haploid_loci[locus_i][species][tip]

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
