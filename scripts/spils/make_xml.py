#!/usr/bin/python2.7

import random
import collections
import os
import concatenation
import starbeast2
import string
import sys
import xml.etree.ElementTree

n_replicates = int(sys.argv[1])

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
		qsub %(analysis)s-%(seed)010d.001.sh
fi
"""

def create_timestring(estimated_runtime):
	estimated_hours, remaining_minutes = divmod(estimated_runtime, 3600)
	estimated_minutes, estimated_seconds = divmod(remaining_minutes, 60)
	timestring = "%02d:%02d:%02d" % (estimated_hours, estimated_minutes, estimated_seconds)

	return timestring

def read_phylip(phylip_path):
	phylip_file = open(phylip_path)
	all_loci = []
	l = phylip_file.readline()
	locus = {}
	while l != "":
		if l.startswith(" "):
			if len(locus) > 0:
				all_loci.append(locus)
				locus = {}
		elif l.startswith("'"):
			label, sequence = l.lstrip("'").rstrip().split("'")
			locus[label] = sequence

		l = phylip_file.readline()

	if len(locus) > 0:
		all_loci.append(locus)

	return all_loci

walltime_seconds = 18000
walltime_string = create_timestring(walltime_seconds)

n_samples = 8192

simulated_sequences_filename = "simulated.sequences"
startup_filename = "start.sh"

limited_configurations = {
	"u": (100, 1, 0, 1, 0, 1),
	"v": (100, 0, 0, 0, 1, 1)
}

random.seed()

for replicate_i in range(n_replicates):
	replicate_folder = "rep-%02d" % (replicate_i)
	simulated_sequences_path = os.path.join(replicate_folder, simulated_sequences_filename)
	startup_path = os.path.join(replicate_folder, startup_filename)

	rep_loci = {}
	ambig_loci = {}
	unstructured_loci = read_phylip(simulated_sequences_path)
	for locus_i, locus in enumerate(unstructured_loci):
		locus_name = "loc%03d" % (locus_i)
		rep_loci[locus_name] = {}
		ambig_loci[locus_name] = {}
		for individual_name in sorted(locus):
			species_name, tip_name = individual_name.split("_")
			sequence = locus[individual_name]
			rep_loci[locus_name][species_name] = {individual_name: sequence}
			ambig_loci[locus_name][species_name] = sequence

	startup_file = open(startup_path, "w")
	startup_file.write("#!/bin/bash\n\n")

	for configuration_code in sorted(limited_configurations):
		analysis_name = "spils-" + configuration_code
		n_loci, analytical_pop_sizes, coordinated_topology, coordinated_heights, concatenate, species_tree_rates = limited_configurations[configuration_code]

		if concatenate == 1: sample_rate = 512
		else: sample_rate = 2048
		chain_length = n_samples * sample_rate

		random_seed = random.randint(0, 2**31 - 1)
		run_name = "%s-%010d" % (analysis_name, random_seed)
		incomplete_filename = run_name + ".incomplete"
		incomplete_path = os.path.join(replicate_folder, incomplete_filename)
		incomplete_file = open(incomplete_path, "w")
		incomplete_file.close()

		pbs_script = pbs_template % {"walltime": walltime_string, "analysis": analysis_name, "seed": random_seed, "seconds": walltime_seconds, "incomplete_filename": incomplete_filename}
		script_filename = run_name + ".000.sh"
		script_path = os.path.join(replicate_folder, script_filename)
		script_file = open(script_path, "w")
		script_file.write(pbs_script)
		script_file.close()
		startup_file.write("qsub " + script_filename + "\n")

		chain_xml_filename = run_name + ".xml"
		chain_xml_path = os.path.join(replicate_folder, chain_xml_filename)
		chain_xml_file = open(chain_xml_path, "w")
		if concatenate:
			chain_xml = concatenation.make_concat_xml(n_loci, ambig_loci, run_name, sample_rate, chain_length,
				relaxed_clock = True,
				continuous_rates = False,
				species_tree_rates = bool(species_tree_rates),
				adaptive_operators = False)
		else:
			chain_xml = starbeast2.make_starbeast_xml(n_loci, rep_loci, run_name, sample_rate, chain_length,
				relaxed_clock = True,
				continuous_rates = False,
				species_tree_rates = bool(species_tree_rates),
				analytical_pop_sizes = bool(analytical_pop_sizes),
				coordinated_topology = bool(coordinated_topology),
				coordinated_heights = bool(coordinated_heights),
				adaptive_operators = False)
		chain_xml.write(chain_xml_file, encoding = "UTF-8", xml_declaration = True)
		chain_xml_file.close()

	startup_file.close()
