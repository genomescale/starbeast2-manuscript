#!/usr/bin/python2.7

import csv
import random
import numpy
import collections
import os
import string
import xml.etree.ElementTree

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

for rep_i in range(n_reps):
	rep_folder = "rep-%02d" % (rep_i)
	if not os.path.isdir(rep_folder):
		os.mkdir(rep_folder)

	startup_path = os.path.join(rep_folder, startup_filename)
	startup_file = open(startup_path, "w")

	for template_name in sorted(all_templates):
		template_xml = all_templates[template_name]
		root = template_xml.getroot()

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
