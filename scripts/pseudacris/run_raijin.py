import csv
import os
import resource
import shutil
import subprocess
import sys

def append_log(input_log_path, output_log_path, sample_modulo):
	temp_log_path = output_log_path + ".new"
	temp_log_file = open(temp_log_path, "w")
	temp_log_writer = csv.writer(temp_log_file, dialect = csv.excel_tab)

	existing_log_file = open(output_log_path)
	existing_log_reader = csv.reader(existing_log_file, dialect = csv.excel_tab)
	next_sample_number = -1
	for existing_row in existing_log_reader:
		temp_log_writer.writerow(existing_row)
		next_sample_number += 1
	existing_log_file.close()

	input_log_file = open(input_log_path)
	input_log_line = input_log_file.readline()
	while input_log_line.startswith("#"):
		input_log_line = input_log_file.readline()

	input_log_reader = csv.reader(input_log_file, dialect = csv.excel_tab)
	for input_row_i, input_row in enumerate(input_log_reader):
		if (input_row_i > 0) and (input_row_i % sample_modulo == 0):
			input_row[0] = str(next_sample_number)
			temp_log_writer.writerow(input_row)
			next_sample_number += 1

	input_log_file.close()
	temp_log_file.close()

	os.rename(temp_log_path, output_log_path)

def copy_log(input_log_path, output_log_path, sample_modulo):
	temp_log_path = output_log_path + ".new"
	temp_log_file = open(temp_log_path, "w")
	temp_log_writer = csv.writer(temp_log_file, dialect = csv.excel_tab)

	input_log_file = open(input_log_path)
	input_log_line = input_log_file.readline()
	while input_log_line.startswith("#"):
		input_log_line = input_log_file.readline()

	header_row = input_log_line.strip().split()
	temp_log_writer.writerow(header_row)

	input_log_reader = csv.reader(input_log_file, dialect = csv.excel_tab)
	next_sample_number = 0
	for input_row_i, input_row in enumerate(input_log_reader):
		if input_row_i % sample_modulo == 0:
			input_row[0] = str(next_sample_number)
			temp_log_writer.writerow(input_row)
			next_sample_number += 1

	input_log_file.close()
	temp_log_file.close()

	os.rename(temp_log_path, output_log_path)

def append_trees(input_trees_path, output_trees_path, sample_modulo):
	temp_trees_path = output_trees_path + ".new"
	temp_trees_file = open(temp_trees_path, "w")

	existing_trees_file = open(output_trees_path)
	next_sample_number = 0
	nexus_header = True
	for existing_line in existing_trees_file.readlines():
		if existing_line.startswith("tree"):
			temp_trees_file.write(existing_line)
			nexus_header = False
			next_sample_number += 1
		elif nexus_header:
			temp_trees_file.write(existing_line)

	existing_trees_file.close()

	input_trees_file = open(input_trees_path)
	input_sample_number = 0
	for input_line in input_trees_file.readlines():
		if input_line.strip().startswith("tree"):
			if (input_sample_number > 0) and ((input_sample_number % sample_modulo) == 0):
				input_newick = input_line.strip().split()[-1]
				new_input_line = "tree STATE_%d = %s\n" % (next_sample_number, input_newick)
				temp_trees_file.write(new_input_line)
				next_sample_number += 1
			input_sample_number += 1

	temp_trees_file.write("End;\n")

	input_trees_file.close()
	temp_trees_file.close()

	os.rename(temp_trees_path, output_trees_path)

def copy_trees(input_trees_path, output_trees_path, sample_modulo):
	temp_trees_path = output_trees_path + ".new"
	temp_trees_file = open(temp_trees_path, "w")

	input_trees_file = open(input_trees_path)
	next_sample_number = 0
	input_sample_number = 0
	nexus_header = True
	for input_line in input_trees_file.readlines():
		if input_line.strip().startswith("tree"):
			if (input_sample_number % sample_modulo) == 0:
				input_newick = input_line.strip().split()[-1]
				new_input_line = "tree STATE_%d = %s\n" % (next_sample_number, input_newick)
				temp_trees_file.write(new_input_line)
				next_sample_number += 1
			nexus_header = False
			input_sample_number += 1
		elif nexus_header:
			temp_trees_file.write(input_line)

	temp_trees_file.write("End;\n")

	input_trees_file.close()
	temp_trees_file.close()

	os.rename(temp_trees_path, output_trees_path)

def create_timestring(estimated_runtime):
	estimated_hours, remaining_minutes = divmod(estimated_runtime, 3600)
	estimated_minutes, estimated_seconds = divmod(remaining_minutes, 60)
	timestring = "%02d:%02d:%02d" % (estimated_hours, estimated_minutes, estimated_seconds)

	return timestring

def read_ess(ess_path):
	ess_file = open(ess_path)
	ess_values = {}
	for l in ess_file.readlines():
		l_split = l.strip().split()
		if len(l_split) >= 9:
			statistic_name = l_split[0]
			statistic_ess = l_split[6]
			if statistic_ess[0].isdigit():
				ess_values[statistic_name] = float(statistic_ess)

	ess_file.close()

	return ess_values

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
		%(next_python_cmd)s
fi

if [ -f %(incomplete_filename)s ]
	then
		%(next_qsub_cmd)s
fi
"""
cmd_template = "python ../run_raijin.py %s %d %d %d"

java_bin = "/path/to/jre1.8.0_101/bin/java"
beast_jar = "/path/to/beast/lib/beast.jar"
loganalyser_bin = "/path/to/BEASTv1.8.4/bin/loganalyser"
n_burnin_states = 1024

chain_modulis = {0: 0} # I ain't stepping out of shit, all my Latin's legit
last_chains = [0]
for exp in range(12):
	rate_mod_exp = exp + 1
	first_chain = 2**exp
	next_first = 2**rate_mod_exp
	chain_range = range(first_chain, next_first)
	last_chains.append(next_first - 1)
	for chain_i in chain_range:
		chain_modulis[chain_i] = rate_mod_exp

analysis_name   = sys.argv[1]
random_seed     = int(sys.argv[2])
chain_number    = int(sys.argv[3])
runtime_seconds = int(sys.argv[4])

chain_modulo_exp = chain_modulis[chain_number]
chain_modulo = 2**chain_modulo_exp

run_name = "%s.%010d" % (analysis_name, random_seed)
beast_xml_filename  = run_name + ".xml"
incomplete_filename = run_name + ".incomplete"
log_filename        = run_name + ".log"
trees_filename      = run_name + ".trees"
state_filename      = run_name + ".state"
resource_filename   = run_name + ".csv"

chain_name = "%s.%03x" % (run_name, chain_number)
time_filename         = chain_name + ".time"
backup_state_filename = chain_name + ".state"
script_filename       = chain_name + ".sh"

modulo_name = "%s.%02d" % (run_name, chain_modulo_exp)
modulo_ess_filename   = modulo_name + ".ess"
modulo_log_filename   = modulo_name + ".log"
modulo_trees_filename = modulo_name + ".trees"

if chain_number == 0:
	beast_args = ["-overwrite", "-seed", str(random_seed), "-statefile", state_filename, "-threads", "1", beast_xml_filename]
else:
	beast_args = ["-resume", "-seed", str(random_seed), "-statefile", state_filename, "-threads", "1", beast_xml_filename]

beast_cmd = [java_bin, "-jar", beast_jar] + beast_args
subprocess.check_call(beast_cmd)
shutil.copy(state_filename, backup_state_filename)

resource_file = open(resource_filename, "a")
resource_writer = csv.writer(resource_file)
if chain_number == 0:
	resource_header = ["ru_utime", "ru_stime", "ru_maxrss", "ru_ixrss", "ru_idrss", "ru_isrss", "ru_minflt", "ru_majflt", "ru_nswap", "ru_inblock", "ru_oublock", "ru_msgsnd", "ru_msgrcv", "ru_nsignals", "ru_nvcsw", "ru_nivcsw"]
	resource_writer.writerow(resource_header)
resource_usage = tuple(resource.getrusage(resource.RUSAGE_CHILDREN))
resource_writer.writerow(resource_usage)
resource_file.close()

if chain_number == 0:
	copy_log(log_filename, modulo_log_filename, 1)
	copy_trees(trees_filename, modulo_trees_filename, 1)
else:
	append_log(log_filename, modulo_log_filename, chain_modulo)
	append_trees(trees_filename, modulo_trees_filename, chain_modulo)
	previous_script_filename = "%s.%03x.sh" % (run_name, chain_number - 1)
	os.remove(previous_script_filename)

insufficient_ess = True
if chain_number in last_chains:
	loganalyser_cmd = [loganalyser_bin, "-burnin", str(n_burnin_states), modulo_log_filename, modulo_ess_filename]
	subprocess.check_call(loganalyser_cmd)

	ess_values = read_ess(modulo_ess_filename).values()
	if min(ess_values) >= 200.0:
		os.remove(incomplete_filename)
		os.remove(script_filename)
		os.remove(beast_xml_filename)
		insufficient_ess = False
	else:
		next_log_filename = "%s.%02d.log" % (run_name, chain_modulo_exp + 1)
		next_trees_filename = "%s.%02d.trees" % (run_name, chain_modulo_exp + 1)
		copy_log(modulo_log_filename, next_log_filename, 2)
		copy_trees(modulo_trees_filename, next_trees_filename, 2)

	compressed_filename = modulo_name + ".tar.xz"
	tar_cmd = ["tar", "--remove-files", "-cJf", compressed_filename, modulo_log_filename, modulo_trees_filename]
	subprocess.check_call(tar_cmd)

os.remove(log_filename)
os.remove(trees_filename)

if insufficient_ess:
	next_script_filename = "%s.%03x.sh" % (run_name, chain_number + 1)
	subsequent_script_filename = "%s.%03x.sh" % (run_name, chain_number + 2)

	python_cmd = cmd_template % (analysis_name, random_seed, chain_number + 1, runtime_seconds)
	qsub_cmd = "qsub " + subsequent_script_filename
	walltime_string = create_timestring(runtime_seconds)
	next_script_content = pbs_template % {"walltime": walltime_string, "next_python_cmd": python_cmd, "next_qsub_cmd": qsub_cmd, "incomplete_filename": incomplete_filename}

	next_script_file = open(next_script_filename, "w")
	next_script_file.write(next_script_content)
	next_script_file.close()
