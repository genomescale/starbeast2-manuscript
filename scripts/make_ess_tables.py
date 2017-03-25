#!/usr/bin/python2.7

import collections
import csv
import os
import sys

analysis_folder = sys.argv[1]
analysis_filenames = os.listdir(analysis_folder)
analysis_filenames.sort()

all_stats = set()
all_keys = set()

per_hour_means = collections.defaultdict(dict)
per_hour_stdevs = collections.defaultdict(dict)

per_mstates_means = collections.defaultdict(dict)
per_mstates_stdevs = collections.defaultdict(dict)

for fn in analysis_filenames:
	file_base, file_ext = os.path.splitext(fn)
	fn_split = file_base.split("-")
	if (len(fn_split) == 3):
		stat_name, ess_name, summary_stat = fn_split
		if ess_name.startswith("ess_per_") and (file_ext == ".csv"):
			print("*" + fn)
			all_stats.add(stat_name)

			summary_path = os.path.join(analysis_folder, fn)
			summary_file = open(summary_path)
			summary_file.readline()

			summary_reader = csv.reader(summary_file)
			for row in summary_reader:
				k = row[0]
				all_keys.add(k)

				v = float(row[1])

				if ess_name == "ess_per_hour":
					if summary_stat == "mean":
						per_hour_means[k][stat_name] = v
					else:
						per_hour_stdevs[k][stat_name] = v
				else:
					if summary_stat == "mean":
						per_mstates_means[k][stat_name] = v
					else:
						per_mstates_stdevs[k][stat_name] = v

			summary_file.close()

def key_sorter(x, y):
	x_method, x_loci, x_config = x.split("-")
	y_method, y_loci, y_config = y.split("-")

	if x_method > y_method:
		return 1
	elif x_method < y_method:
		return -1
	elif x_loci > y_loci:
		return 1
	elif x_loci < y_loci:
		return -1
	elif int(x_config) > int(y_config):
		return 1
	elif int(x_config) < int(y_config):
		return -1
	else:
		return 0

ess_per_hour_filename = analysis_folder.rstrip("/") + ".ess_per_hour.tex"
ess_per_hour_file = open(ess_per_hour_filename, "w")

begin_array_line = "$\\begin{array}{|c|"
for i in range(len(all_stats)):
	begin_array_line += "c|"

begin_array_line += "}\n"
ess_per_hour_file.write(begin_array_line)

header_line = "\\hline\n\\text{key}"
for stat_name in sorted(all_stats):
	header_line += " & \\text{%s}" % (stat_name)

header_line += "\\tabularnewline\n\\hline\n"
ess_per_hour_file.write(header_line)

for k in sorted(all_keys, cmp = key_sorter):
	k_line = "\\text{%s}" % (k)
	for stat_name in sorted(all_stats):
		if stat_name in per_hour_means[k]:
			ess_per_hour_mean = per_hour_means[k][stat_name]
			ess_per_hour_stdev = per_hour_stdevs[k][stat_name]
			k_line += " & %.3f\\pm%.3f" % (ess_per_hour_mean, ess_per_hour_stdev)
		else:
			k_line += " & NA"

	k_line += "\\tabularnewline\n\\hline\n"
	ess_per_hour_file.write(k_line)

ess_per_hour_file.write("\\end{array}$\n")

ess_per_hour_file.close()

# another verse, same as the first

ess_per_mstates_filename = analysis_folder.rstrip("/") + ".ess_per_mstates.tex"
ess_per_mstates_file = open(ess_per_mstates_filename, "w")

begin_array_line = "$\\begin{array}{|c|"
for i in range(len(all_stats)):
	begin_array_line += "c|"

begin_array_line += "}\n"
ess_per_mstates_file.write(begin_array_line)

header_line = "\\hline\n\\text{key}"
for stat_name in sorted(all_stats):
	header_line += " & \\text{%s}" % (stat_name)

header_line += "\\tabularnewline\n\\hline\n"
ess_per_mstates_file.write(header_line)

for k in sorted(all_keys, cmp = key_sorter):
	k_line = "\\text{%s}" % (k)
	for stat_name in sorted(all_stats):
		if stat_name in per_mstates_means[k]:
			ess_per_mstates_mean = per_mstates_means[k][stat_name]
			ess_per_mstates_stdev = per_mstates_stdevs[k][stat_name]
			k_line += " & %.3f\\pm%.3f" % (ess_per_mstates_mean, ess_per_mstates_stdev)
		else:
			k_line += " & NA"

	k_line += "\\tabularnewline\n\\hline\n"
	ess_per_mstates_file.write(k_line)

ess_per_mstates_file.write("\\end{array}$\n")

ess_per_mstates_file.close()
