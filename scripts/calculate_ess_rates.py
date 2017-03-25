#!/usr/bin/python2.7

import csv
import os
import string
import sys
import subprocess

input_folder = sys.argv[1]

if input_folder.startswith("21-species"):
	dataset_name = "Simulation"
else:
	dataset_name = input_folder.rstrip("/").capitalize()

starbeast2_template = """library(scales)
library(ggplot2)
library(cowplot)

ess_rates = read.csv("%(safe_name)s-ess_rates.csv")

ess_rates$popsize_integration = factor(ess_rates$popsize_integration, levels = c("mcmc", "analytical"))
ess_rates$clock = factor(ess_rates$clock, levels = c("strict", "gt-ucln", "st-ucln"))
ess_rates$operators = factor(ess_rates$operators, levels = c("neither", "topology_changing", "height_changing", "both"))

ess_rates$coordinated_height_operators = ess_rates$operators == "both" | ess_rates$operators == "height_changing"
ess_rates$coordinated_topology_operators = ess_rates$operators == "both" | ess_rates$operators == "topology_changing"

starbeast2_rates = subset(ess_rates, method == "starbeast2" & loci == "1x")

sink("%(safe_name)s-%(ess_rate)s-starbeast2.txt")
summary(lm(log(%(ess_rate)s) ~ coordinated_topology_operators + coordinated_height_operators + popsize_integration, data = subset(starbeast2_rates, clock == "strict")))
summary(lm(log(%(ess_rate)s) ~ coordinated_topology_operators + coordinated_height_operators + popsize_integration, data = subset(starbeast2_rates, clock == "gt-ucln")))
summary(lm(log(%(ess_rate)s) ~ coordinated_topology_operators + coordinated_height_operators + popsize_integration, data = subset(starbeast2_rates, clock == "st-ucln")))
sink()

levels(starbeast2_rates$popsize_integration) = c("MCMC population size integration", "Analytical population size integration")
levels(starbeast2_rates$clock) = c("Strict clock", "Gene tree relaxed clocks", "Species tree relaxed clock")
levels(starbeast2_rates$operators) = c("Neither", "Topology changing", "Height changing", "Both")

ggplot(starbeast2_rates, aes(y = %(ess_rate)s, x = operators, fill = popsize_integration)) +
	geom_boxplot() +
	xlab("Coordinated MCMC operators") +
	ylab("%(ess_title)s") +
	scale_y_log10(breaks = 2^seq(-20, 20, 2)) +
	facet_grid(. ~ clock) +
	theme(axis.text.x = element_text(angle = 45, hjust = 1.0),
		legend.position = "top",
		legend.title = element_blank()) +
	background_grid(major = "y", minor = "none")

ggsave("%(safe_name)s-%(ess_rate)s-starbeast2.pdf", units = "mm", height = 90, width = 180)
ggsave("%(safe_name)s-%(ess_rate)s-starbeast2.png", units = "mm", height = 90, width = 180)
"""

comparison_template = """library(scales)
library(ggplot2)
library(cowplot)

ess_rates = read.csv("%(safe_name)s-ess_rates.csv")
ess_rates = subset(ess_rates, sequence == "haploid")

starbeast_comparison = subset(ess_rates, method == "starbeast1" | method == "beast" | (operators == "height_changing" & popsize_integration == "analytical"))

starbeast_comparison$clock = factor(starbeast_comparison$clock, levels = c("strict", "gt-ucln", "st-ucln"))
levels(starbeast_comparison$clock) = c("Strict clock", "GT-UCLN", "ST-UCLN")

starbeast_comparison$method = factor(starbeast_comparison$method, levels = c("starbeast1", "starbeast2", "beast"))
levels(starbeast_comparison$method) = c("*BEAST", "StarBEAST2", "Concat'ion")

starbeast_comparison$loci = factor(starbeast_comparison$loci, levels = c("1x", "2x", "4x", "5x"))
levels(starbeast_comparison$loci) = c("%(onex_loci)d loci", "%(twox_loci)d loci", "%(fourx_loci)d loci", "%(fivex_loci)d loci")

starbeast_comparison$category = paste0(starbeast_comparison$method, "\n", starbeast_comparison$loci)
starbeast_comparison$category = factor(starbeast_comparison$category, levels = paste0(rep(levels(starbeast_comparison$method), c(4, 4, 4)), "\n", levels(starbeast_comparison$loci)))

ggplot(starbeast_comparison, aes(y = %(ess_rate)s, x = category, fill = clock)) +
	geom_boxplot() +
	xlab("Method and data") +
	ylab("%(ess_title)s") +
	scale_y_log10(limits = c(0.1, 1000.0), breaks = 10^seq(-10, 10, 1), labels = comma) +
	facet_grid(. ~ dataset) +
	theme(axis.text.x = element_text(angle = 45, hjust = 1.0),
		legend.position = "top",
		legend.title = element_blank()) +
	background_grid(major = "y", minor = "none")

ggsave("%(safe_name)s-%(ess_rate)s-comparison.pdf", units = "mm", height = 90, width = 90)
ggsave("%(safe_name)s-%(ess_rate)s-comparison.png", units = "mm", height = 90, width = 90)
"""

tables_template = """library(ggplot2)
ess_rates = read.csv("%(safe_name)s-ess_rates.csv")

mean_table = aggregate(log(%(ess_rate)s) ~ config, ess_rates, mean)
write.csv(mean_table, "%(safe_name)s-%(ess_rate)s-mean.csv", row.names = FALSE)

sd_table = aggregate(log(%(ess_rate)s) ~ config, ess_rates, sd)
write.csv(sd_table, "%(safe_name)s-%(ess_rate)s-sd.csv", row.names = FALSE)
"""

starbeast2_configurations = {}
for clock_i, clock in enumerate(["strict", "gt-ucln", "st-ucln"]):
	for coordinated_topology in [4, 0]:
		for coordinated_height in [2, 0]:
			for analytical_popsizes in [1, 0]:
				config = (10 * clock_i) + coordinated_topology + coordinated_height + analytical_popsizes
				config_str = "%02d" % config
				if coordinated_topology == 0:
					if coordinated_height == 0:
						operators = "neither"
					else:
						operators = "height_changing"
				else:
					if coordinated_height == 0:
						operators = "topology_changing"
					else:
						operators = "both"

				if analytical_popsizes == 0:
					popsize_integration = "mcmc"
				else:
					popsize_integration = "analytical"

				starbeast2_configurations[config_str] = [clock, operators, popsize_integration, "haploid"]

beast_configurations = {}
beast_configurations["00"] = ["strict", "NA", "NA", "haploid"]
beast_configurations["01"] = ["strict", "NA", "NA", "diploid"]
beast_configurations["20"] = ["st-ucln", "NA", "NA", "haploid"]
beast_configurations["21"] = ["st-ucln", "NA", "NA", "diploid"]

n_reps = 30
all_stats = {"minimum": []}
for rep_i in range(n_reps):
	rep_subfolder = "rep-%02d" % (rep_i)
	rep_folder = os.path.join(input_folder, rep_subfolder)

	rep_filenames = os.listdir(rep_folder)
	rep_filenames.sort()

	rep_group_numbers = {}
	rep_random_seeds = {}

	for fn in rep_filenames:
		if fn.endswith(".ess"):
			fn_split = fn.split(".")
			analysis_name = fn_split[0]
			random_seed = int(fn_split[1])
			group_i = int(fn_split[2])

			if (analysis_name not in rep_group_numbers) or (group_i > rep_group_numbers[analysis_name]):
				rep_random_seeds[analysis_name] = random_seed
				rep_group_numbers[analysis_name] = group_i

	for analysis_name in sorted(rep_group_numbers):
		method, loci, config = analysis_name.split("-")

		random_seed = rep_random_seeds[analysis_name]
		group_i = rep_group_numbers[analysis_name]
		n_chains = 2**group_i

		times_filename = "%s.%010d.csv" % (analysis_name, random_seed)
		times_path = os.path.join(rep_folder, times_filename)
		times_file = open(times_path)
		times_reader = csv.reader(times_file)

		total_cpu_time = 0.0
		row_i = 0
		for row in times_reader:
			if row[0] != "ru_utime":
				row_i += 1
				if row_i <= n_chains:
					total_cpu_time += float(row[0])
					total_cpu_time += float(row[1])

		times_file.close()

		post_burnin_hours = (total_cpu_time / 3600.0) * 0.875

		if method == "beast":
			if loci == "1x":
				n_states_exp = 22
			else:
				n_states_exp = 20

			if config[1] == "0":
				n_states_exp += 1
		else:
			if loci == "1x":
				n_states_exp = 24
			else:
				n_states_exp = 23

		n_states = 2**(n_states_exp + group_i)

		post_burnin_mstates = (float(n_states) / 1000000.0) * 0.875

		ess_filename = "%s.%010d.%02d.ess" % (analysis_name, random_seed, group_i)
		ess_path = os.path.join(rep_folder, ess_filename)
		ess_file = open(ess_path)

		minimum_ess_per_hour = 100000.0
		minimum_ess_per_mstates = 100000.0
		for line in ess_file.readlines():
			l_split = line.strip().split()
			if len(l_split) >= 9 and l_split[6][0].isdigit():
				statistic_name = l_split[0]
				if "." in statistic_name and ":" in statistic_name:
					statistic_name = statistic_name[:min(statistic_name.find("."), statistic_name.find(":"))]
				elif "." in statistic_name:
					statistic_name = statistic_name[:statistic_name.find(".")]
				elif ":" in statistic_name:
					statistic_name = statistic_name[:statistic_name.find(":")]
				statistic_name = statistic_name.lower()

				if statistic_name.startswith("popmean") or statistic_name.startswith("constpopmean"):
					statistic_name = "popmean"
				elif statistic_name.startswith("birthdeath") or statistic_name == "yulemodel":
					statistic_name = "speciestree"
				elif statistic_name == "birthrate2" or statistic_name == "speciationrate":
					statistic_name = "netdiversificationrate"
				elif statistic_name == "relativedeathrate2":
					statistic_name = "extinctionfraction"

				statistic_ess = float(l_split[6])
				ess_per_hour = statistic_ess / post_burnin_hours
				ess_per_mstates = statistic_ess / post_burnin_mstates
				minimum_ess_per_hour = min(minimum_ess_per_hour, ess_per_hour)
				minimum_ess_per_mstates = min(minimum_ess_per_mstates, ess_per_mstates)

				if method == "beast":
					output_row = [dataset_name, rep_i, method, analysis_name] + beast_configurations[config] + [loci, ess_per_hour, ess_per_mstates]
				else:
					output_row = [dataset_name, rep_i, method, analysis_name] + starbeast2_configurations[config] + [loci, ess_per_hour, ess_per_mstates]

				if statistic_name in all_stats:
					all_stats[statistic_name].append(output_row)
				else:
					all_stats[statistic_name] = [output_row]

		ess_file.close()

		if method == "beast":
			output_row = [dataset_name, rep_i, method, analysis_name] + beast_configurations[config] + [loci, minimum_ess_per_hour, minimum_ess_per_mstates]
		else:
			output_row = [dataset_name, rep_i, method, analysis_name] + starbeast2_configurations[config] + [loci, minimum_ess_per_hour, minimum_ess_per_mstates]

		all_stats["minimum"].append(output_row)

def make_safe(original_name):
	new_name = ""
	for c in original_name:
		if c in string.ascii_uppercase + string.ascii_lowercase:
			new_name += c
		else:
			new_name += "."

	return new_name

for stat_name in sorted(all_stats):
	safe_name = make_safe(stat_name)

	output_filename = safe_name + "-ess_rates.csv"
	output_path = os.path.join(input_folder, output_filename)
	output_file = open(output_path, "w")
	output_writer = csv.writer(output_file)

	header_row = ["dataset", "rep_i", "method", "config", "clock", "operators", "popsize_integration", "sequence", "loci", "ess_per_hour", "ess_per_mstates"]
	output_writer.writerow(header_row)

	for row in all_stats[stat_name]:
		output_writer.writerow(row)

	output_file.close()

	for ess_rate, ess_title in [["ess_per_hour", "ESS per hour"], ["ess_per_mstates", "ESS per million states"]]:
		if "crocidura" in input_folder:
			onex_loci = 50
		else:
			onex_loci = 26

		twox_loci  = onex_loci * 2
		fourx_loci = onex_loci * 4
		fivex_loci  = onex_loci * 5

		if "crocidura" not in input_folder:
			starbeast2_script = starbeast2_template % {"safe_name": safe_name, "ess_rate": ess_rate, "ess_title": ess_title}
			starbeast2_filename = "%s-%s-starbeast2.R" % (safe_name, ess_rate)
			starbeast2_path = os.path.join(input_folder, starbeast2_filename)
			starbeast2_file = open(starbeast2_path, "w")
			starbeast2_file.write(starbeast2_script)
			starbeast2_file.close()

			subprocess.check_call(["Rscript", starbeast2_filename], cwd = input_folder)

		comparison_script = comparison_template % {"safe_name": safe_name, "ess_rate": ess_rate, "ess_title": ess_title, "onex_loci": onex_loci, "twox_loci": twox_loci, "fourx_loci": fourx_loci, "fivex_loci": fivex_loci}
		comparison_filename = "%s-%s-comparison.R" % (safe_name, ess_rate)
		comparison_path = os.path.join(input_folder, comparison_filename)
		comparison_file = open(comparison_path, "w")
		comparison_file.write(comparison_script)
		comparison_file.close()

		subprocess.check_call(["Rscript", comparison_filename], cwd = input_folder)

		tables_script = tables_template % {"safe_name": safe_name, "ess_rate": ess_rate, "ess_title": ess_title, "onex_loci": onex_loci, "twox_loci": twox_loci, "fourx_loci": fourx_loci, "fivex_loci": fivex_loci}
		tables_filename = "%s-%s-tables.R" % (safe_name, ess_rate)
		tables_path = os.path.join(input_folder, tables_filename)
		tables_file = open(tables_path, "w")
		tables_file.write(tables_script)
		tables_file.close()

		subprocess.check_call(["Rscript", tables_filename], cwd = input_folder)
