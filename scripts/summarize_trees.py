#!/usr/bin/python2.7

import csv
import dendropy
import dendropy.calculate
import os
import subprocess
import sys

treeannotator_bin = sys.argv[1]
simulation_folder = "21-species"
reference_filename = "biopy.species.tree"

n_burnin = 1026
n_reps = 30

for rep_i in range(n_reps):
	rep_subfolder = "rep-%02d" % (rep_i)
	rep_folder = os.path.join(simulation_folder, rep_subfolder)

	rep_filenames = os.listdir(rep_folder)
	rep_filenames.sort(reverse = True)

	archive_filenames = [fn for fn in rep_filenames if fn.endswith(".tar.xz")]
	incomplete_run_names = set()

	for archive_filename in archive_filenames:
		run_name, random_seed, group_number_str, extension = archive_filename.split(".", 3)
		alt_name, n_loci, analysis_code = run_name.split("-")

		if run_name not in incomplete_run_names and (alt_name != "starbeast2" or analysis_code.endswith("3")):
			incomplete_run_names.add(run_name)
			group_number = int(group_number_str)

			truth_heights_filename = "%s.%s.truth_heights.tree" % (run_name, random_seed)
			mcc_heights_filename = "%s.%s.mcc_heights.tree" % (run_name, random_seed)

			truth_heights_path = os.path.join(rep_folder, truth_heights_filename)
			mcc_heights_path = os.path.join(rep_folder, mcc_heights_filename)

			if not (os.path.isfile(truth_heights_path) and os.path.isfile(mcc_heights_path)):
				posterior_filename = "%s.%s.%02d.trees" % (run_name, random_seed, group_number)
				extract_cmd = ["tar", "xJvf", archive_filename, posterior_filename]
				subprocess.check_call(extract_cmd, cwd = rep_folder)
				posterior_path = os.path.join(rep_folder, posterior_filename)
				posterior_file = open(posterior_path)

				post_burnin_filename = "%s.%s.post_burnin.trees" % (run_name, random_seed)
				post_burnin_path = os.path.join(rep_folder, post_burnin_filename)
				post_burnin_file = open(post_burnin_path, "w")

				tree_i = 0
				for l in posterior_file.readlines():
					l_split = l.strip().split(None, 3)
					if len(l_split) == 4 and l_split[0] == "tree":
						newick_string = l_split[3]
						tree_i += 1
						if tree_i >= n_burnin:
							tree_line = "tree STATE_%d = %s\n" % (tree_i - n_burnin, newick_string)
							post_burnin_file.write(tree_line)
					else:
						post_burnin_file.write(l)

				post_burnin_file.close()

				mcc_heights_cmd = [treeannotator_bin, "-heights", "ca", "-burnin", "0", post_burnin_filename, mcc_heights_filename]
				truth_heights_cmd = [treeannotator_bin, "-heights", "ca", "-burnin", "0", "-target", reference_filename, post_burnin_filename, truth_heights_filename]

				subprocess.check_call(mcc_heights_cmd, cwd = rep_folder)
				subprocess.check_call(truth_heights_cmd, cwd = rep_folder)

				os.remove(posterior_path)
				os.remove(post_burnin_path)
