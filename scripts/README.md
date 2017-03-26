# Description of scripts in the root folder

These are all scripts used to analyse the computational performance and
statistical accuracy of StarBEAST2. There are most scripts in each subfolder
which are specific to the corresponding analysis.

## calculate_branch_length_accuracy.py

Extracts the species tree branch lengths inferred by StarBEAST2 and other
methods, summarised using the common ancestor method, from simulated data in
the `21-species` folder. Also extracts the boundaries of highest posterior
density (HPD) credible intervals for branch lengths. `summarize_trees.py` must
be run first.

## calculate_branch_rates.py

Extracts the species tree branch clock rates inferred by StarBEAST2 and other
methods, summarised using the common ancestor method, from simulated data in
the `21-species` folder. Also extracts the boundaries of highest posterior
density (HPD) credible intervals for branch rates. `summarize_trees.py` must
be run first.

## calculate_ess_rates.py

Calculates ESS rates based on ESS values and CPU time logged for each chain.
Also generates R scripts to make plots of ESS rates, for each statistic, to
compare the computational performance of StarBEAST2 with other methods.

## calculate_topology_accuracy.py

Collates topological accuracy and coverage from simulated data in the
`21-species` folder. Accuracy is based on the Robinson-Foulds distance between
the true simulated topology and the maximum clade credibility (MCC) tree from
the posterior distribution. Coverage is based on presence of the true
simulated topology in the 95% credible set of topologies. `measure_trees.sh`
must be run first.

## make_ess_tables.py

Generates LaTeX tables of ESS rate means and standard deviations. Used for
simulated data in the `21-species` folder as well as empirical data sets
in the `pseudacris` and `crocidura` folders. `calculate_ess_rates.py` must
be run first.

## measure_trees.py

Summarizes the posterior probability of each species tree topology present in
a posterior distribution inferred from simulated data. Also records the
distance of each present topology from the maximum clade credibility (MCC)
topology, and from the true simulated topology.  Uses a lot of memory so only
run on a single replicate at a time.

## measure_trees.sh

Runs `measure_trees.py` sequentially for every replicate.

## minimum-ess_per_hour-comparison.R

Generates Figure 3. `calculate_ess_rates.py` must be run first.

## minimum-ess_per_hour-starbeast2.R

Generates Figure 4. `calculate_ess_rates.py` must be run first.

## mstates_per_hour.R

Generates Figure S9. `calculate_ess_rates.py` must be run first.

## summarize_trees.py

Runs `treeannotator` twice for each posterior distribution of species trees
inferred from simulated data. One run is to calculate point estimates and HPDs
of branch lengths, conditional on the true species tree topology, using the
common ancestor method. The other run is to find the maximum clade credibility
(MCC) tree.
