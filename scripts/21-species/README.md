# Description of simulated data scripts

These are all scripts used to analyse the computational performance and
statistical accuracy of StarBEAST2 using simulated data.

## branch_length_accuracy_phased.R

Generates Figure 5. `calculate_branch_length_accuracy.py` in the root folder must be run first.

## branch_length_accuracy_unphased.R

Generates Figure S10. `calculate_branch_length_accuracy.py` in the root folder must be run first.

## branch_rates_phased.R

Generates Figure 7. Also the R-squared values for estimated branch rates using
phased data for concatenation. `calculate_branch_rates.py` in the root folder
must be run first.

## branch_rates_unphased.R

Generates Figure S12. Also the R-squared values for estimated branch rates
using unphased data with ambiguity codes for concatenation.
`calculate_branch_rates.py` in the root folder must be run first.

## make_xml.py

Makes BEAST XML files with simulated sequences, based on every XML file in the
`templates` folder. It also generates a `start.sh` script for every replicate
used to begin all chains in that folder. `simulate.sh` must be run first.

## run_raijin.py

This script was used to run chains on the NCI Raijin cluster, which uses the
PBS scheduling software. It runs a chain for a given length and checks ESS
values. If any value is below 200, it doubles the length of the chain and
checkes again.

## simulate_gene_trees.py

Simulates gene trees evolving within species trees. Branch lengths of the
gene trees are shrunk and stretched according to the species tree branch rates.

## simulate_sequences.py

Simulates multiple sequence alignments based on simulated gene trees.

## simulate_species_trees.py

Simulates species trees according to a birth-death model.

## simulate.sh

Runs each simulate python script with the settings used for our simulation
study.

## topology_accuracy_phased.R

Generates Figure 6. `calculate_topology_accuracy.py` in the root folder must be run first.

## topology_accuracy_unphased.R

Generates Figure S11. `calculate_topology_accuracy.py` in the root folder must be run first.