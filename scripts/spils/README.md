# Description of SPILS data scripts

These are all scripts used to analyse the effect of SPILS on concatenation
and StarBEAST2 using simulated data.

## branch_rates.R

Generates Figure 2. `spils_branch_rates.py` must be run first.

## concatenation.py

Not a separate script, but imported by `make_xml.py` to generate concatenation
XML files.

## make_xml.py

Makes concatenation and StarBEAST2 XML files with simulated sequences.
`simulate.sh` must be run first.

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

Runs each simulate python script with the settings used for our SPILS study.

## spils_branch_rates.py

Calculates estimated branch rates and collates with true rates. Also calculates
false positive branch rates used in Figures S7 and S8. Path to the BEAST
treeannotator binary must be supplied as a command line argument.

## starbeast2.py

Not a separate script, but imported by `make_xml.py` to generate StarBEAST2
XML files.