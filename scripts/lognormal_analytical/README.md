# Scripts for testing correctness of analytical population sizes

These are all scripts used to test the correctness of the StarBEAST2
implementation of the multispecies coalescent and species tree relaxed clocks.

## complete_gene_trees.py

Adds simulated clock rates based on containing species tree branch rates,
multiplied by some factor.

## complete_species_trees.py

Adds simulated lognormal relaxed clock branch rates to each simulated species
tree.

## plot_gene_branch_rates.R

Generates Figure S5. `read_trees.py` must be run first on both gene trees files.

## plot_gene_node_heights.R

Generates Figure S4. `read_trees.py` must be run first on both gene trees files.

## plot_gene_topology_frequencies.R

Generates Figure S2. `read_trees.py` must be run first on both gene trees files.

## plot_species_node_heights.R

Generates Figure S3. `read_trees.py` must be run first on the species trees.

## plot_species_topology_frequencies.R

Generates Figure S1. `read_trees.py` must be run first on the species trees.

## read_trees.py

Extracts and collates sampled topologies, node heights and branch rates for a
posterior distribution of species trees or gene trees. Must be run separately
for each trees file.

## simulate.sh

Simulates 100,000 species trees and 200,000 gene trees using biopy, then runs
the various "complete" scripts with the parameters used in the paper.

## test.xml

StarBEAST2 XML file to sample from the prior, in order to test correctness of
its implementation.