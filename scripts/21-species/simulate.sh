#!/bin/bash

echo "Simulating species trees..."
python ./simulate_species_trees.py 21 30

echo "Simulating gene trees..."
python ./simulate_gene_trees.py 30 260

echo "Simulating HKY sequences..."
python ./simulate_sequences.py 30 600
