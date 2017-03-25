#!/bin/bash

echo "Simulating species trees..."
python simulate_species_trees.py 96

echo "Simulating gene trees..."
python simulate_gene_trees.py 96 100

echo "Simulating HKY sequences..."
python simulate_sequences.py 96 1000

echo "Generating BEAST XML"
python make_xml.py 96
