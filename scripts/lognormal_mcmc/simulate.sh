#!/bin/bash

# rm biopy.species.trees
# rm biopy.gene0.trees
# rm biopy.gene1.trees

# echo "Simulating species trees..."
# sample_species_tree -o biopy.species.trees -p i,3,0.004 -b 200.0 -d 100.0 -n 100000 species0,species1,species2,species3,species4

# echo "Simulating half rate gene trees..."
# genetree_in_sptree_sim -t 1 -o biopy.gene0.trees biopy.species.trees

# echo "Simulating double rate gene trees..."
# genetree_in_sptree_sim -t 1 -o biopy.gene1.trees biopy.species.trees

echo "Adding rates to species trees..."
python complete_species_trees.py

echo "Adding rates to half rate gene trees..."
python complete_gene_trees.py 0 0.5

echo "Adding rates to double rate gene trees..."
python complete_gene_trees.py 1 2.0
