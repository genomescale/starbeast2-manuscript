# Description of _Crocidura_ files

These are all scripts used to analyse the computational performance
of StarBEAST2 using _Crocidura_ data.

## convert_to_fasta.py

Used to convert _Crocidura_ nexus files into fasta format. A zip file
containing _Crocidura_ sequences can be downloaded from
http://dx.doi.org/10.5061/dryad.b7156 under "1112 Empirical UCE Nexus Files".
This file should be extracted inside this folder before this script is run.

## make_xml.py

Makes BEAST XML files with _Crocidura_ sequences, based on every XML file in
the `templates` folder. It also generates a `start.sh` script for every
replicate used to begin all chains in that folder. `convert_to_fasta.py` must
be run first.

## run_raijin.py

This script was used to run chains on the NCI Raijin cluster, which uses the
PBS scheduling software. It runs a chain for a given length and checks ESS
values. If any value is below 200, it doubles the length of the chain and
checkes again.

## site_models.csv

List of site models (e.g. JC69, HKY, GTR) for each locus.