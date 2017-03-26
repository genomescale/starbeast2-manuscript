# Description of _Pseudacris_ files

These are all scripts used to analyse the computational performance
of StarBEAST2 using _Pseudacris_ data.

## make_xml.py

Makes BEAST XML files with _Pseudacris_ sequences, based on every XML file in
the `templates` folder. It also generates a `start.sh` script for every
replicate used to begin all chains in that folder.

## run_raijin.py

This script was used to run chains on the NCI Raijin cluster, which uses the
PBS scheduling software. It runs a chain for a given length and checks ESS
values. If any value is below 200, it doubles the length of the chain and
checkes again.