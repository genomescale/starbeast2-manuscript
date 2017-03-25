#!/bin/bash
for i in `seq 0 29`;
do
	python measure_trees.py $i
done
