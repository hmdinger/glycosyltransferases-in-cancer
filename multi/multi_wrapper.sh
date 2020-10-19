#!/bin/bash

#This is a wrapper around scripts to analyze multimodal data for gts

#Step 1 - Configure paths and working directories

python multi_config.py

#Step 2 - Download all inputs

python download.py

#Step 3 - Filter files appropriatel

python filter.py

# Step 4 - miRNA expression

python mirna_expression.py

#Step 5 - scRNA expression

python scrna_expression.py

#Step 6 - Literature

python lit.py

#Stp 7 - Glycan

python glycans.py

#Step 8 - Residues

sh residue.sh


