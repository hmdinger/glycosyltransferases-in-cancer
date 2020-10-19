#!/bin/bash

#This is a wrapper around the differential conservation scripts. 

#Step 1 - Configuration of working directories

python cons_config.py

#Step 2 - Download inputs from OMA and GlyGen

python download.py

#Step 3 - Filter OMA list for 1:1 orthologs

python filter.py

#Step 4 - Create table of ortholog pairs

python orthologs.py

#Step 5 - Create table of expression conservation

python conservation.py

