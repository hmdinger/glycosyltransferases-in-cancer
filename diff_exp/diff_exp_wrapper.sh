#!/bin/bash

#This is a wrapper around the differential expression scripts. 

#Step 1 - Configuration of working directories

python diff_exp_config.py

#Step 2 - Download inputs from BioXpress and GlyGen

python download.py

#Step 3 - Define lists for GTs and cancers

python define_vars.py

#Step 4 - Extract and summarize differential expression from GTs

sh parse.sh

