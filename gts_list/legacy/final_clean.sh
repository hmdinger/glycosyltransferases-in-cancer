#!/bin/bash

#This is a simple script to clean the extra quotations from the master gts and details file.

path_intermediate=$(python -c "import gt_list_config,os; path_intermediate = gt_list_config.path_intermediate; print(path_intermediate)")

cd $path_intermediate

sed 's/"//g' gts_with_evidence.tsv > gts_partial_clean.tsv

sed -e 's/[\r\n]$//g' gts_partial_clean.tsv > gts_formatted.tsv

header="entry\tgene\tprotein\tec\tcazy"

echo -e $header > header.txt

cat header.txt gts_formatted.tsv > gts_master_final.tsv

