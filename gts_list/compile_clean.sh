#!/bin/bash

#This is a simple script to remove noncanonical entries from the list of unique UniProt KB accessions.

path_intermediate=$(python -c "import gt_list_config,os; path_intermediate = gt_list_config.path_intermediate; print(path_intermediate)")

cd $path_intermediate

#This step no longer needed since noncanonical accessions coming from GO were fixed in the updated GO reordered script
#sed 's/-.*$//g' unique_gts.tsv > unique_gts_canonical.tsv 

sed -e 's/[\r\n]$//g' unique_gts.tsv > unique_gts_formatteid.tsv

#This step shouldn't be needed due to sorted(set()) used in compile_uniq.py, but leaving it in shouldn't hurt unless the input grows substantially.
sort -u unique_gts_formatted.tsv > unique_gts_final.tsv
