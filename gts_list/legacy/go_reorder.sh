#!/bin/bash

#This is a simple script to remove the first column of the go slim output such that UniProt accession is in the first column and use only canonical accessions.

echo "Reordering GO columns and generating canonical accessions only"

go_now=$(python -c "import gt_list_config,os; go_now = gt_list_config.go_now; print(go_now)")
path_intermediate=$(python -c "import gt_list_config,os; path_intermediate = gt_list_config.path_intermediate; print(path_intermediate)")

cd $go_now

pwd

cut -d$'\t' -f2- go.tsv > $path_intermediate/go_reordered.tsv

cd $path_intermediate

#awk -F'\t' '{print $1}' go_reordered2.tsv | grep - | head

##awk -F'\t' '{sub(/-.*$/,"",$1); print $0}' go_reordered2.tsv > go_reordered3.tsv

awk -F'\t' '{sub(/-.*$/,"",$1)}1' OFS='\t' go_reordered.tsv > go_canonical.tsv

#awk -F'\t' '{print $1}' go_reordered3.tsv | grep - | head 

echo "Done"
