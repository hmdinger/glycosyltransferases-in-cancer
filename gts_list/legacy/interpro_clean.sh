#!/bin/bash

#This is a simple script to clean the extra quotations out of the interpro_proteins file.

interpro_now=$(python -c "import gt_list_config,os; interpro_now = gt_list_config.interpro_now; print(interpro_now)")
path_intermediate=$(python -c "import gt_list_config,os; path_intermediate = gt_list_config.path_intermediate; print(path_intermediate)")

cd $interpro_now

sed 's/"""//g' interpro_proteins.tsv > $path_intermediate/interpro_formatted.tsv

cd $path_intermediate

grep -v accession interpro_formatted.tsv > interpro_partial_clean.tsv

grep -v timed interpro_partial_clean.tsv > interpro_clean.tsv

header="accession	length	name	source_database	source_organism	id"

echo $header > header.txt

cat header.txt interpro_clean.tsv > interpro_final.tsv

