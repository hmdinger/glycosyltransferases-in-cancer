#!/bin/bash

#This is a script to summarize residue impact

short_file="/home/hmhamilt/phd_dissertation/multi/generated/short_3.0.tsv"
int_path="/home/hmhamilt/phd_dissertation/multi/intermediate"
gen_path="/home/hmhamilt/phd_dissertation/multi/generated"

cd ${gen_path}

pwd 

parents=$(awk -F'\t' '{print $9}' ${short_file} | sort -u)
echo "${parents[*]}"

for parent_id in $parents; do
	echo ${parent_id}
	residues=$(grep -w ${parent_id} ${short_file})
	echo ${residues[*]}
	for residue in $residues; do
		uniprot=$(${residue[1]})
		echo ${uniprot}

	done
done			



echo 'Done' 
