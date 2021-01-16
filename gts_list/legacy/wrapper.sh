#!/bin/bash

#This is a wrapper around the gt list scripts. 

#Step 1 - Configuration of working directories

python gt_list_config.py

#Step 2 - Retrieve UniProt enzymes with KW-0328

python uniprot_kw_0328.py

#Step 3 - Retrieve enzymes with CAZy annotations in UniProt

python uniprot_cazy.py

#Step 4 - Retrieve enzyme list from CFG

python cfg_table_parse.py

#Step 5 - Retrieve UniProt accessions from CFG html

python cfg_map.py

#Step 6 - Retrieve go slim of relevant terms under 0016757

python go_custom_slim_0016757.py

#Step 7 - Cut first column of go output

sh go_reorder.sh

#Step 8 - Retrieve ids and corresponding enzymes from InterPro using search term "glycosyltransferase"

python interpro_glycosyltransferase.py

#Step 9 - Clean up interpro formattings

sh interpro_clean.sh

#Step 10 - Compile unique list of GTs from all sources

python compile_uniq.py

#Step 11 - Clean up GTs list

sh compile_clean.sh

#Step 12 - Retrieve master table from UniProt

python uniprot_all_gts.py

#Step 13 - Final clean up

sh final_clean.sh

#Step 14 - Append additional evidence (GO, CFG, InterPro)

python append.py
