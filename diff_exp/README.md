# README for GT list generation pipeline

### Last updated: September 10, 2020 10:25PM
**************************************

The following steps are each accomplished by a single script, which are joined by a wrapper to execute all. More details for each step can be found commented in the script for that step.

NOTE: *** denotes some manual curation process was applied to that step

## CONFIGURATION
1. Configure and create working and download directories 
	- diff_exp_config.py 

## DOWNLOADS
2. Retrieve BioXpress v3 and v4

## COMPILE
6. Create a list of accessions from all sources, and produce a list of the unique ones
	- compile_uniq.py


### GENERATE FINAL TABLE
Retrieve list of enzymes by UniProt keyword
	- Columns should be uniprotkb_ac | protein_name | gene_name/gene_symbol | uniprotkb_kw
	
3. Retrieve list of enzymes corresponding to slim
	- Columns should be uniprotkb_ac | go | go_eco
	
4. Retrieve CAZY list
	- Columns should be uniprotkb_ac | cazy_family

5. Retrieve CFG list
	- Columns should be uniprotkb_ac | cfg_id

6. Retrieve InterPro families
	- Columns should be uniprotkb_ac | interpro_families
	- NOTE: may require curation

7. Retrieve Pfam domains
	- Columns should be uniprotkb_ac | pfam_domain
	- NOTE: may require curation

8. Retrieve Kelley Moremen's list
	- Check GlyGen

9. Retrieve Ensembl human GTs

10. Retrieve Ensembl mouse GTs

11. Retrieve 1:1 orthologs from OMA


### JOIN TABLES

Output and interim files will be hosted at 
/home/hmhamilt/phd_dissertation/gt_list_generation/generated 

1. Join all tables using uniprotkb_ac as the key
2. May need to manually add reactions


### MANUAL CURATION

Final table will be hosted at /home/hmhamilt/phd_dissertation/gt_list_generation/generated/reviewed

1. Review by Jeet/GlyGen
2. Review by Will/Kelley?
3. Review by Nathan?


### CHART GENERATION

Charts and tables based on this information were generated for the publication using excel, ppt, R, etc.
