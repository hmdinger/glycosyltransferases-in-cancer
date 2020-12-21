# GT list generation pipeline
#### Last uploaded: September 6, 2020 9:53AM

The following steps are each accomplished by a single script, which are joined by a wrapper to execute all. More details for each step can be found commented in the script for that step.

NOTE: *** denotes some manual curation process was applied to that step

## CONFIGURATION
1. Configure and create working and download directories. From this directory run: 
	- `python gt_list_config.py` 
	
	This will create a set of working directories in `/working_dir_`

## DOWNLOADS
2. Retrieve UniProt proteins with keyword 0328
	- `uniprot_kw_0328.py`

3. Retrieve CAZy enzymes
	- `uniprot_cazy.py`

4. Retrieve CFG enzymes and map entries to canonical UniProtKB Ac
	- `cfg_table_parse.py`
	- `cfg_map.py`

5. Retrieve GO enzymes involved in sugar transfer using custom slim under 0016757 and reorder to make UniProtKB Ac first column***
	- Manually curate to exclude transfer to/from non-protein substrates
	- `go_custom_slim_0016757.py`
	- `go_reorder.sh`

6. Retrieve InterPro domains for search term "glycosylatransferase" and retrieve corresponding proteins
	- `interpro_glycosyltransferase.py`

## COMPILE
7. Create a list of accessions from all sources, and produce a list of the unique ones
	- `compile_uniq.py`


## GENERATE FINAL TABLE
1. Retrieve list of enzymes by UniProt keyword
  - Columns should be: `uniprotkb_ac | protein_name | gene_name/gene_symbol | uniprotkb_kw`

2. Retrieve list of enzymes corresponding to slim
	- Columns should be uniprotkb_ac | go | go_eco
3. Retrieve CAZY list
	- Columns should be uniprotkb_ac | cazy_family
4. Retrieve CFG list
	- Columns should be uniprotkb_ac | cfg_id
5. Retrieve InterPro families
	- Columns should be uniprotkb_ac | interpro_families
	- NOTE: may require curation
6. Retrieve Pfam domains
	- Columns should be uniprotkb_ac | pfam_domain
	- NOTE: may require curation
7. Retrieve Kelley Moremen's list
	- Check GlyGen
8. Retrieve Ensembl human GTs
9. Retrieve Ensembl mouse GTs
10. Retrieve 1:1 orthologs from OMA


## JOIN TABLES

Output and interim files will be hosted at /home/hmhamilt/phd_dissertation/gt_list_generation/generated 

1. Join all tables using uniprotkb_ac as the key
2. May need to manually add reactions


## MANUAL CURATION

Final table will be hosted at /home/hmhamilt/phd_dissertation/gt_list_generation/generated/reviewed

1. Review by Jeet/GlyGen
2. Review by Will/Kelley?
3. Review by Nathan?


CHART GENERATION

Charts and tables based on this information were generated for the publication using excel, ppt, R, etc.
