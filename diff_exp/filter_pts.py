#This script filters data from the master output differential expression for studies with patient number > 10.

import sys
import urllib
import os
import diff_exp_config
import bs4 as bs
import csv
import pandas as pd
import re

#Assign variable for datetime
from datetime import datetime
today = datetime.now()
now = today.strftime('%Y_%m_%d')

###############################################
###                                         ###
### Set working directory to downloads path ###
###                                         ###
###############################################

path_downloads = diff_exp_config.path_downloads

path_intermediate = diff_exp_config.path_intermediate

path_generated = diff_exp_config.path_generated

bioxpress_now = diff_exp_config.bioxpress_now

glygen_now = diff_exp_config.glygen_now

os.chdir(path_generated)

##########################
###                    ###
### Print current task ###
###                    ###
##########################

print('Filtering tables')

##################################
###			       ###
### Extracting gts from tables ###
###			       ###
##################################


#EXAMPLE ROW:
#[0'UniProtKB_AC', 1'RefSeq', 2'Gene', 3'log2FoldChange', 4'p_value', 5'adjusted p_value', 6'Significant', 7'Trend', 8'TCGA Cancer', 9'Cancer Ontology', '10#Patients', 11'Data Source', 12'PMID', 13'UBERON_ID', 14'PtsTrend ', 15'PtsOpp', 16'TotPts ', 17'Percent', 18'PtsUp', 19'PtsDown']
#['A0PJZ3', 'NP_001073862.1', 'GXYLT2', '-1.0', '1.23e-02', '3.19e-02', 'Yes', 'Down', 'BLCA', 'DOID:11054 / Urinary bladder cancer [UBC]', '12/19(63.16)', 'RNASeqV2', '-', 'UBERON:0001255', '12', '7', '19', '63.16', '7', '12']

print('Filtering TCGA studies with patients < 10')
output_rows=[]

with open('de_gts_master_reordered.tsv', "r") as de_file:
	de_file_csv = csv.reader(de_file, delimiter="\t")
	needheader=True
	for row in de_file_csv:
		if needheader:
			needheader=False
			output_rows.append(row)
			continue
		if row[8] == "CESC":
			continue
		elif row[8] == "PAAD":
			continue
		elif row[8] == "READ":
			continue
		elif row[8] == "SARC":
			continue
		elif row[8] == "UCEC":
			continue
		else: 
			output_rows.append(row)

with open("de_gts_master_filtered.tsv", "w") as out_file:
	writer = csv.writer(out_file, delimiter='\t')
	writer.writerows(output_rows)

output_rows=[]
print('Filtering for increased foldchange')
with open('de_gts_master_filtered.tsv', "r") as de_file:
        de_file_csv = csv.reader(de_file, delimiter="\t")
        needheader=True
        for row in de_file_csv:
                if needheader:
                        needheader=False
                        output_rows.append(row)
                        continue
		elif (((float(row[3])) >= 1.0) == True):
			output_rows.append(row)
		elif (((float(row[3])) <= -1.0) == True):
			output_rows.append(row)
		else:
			continue

with open("de_gts_master_fc_filtered.tsv", "w") as out_file:
        writer = csv.writer(out_file, delimiter='\t')
        writer.writerows(output_rows)

output_rows=[]
print('Filtering for higher threshold fold change')
with open('de_gts_master_filtered.tsv', "r") as de_file:
        de_file_csv = csv.reader(de_file, delimiter="\t")
        needheader=True
        for row in de_file_csv:
                if needheader:
                        needheader=False
                        output_rows.append(row)
                        continue
                elif (((float(row[3])) >= 2.0) == True):
                        output_rows.append(row)
                elif (((float(row[3])) <= -2.0) == True):
                        output_rows.append(row)
                else:
                        continue

with open("de_gts_master_fc_filtered_high.tsv", "w") as out_file:
        writer = csv.writer(out_file, delimiter='\t')
        writer.writerows(output_rows)

output_rows=[]
print('Filtering for even higher threshold fold change')
with open('de_gts_master_filtered.tsv', "r") as de_file:
        de_file_csv = csv.reader(de_file, delimiter="\t")
        needheader=True
        for row in de_file_csv:
                if needheader:
                        needheader=False
                        output_rows.append(row)
                        continue
                elif (((float(row[3])) >= 2.5) == True):
                        output_rows.append(row)
                elif (((float(row[3])) <= -2.5) == True):
                        output_rows.append(row)
                else:
                        continue

with open("de_gts_master_fc_filtered_higher.tsv", "w") as out_file:
        writer = csv.writer(out_file, delimiter='\t')
        writer.writerows(output_rows)

output_rows=[]
print('Filtering for yet even higher threshold fold change')
with open('de_gts_master_filtered.tsv', "r") as de_file:
        de_file_csv = csv.reader(de_file, delimiter="\t")
        needheader=True
        for row in de_file_csv:
                if needheader:
                        needheader=False
                        output_rows.append(row)
                        continue
                elif (((float(row[3])) >= 3.0) == True):
                        output_rows.append(row)
                elif (((float(row[3])) <= -3.0) == True):
                        output_rows.append(row)
                else:
                        continue

with open("de_gts_master_fc_filtered_highest.tsv", "w") as out_file:
        writer = csv.writer(out_file, delimiter='\t')
        writer.writerows(output_rows)

#Print completion of current task
print('Done')

