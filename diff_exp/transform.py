#This script transforms data from the master output differential expression file for compatibility with visual tools.

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

print('Transforming tables')

##################################
###			       ###
### Extracting gts from tables ###
###			       ###
##################################


#EXAMPLE ROW:
#UniProtKB_AC    RefSeq  Gene    log2FoldChange  p_value adjusted p_value        Significant     Trend   TCGA Cancer     Cancer Ontology #Patients       Data Source    PMID     UBERON_ID       PtsTrend        TotPts  Percent
#A0PJZ3  NP_001073862.1  GXYLT2  -1.0    1.23e-02        3.19e-02        Yes     Down    BLCA    DOID:11054 / Urinary bladder cancer [UBC]       12/19(63.16)    RNASeqV2-       UBERON:0001255  12      19      63.16

print('Reordering patient counts')
output_rows=[]

with open('Differential_expression_gts_master.tsv', "r") as de_file:
	de_file_csv = csv.reader(de_file, delimiter="\t")
#	next(de_file_csv)
	patient_data = {}
	i=0
	for row in de_file_csv:
		if i == 0:
			#row[10]=row[10].replace("#Patients","Patients_trend, Patients_opp, Patients_percent") 
			row.insert(15, 'PtsOpp')
			row.append('PtsUp')
			row.append('PtsDown')
			order = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 13, 14, 18, 19, 16, 17]
			row = [row[j] for j in order]
			output_rows.append(row)
		else:
			#row[10]=row[10].replace(")","")
			#row[10]=re.split('\/|\(', row[10])
			opp=int(row[15])-int(row[14])
			opp=str(opp)
			row.insert(15, opp)
			if row[7] == "Up":
				row.append(row[14])
				row.append(row[15])
			else:
				row.append(row[15])
				row.append(row[14])
			order = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 13, 14, 18, 19, 16, 17]
                        row = [row[j] for j in order]
			output_rows.append(row)
		i+=1	

#Example with new rows - [0'UniProtKB_AC', 1'RefSeq', 2'Gene', 3'log2FoldChange', 4'p_value', 5'adjusted p_value', 6'Significant', 7'Trend', 8'TCGA Cancer', 9'Cancer Ontology', '10#Patients', 11'Data Source', 12'PMID', 13'UBERON_ID', 14'PtsTrend ', 15'PtsOpp', 16'TotPts ', 17'Percent', 18'PtsUp', 19'PtsDown']
['A0PJZ3', 'NP_001073862.1', 'GXYLT2', '-1.0', '1.23e-02', '3.19e-02', 'Yes', 'Down', 'BLCA', 'DOID:11054 / Urinary bladder cancer [UBC]', '12/19(63.16)', 'RNASeqV2', '-', 'UBERON:0001255', '12', '7', '19', '63.16', '7', '12']

with open("de_gts_master_reordered.tsv", "w") as out_file:
	writer = csv.writer(out_file, delimiter='\t')
	writer.writerows(output_rows)


#Print completion of current task
print('Done')

