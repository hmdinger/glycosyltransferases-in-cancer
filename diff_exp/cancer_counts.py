#This script summarizes cancer counts for GTs.

import sys
import urllib
import os
import diff_exp_config
import csv
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

#Define inputs

filt = 'de_gts_master_filtered.tsv'
fc_filt = 'de_gts_master_fc_filtered.tsv'
fc_high = 'de_gts_master_fc_filtered_high.tsv'
fc_higher = 'de_gts_master_fc_filtered_higher.tsv'
fc_highest = 'de_gts_master_fc_filtered_highest.tsv'

filtered = [filt, fc_filt, fc_high, fc_higher, fc_highest]

#Define variables

for file in filtered:

	#cmd = """awk -F',' '$1 != "-" {print $0}' cancers_clean.tsv > cancers_no_manual.tsv"""
	cmd = """awk -F'\t' 'NR!=1{print $1}' """ + file + """ | sort -u > temp.tsv""" 
	print(cmd)
	os.system(cmd)
	os.system('head temp.tsv')
	
	os.chdir(path_intermediate)
	
	cancers_clean = []
	with open('cancers_no_manual.tsv','r') as csv_file:
	    cancers = csv_file.readlines()
	    print(cancers)
	    for cancer in cancers:
		cancers_clean.append(cancer.strip())
	    print(cancers_clean)

	gts_clean = []
	with open(path_generated + '/temp.tsv','r') as csv_file:
	    gts = csv_file.readlines()
	    for gt in gts:
		gts_clean.append(gt.strip())

	cancer_total = len(cancers_clean)

	print('Cancers ' + str(cancer_total))
	print('GTs ' + str(len(gts_clean)))

	output_rows = []
	header = ['uniprotkb_ac', 'gene', 'cancers_up', 'cancers_down', 'cancers_no_sig']

	output_rows.append(header)

	print('Compiling cancer counts for GTs')

	for gt in gts_clean:
		print(gt)
		output_row = []
		output_row.append(gt)
		with open(path_generated + '/' + file, "r") as gene_file:
			gene_list = csv.reader(gene_file, delimiter="\t")
			for genes in gene_list:
				if genes[0] == gt:
					gene = genes[2]
					#print(gene)
		output_row.append(gene)
		count_up = 0
		count_down = 0
		for cancer in cancers_clean:
		#	print(cancer)
			with open(cancer + '_05_up_list.tsv', "r") as cancer_up_file:
				cancer_up = csv.reader(cancer_up_file, delimiter="\t")
				for line in cancer_up:
					if (gt in line):
						count_up += 1
					else:
						continue
			with open(cancer + '_05_down_list.tsv', "r") as cancer_down_file:
				cancer_down = csv.reader(cancer_down_file, delimiter="\t")
				for line in cancer_down:
					if (gt in line):
						count_down += 1
					else:
						continue
			#print(count_up)
			#print(count_down)
		count_no_sig = cancer_total - (count_up + count_down)
		output_row.append(count_up)
		output_row.append(count_down)
		output_row.append(count_no_sig)
		print(output_row)
		output_rows.append(output_row)

	os.chdir(path_generated)

        file = file.replace("de_gts_master_", "")
        file = file.replace(".tsv", "")

	with open("cancer_counts_by_gt_" + file + ".tsv", "w") as out_file:
		writer = csv.writer(out_file, delimiter='\t')
		writer.writerows(output_rows)

