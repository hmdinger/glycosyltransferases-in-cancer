#This program generates a table for with gt differential expression and differential expression for miRNAs targeting those gts 

import sys
import urllib
import os
import multi_config
import glob
import csv
import pandas as pd
import numpy as np

#Assign variable for datetime
from datetime import datetime
today = datetime.now()
now = today.strftime('%Y_%m_%d')

###############################################
###			  		    ###
### Set working directory to downloads path ###
###			  		    ###
###############################################

path_downloads = multi_config.path_downloads

pmc_now = multi_config.pmc_now

glygen_now = multi_config.glygen_now

oncomx_now = multi_config.oncomx_now

hgnc_now = multi_config.hgnc_now

path_intermediate = multi_config.path_intermediate

path_diff_gen = '/home/hmhamilt/phd_dissertation/diff_exp/generated'

path_generated = multi_config.path_generated

os.chdir(path_intermediate)

##########################
###		       ###
### Print current task ###
###		       ###
##########################

print('Downloading UniProtKB gene info for human and mouse and creating ortholog table')

##############################
###			   ###
### Retrieve file from URL ###
###			   ###
##############################

#Define gts from miRNA target pairs
print('Defining GTs')

output_rows = []
with open('gt_target_pairs.tsv', "r") as gt_file:
        gt_rows = csv.reader(gt_file, delimiter="\t")
        next(gt_rows)
        mirna_gts = []
	gt_mirnas = []
        for row in gt_rows:
                if row[1] in mirna_gts :
                        continue
                else:
                        mirna_gts.append(row[1])
		if row[0] in gt_mirnas:
			continue
		else:
			gt_mirnas.append(row[0])

print('There are ' + str(len(mirna_gts)) + ' gts targeted by ' + str(len(gt_mirnas)) + ' miRNAs differentially expressed in cancer.')

for gt in mirna_gts :
        output_row = []
        output_row.append(gt)
        output_rows.append(output_row)

with open('mirna_gts.tsv', 'wb') as csvfile :
        writer = csv.writer(csvfile, delimiter = '\t')
        writer.writerows(output_rows)

#Make miRNA expression file just for the gts
#Skipping this step for now

#print('Filtering miRNA expression data for gts with target info')
#Copying expression file to intermediate path
copy = """cp """ + oncomx_now + """/human_cancer_miRNA_expression.csv """ + path_intermediate + """/human_cancer_miRNA_expression.csv"""
os.system(copy)
copy = """cp """ + hgnc_now + """/mirna_accessions.tsv """ + path_intermediate + """/mirna_accessions.tsv"""
os.system(copy)

#Define cancers
#cancers_clean = []
#with open(path_diff_int + '/cancers_no_manual.tsv','r') as csv_file:
#            cancers = csv_file.readlines()
#            print(cancers)
#            for cancer in cancers:
#                cancers_clean.append(cancer.strip())
#            print(cancers_clean)

#Make master miRNA expression table
print('Generating miRNA expression table')

degts_file = pd.read_csv(path_diff_gen + '/de_gts_master_filtered.tsv', sep="\t")
targets_file = pd.read_csv(path_intermediate + '/gt_target_pairs.tsv', sep="\t")
map_file = pd.read_csv(path_intermediate + '/mirna_accessions.tsv', sep="\t")
mirna_file= pd.read_csv(path_intermediate + '/human_cancer_miRNA_expression.csv', sep=",")

print(degts_file.head())
print(targets_file.head())
print(map_file.head())
print(mirna_file.head())

common1 = pd.merge(degts_file, targets_file, how='left', left_on=['UniProtKB_AC'], right_on=['target_gene'])
print(common1.head())

print(len(degts_file.index))
print(len(common1.index))

#common1_null = common1.loc[common1['mirna'].isnull()]
#print(common1_null.head())
#print(len(common1_null.index))

common1_full = common1.dropna(subset=['mirna'])
print(common1_full.head())
print(len(common1_full.index))

common2 = pd.merge(common1_full, map_file, how='left', left_on=['mirna'], right_on=['Alias symbols'])
print(common2.head())
print(len(common2.index))

common2_full = common2.dropna(subset=['Alias symbols'])
print(common2_full.head())
print(len(common2_full.index))

common3 = pd.merge(common2_full, mirna_file, how='left', left_on=['Approved symbol', 'TCGA Cancer'], right_on=['gene_symbol', 'tcga_cancer_type'])
print(common3.head())
print(len(common3.index))

common3_full = common3.dropna(subset=['gene_symbol'])
print(common3_full.head())
print(len(common3_full.index))

#common1['tissue'] = common1['TCGA Cancer'].map(map)
#print(common1.head())

#for col in common3.columns:
#	print(col)

common3_better_headers = common3_full.rename(columns = {'UniProtKB_AC': 'uniprotkb_ac', 'Gene': 'gene', 'TCGA Cancer': 'cancer', 'Trend': 'gt_cancer_trend', 'adjusted p_value': 'gt_adjusted_p_value', 'log2FoldChange': 'gt_log2fc', 'Approved symbol': 'mirna_symbol', 'adjpvalue': 'mirna_adjusted_p_value', 'log2fc': 'mirna_log2fc', 'expression_change_direction': 'mirna_cancer_trend'})

#score = {
#        "HIGH": 1,
#        "MEDIUM": 0,
#        "LOW": -1,
#        "ABSENT": -10,
#	"Up": -2,
#	"Down": 2
#}

#print(score)
#print(score["ABSENT"])

#common3_better_headers['conservation_score'] = (common3_better_headers['human_healthy_expression'].map(score)) + (common3_better_headers['mouse_healthy_expression'].map(score))

#common3_better_headers['impact_score'] = (common3_better_headers['conservation_score']) + (common3_better_headers['cancer_trend'].map(score))

print(common3_better_headers.head())

common3_better_headers.to_csv(path_generated + '/gt_mirna_pairs_expression.tsv', sep='\t', columns=['uniprotkb_ac', 'gene', 'cancer', 'gt_cancer_trend', 'gt_adjusted_p_value', 'gt_log2fc', 'mirna_symbol', 'mirna_cancer_trend', 'mirna_adjusted_p_value', 'mirna_log2fc'], index=False)

#print('Generating table filtered for absent calls')
#output_rows = []
#with open('conservation_master.tsv', "r") as master_file:
#        rows = csv.reader(master_file, delimiter="\t")
#        i = 0
#        for row in rows:
#                if i == 0:
#                        output_rows.append(row)
#                else:
#                        if int(row[15]) < -4:
#                                continue
#                        else:
#                                output_rows.append(row)
#                i += 1

#with open('no_absent_conservation.tsv', "w") as out_file:
#        writer = csv.writer(out_file, delimiter = "\t")
#        writer.writerows(output_rows)

os.chdir(path_generated)

print('Generating subset of GT/miRNA pairs with opposing expression trends')
output_rows = []
with open('gt_mirna_pairs_expression.tsv', "r") as master_file:
	rows = csv.reader(master_file, delimiter="\t")
	i = 0
	for row in rows:
		if i == 0:
			output_rows.append(row)
		else:
			if (row[3] == 'Down') and (row[7] == 'Up'):
				output_rows.append(row)
			elif (row[3] == 'Up') and (row[7] == 'Down'):
				output_rows.append(row)
			else:
				continue
		i += 1

with open('opposing_gt_mirna_pairs_expression.tsv', "w") as out_file:
	writer = csv.writer(out_file, delimiter = "\t") 
	writer.writerows(output_rows)

print('Filtering for increased foldchange to identify high-value pairs')
thresholds = [1.0, 2.0, 3.0]
for threshold in thresholds:
	output_rows = []
	with open('opposing_gt_mirna_pairs_expression.tsv', "r") as opposites_file:
        	opposites = csv.reader(opposites_file, delimiter="\t")
	        needheader=True
   		for row in opposites:
                	if needheader:
                        	needheader=False
	                        output_rows.append(row)
        	                continue
                	elif  ( ((((float(row[5])) >= threshold) == True) or (((float(row[5])) <= -threshold) == True)) and ((((float(row[9])) >= threshold) == True) or (((float(row[9])) <= -threshold) == True))) :
                        	output_rows.append(row)
	                else:
        	                continue
	with open("high_value_gt_mirna_opposite_pairs_" + str(threshold) + ".tsv", "w") as out_file:
        	writer = csv.writer(out_file, delimiter='\t')
	        writer.writerows(output_rows)

#Summing miRNAs per gt

output_rows = []
header = ['uniprotkb_ac', 'gene', 'number_cancers', 'cancers_same', 'cancers_different', 'number_mirnas', 'mirnas_same', 'mirnas_different']
output_rows.append(header)
#gt = 'Q9H1C3'


check = []
gt_dict = {}
with open('gt_mirna_pairs_expression.tsv', "r") as gt_file:
        gt_rows = csv.DictReader(gt_file, delimiter="\t")
        next(gt_rows)
        for d in gt_rows:
		if d in check:
			continue
		else:	
			check.append(d)
                	gt_dict.setdefault(d['uniprotkb_ac'],d['gene'])


print(gt_dict)
gts = gt_dict.keys()
genes = gt_dict.values()

summary = pd.read_csv('gt_mirna_pairs_expression.tsv', sep="\t")

def same(s):
	if (s['gt_cancer_trend'] == s['mirna_cancer_trend']):
		return 'same'
	else:
		return 'different'

summary['is_same'] = summary.apply(same, axis=1)

for gt in gts:
	output_row = []
	gene = gt_dict[gt]
	if gt in summary.values:
		summary_gt = summary.loc[summary['uniprotkb_ac'] == gt]
		cancer_list = summary_gt['cancer'].unique()
		mirna_list = summary_gt['mirna_symbol'].unique()
		cancer_count = len(cancer_list)
		mirna_count = len(mirna_list)
		summary_same = summary_gt.loc[summary_gt['is_same'] == 'same']
		cancer_same = summary_same['cancer'].unique()
        	cancer_same_count = len(cancer_same)
		mirna_same = summary_same['mirna_symbol'].unique()
		mirna_same_count = len(mirna_same)
		summary_diff = summary_gt.loc[summary_gt['is_same'] == 'different']
                cancer_diff = summary_diff['cancer'].unique()
                cancer_diff_count = len(cancer_diff)
                mirna_diff = summary_diff['mirna_symbol'].unique()
                mirna_diff_count = len(mirna_diff)
		output_row.append(gt)
		output_row.append(gene)
		output_row.append(cancer_count)
		output_row.append(cancer_same_count)
		output_row.append(cancer_diff_count)
		output_row.append(mirna_count)
		output_row.append(mirna_same_count)
		output_row.append(mirna_diff_count)
		output_rows.append(output_row)
	else:
		continue

with open('gt_mirna_summary.tsv', "w") as out_file:
	writer = csv.writer(out_file, delimiter = "\t")
	writer.writerows(output_rows)

diff_rows = []
for row in output_rows:
	if (row[3] == '0') and (row[7] == '0'):
		diff_rows.append(row)

print(diff_rows)
	



print('Done')

