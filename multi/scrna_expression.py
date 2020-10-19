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
with open('scrna_filtered.csv', "r") as gt_file:
        gt_rows = csv.reader(gt_file, delimiter="\t")
        next(gt_rows)
        gts = []
        for row in gt_rows:
                if row[0] in gts :
                        continue
                else:
                        gts.append(row[0])

print('There are ' + str(len(gts)) + ' DEGTs that are also present in the scRNA-seq preference measure analysis.')

for gt in gts :
        output_row = []
        output_row.append(gt)
        output_rows.append(output_row)

with open('scrna_gts.tsv', 'wb') as csvfile :
        writer = csv.writer(csvfile, delimiter = '\t')
        writer.writerows(output_rows)

#Make master scRNA expression table
print('Generating scRNA expression table')

de_file = pd.read_csv(path_diff_gen + '/de_gts_master_filtered.tsv', sep="\t")
sc_file = pd.read_csv(path_intermediate + '/scrna_filtered.csv', sep=",")

print(de_file.head())
print(sc_file.head())

#Define cancer-tissue map from unique tissues in bgee file
print('Making a dictionary for cancer terms and their long names')
map = {
        "LUAD": "Lung cancer",
        "LUSC": "Lung cancer",
}
print(map)
print(map["LUAD"])
#print(map.values())

de_file['cancer_name'] = de_file['TCGA Cancer'].map(map)

common1 = pd.merge(de_file, sc_file, how='left', left_on=['UniProtKB_AC','cancer_name'], right_on=['uniprotkb_ac','do_name'])
print(common1.head())

print(len(de_file.index))
print(len(common1.index))

#common1_null = common1.loc[common1['mirna'].isnull()]
#print(common1_null.head())
#print(len(common1_null.index))

common1_full = common1.dropna(subset=['uniprotkb_ac'])
print(common1_full.head())
print(len(common1_full.index))

for col in common1_full.columns:
	print(col)

better_headers = common1_full.rename(columns = {'Gene': 'gene', 'TCGA Cancer': 'cancer', 'Trend': 'cancer_trend', 'adjusted p_value': 'adjusted_p_value', 'log2FoldChange': 'log2fc'})

print(better_headers.head())

better_headers.to_csv(path_generated + '/gt_scrna.tsv', sep='\t', columns=['uniprotkb_ac', 'gene', 'cancer', 'cancer_trend', 'adjusted_p_value', 'log2fc', 'cell_type', 'pem_score', 'gene_preference', 'is_marker_gene'], index=False)

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

print('Generating subset of GT/miRNA pairs with similar expression trends')
output_rows = []
with open('gt_scrna.tsv', "r") as master_file:
	rows = csv.reader(master_file, delimiter="\t")
	i = 0
	for row in rows:
		if i == 0:
			output_rows.append(row)
		else:
			if (row[3] == 'Down') and (row[8] == 'low'):
				output_rows.append(row)
			elif (row[3] == 'Up') and (row[8] == 'high'):
				output_rows.append(row)
			else:
				continue
		i += 1

with open('similar_gt_scrna.tsv', "w") as out_file:
	writer = csv.writer(out_file, delimiter = "\t") 
	writer.writerows(output_rows)

print('Filtering for increased foldchange to identify high-value pairs')
thresholds = [1.0, 2.0, 3.0]
for threshold in thresholds:
	output_rows = []
	with open('similar_gt_scrna.tsv', "r") as similar_file:
        	sims = csv.reader(similar_file, delimiter="\t")
	        needheader=True
   		for row in sims:
                	if needheader:
                        	needheader=False
	                        output_rows.append(row)
        	                continue
                	elif (((float(row[5])) >= threshold) == True) or (((float(row[5])) <= -threshold) == True) :
                        	output_rows.append(row)
	                else:
        	                continue
	with open("high_value_gt_scrna_similar_pairs_" + str(threshold) + ".tsv", "w") as out_file:
        	writer = csv.writer(out_file, delimiter='\t')
	        writer.writerows(output_rows)

quit()

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

