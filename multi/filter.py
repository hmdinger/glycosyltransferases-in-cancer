#This program prints a dictionary of miRNAs and single protein targets, and then filters the output for GTs. 

import sys
import urllib
import os
import multi_config
import ssl
import csv
from os.path import basename
import pandas as pd

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

os.chdir(path_intermediate)

##########################
###		       ###
### Print current task ###
###		       ###
##########################

print('Downloading mouse and ortholog data')

##############################
###			   ###
### Retrieve file from URL ###
###			   ###
##############################

#Rename files
os.system('cp NIHMS1511243-supplement-9.csv single_miRNAs.csv')
os.system('cp NIHMS1511243-supplement-11.csv eight_miRNAs.csv')
os.system('cp NIHMS1511243-supplement-12.csv targets.csv')

#Define variables

cancer_gts = []
with open(path_diff_gen + '/temp.tsv','r') as csv_file:
	gts = csv_file.readlines()
        for gt in gts:
        	cancer_gts.append(gt.strip())

print('GTs: ' + str(len(cancer_gts)))

##############################
###    miRNA filtering    ####
##############################

cancer_mirnas = []
target_pairs = []
#output_rows = []
#header = ['mirna','target_gene']
#output_rows.append(header)
comma = ","
with open('targets.csv', 'r') as targets_file:
	lines = csv.reader(targets_file)	
	next(lines)
	for line in lines:
		#output_row.append(line[0])
		if comma in line[2]:
			targets = line[2]
			targets = targets.split(',')
			for x in range(len(targets)):
				target_pair = []
				target_pair.append(line[0])
				target_pair.append(targets[x])
				if target_pair in target_pairs:
					continue
				else:
					target_pairs.append(target_pair)
		else:
			target_pair = []
			target_pair.append(line[0])
			target_pair.append(line[2])
			if target_pair in target_pairs:	
				continue
			else:
				target_pairs.append(target_pair)

#with open('target_pairs.tsv', 'w') as out_file:
#	writer = csv.writer(out_file, delimiter = "\t")
#	writer.writerows(output_rows)

print('There are ' + str(len(target_pairs)) + ' unique miRNA/target pairs.')

output_rows = []
header = ['mirna','target_gene']
output_rows.append(header)

for pair in target_pairs:
	ac = pair[1]
	for gt in cancer_gts:
		if ac == gt:
			print(pair)
			output_rows.append(pair)
		else:
			continue

print('There are ' + str(len(output_rows)-1) + ' unique miRNA/gt target pairs.')

with open('gt_target_pairs.tsv', 'w') as out_file:
	writer = csv.writer(out_file, delimiter='\t')
        writer.writerows(output_rows)

##############################
###    scRNA filtering     ###
##############################

copy = """cp """ + oncomx_now + """/human_cancer_scRNA_preferential_expression.csv scRNA_pref_exp.csv"""
os.system(copy)

output_rows = []
#Filter for rows containing gts and for lung cancer
print('Filtering scRNA-seq data')
with open('scRNA_pref_exp.csv', "r") as sc_file:
	rows = csv.reader(sc_file)
	i = 0
	for row in rows:
		if i == 0:
			output_rows.append(row)
		else:
			gt = row[1].replace('"','')
			cancer = row[3].replace('"','')
			if (gt in cancer_gts) and (cancer == 'Lung cancer'):
				output_rows.append(row)
			else:
				continue
		i += 1

with open ('scrna_filtered.csv', "w") as out_file:
	writer = csv.writer(out_file)
	writer.writerows(output_rows)

#Print completion of current task
print('Done')


