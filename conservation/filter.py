#This program downloads the BioXpress v3 and v4 results for differential expression in cancer and the list of gts from GlyGen.

import sys
import urllib
import os
import cons_config
import ssl
import csv
from os.path import basename

#Assign variable for datetime
from datetime import datetime
today = datetime.now()
now = today.strftime('%Y_%m_%d')

###############################################
###			  		    ###
### Set working directory to downloads path ###
###			  		    ###
###############################################

path_downloads = cons_config.path_downloads

oma_now = cons_config.oma_now

glygen_now = cons_config.glygen_now

bgee_now = cons_config.bgee_now

oncomx_now = cons_config.oncomx_now

path_intermediate = cons_config.path_intermediate

path_diff_gen = '/home/hmhamilt/phd_dissertation/diff_exp/generated'

os.chdir(oma_now)

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


#Define variables

cancer_gts = []
with open(path_diff_gen + '/temp.tsv','r') as csv_file:
	gts = csv_file.readlines()
        for gt in gts:
        	cancer_gts.append(gt.strip())

print('GTs: ' + str(len(cancer_gts)))

os.system('head -n -1 PairwiseOrthologs.txt > ' + path_intermediate + '/oma_pairs.txt')

os.chdir(path_intermediate)

recip_pairs = []
with open('oma_pairs.txt', 'r') as oma_file:
	lines = csv.reader(oma_file, delimiter="\t")	
	for line in lines:
		#print(line)
		if line == "":
			print(line)
			print("skipping this line")		
		elif line[2] == "1:1": 
			recip_pairs.append(line)
		else:
			continue

print(str(len(recip_pairs)))

output_rows = []
header = ['human_uniprotkb_ac', 'mouse_uniprotkb_ac', 'relationship type', 'oma_group_id']
output_rows.append(header)
for pair in recip_pairs:
	ac = pair[0]
	for gt in cancer_gts:
		if ac == gt:
			print(pair)
			output_rows.append(pair)
		else:
			continue

print('GTs with cancer data and 1:1 orthologs: ' + str(len(output_rows)-1))

with open('gt_human_mouse_orthologs.tsv', 'w') as out_file:
	writer = csv.writer(out_file, delimiter='\t')
        writer.writerows(output_rows)

#Print completion of current task
print('Done')


