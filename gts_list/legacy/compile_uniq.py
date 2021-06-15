#This script generates a unique list of UniProt accessions from across all sources.

import sys
import urllib
import os
import gt_list_config
import bs4 as bs
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

path_downloads = gt_list_config.path_downloads
path_intermediate = gt_list_config.path_intermediate
path_generated = gt_list_config.path_generated
keyword_now = gt_list_config.keyword_now
cazy_now = gt_list_config.cazy_now
cfg_now = gt_list_config.cfg_now
go_now = gt_list_config.go_now
interpro_now = gt_list_config.interpro_now

os.chdir(path_downloads)

##########################
###                    ###
### Print current task ###
###                    ###
##########################

print('Compiling list of unique protein accessions from all sources')

keyword = keyword_now + '/uniprot_kw_0328_proteome.tab'
cazy = cazy_now + '/uniprot_cazy_proteome.tab'
cfg = path_intermediate + '/cfg_map.tsv'
go = path_intermediate + '/go_canonical.tsv'
interpro = path_intermediate + '/interpro_final.tsv'

annotation_files = [keyword, cazy, cfg, go, interpro]

total_gts = []
output_rows = []

for afile in annotation_files :
	print(afile)
	with open(afile, 'r') as tabfile:
    		lines = tabfile.readlines()
        i = 0
	for line in lines :
    	        i +=1
		if i == 1 :
			continue
		else : 
			data = line.split('\t')
			total_gts.append(data[0])

total_gts = sorted(set(total_gts))

for gt in total_gts :
	output_row = []
	output_row.append(gt)
	output_rows.append(output_row)

os.chdir(path_intermediate)

with open('unique_gts.tsv', 'wb') as csvfile :
	writer = csv.writer(csvfile, delimiter = '\t')
	writer.writerows(output_rows)

print('Done')
