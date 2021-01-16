#This program downloads the table of human UniProtKB accessions for all proteins with evidence of glycosyltransferase activity compiled across resources

import sys
import urllib
import os
import gt_list_config
import glob
import csv

#Assign variable for datetime
from datetime import datetime
today = datetime.now()
now = today.strftime('%Y_%m_%d')

###############################################
###			  		    ###
### Set working directory to downloads path ###
###			  		    ###
###############################################

path_downloads = gt_list_config.path_downloads

path_intermediate = gt_list_config.path_intermediate

path_generated = gt_list_config.path_generated

all_path = gt_list_config.all_path

all_now = gt_list_config.all_now

os.chdir(path_intermediate)

##########################
###		       ###
### Print current task ###
###		       ###
##########################

print('Downloading UniProtKB table for all proteins with glycosyltransferase evidence across resources')

##############################
###			   ###
### Retrieve file from URL ###
###			   ###
##############################

#for file in glob.glob("*"):
#	print(file)

output_rows = []

with open('unique_gts_final.tsv', 'r') as uniprot_acs:
	os.chdir(all_now)
	for uniprot_ac in uniprot_acs:
		print(uniprot_ac)
		uniprot_ac = uniprot_ac.rstrip('\n')
		url = "https://www.uniprot.org/uniprot/?query=organism:9606+AND+reviewed:yes+AND+accession:" + uniprot_ac + "&format=tab&columns=id,genes(PREFERRED),protein names,ec,database(CAZy)"
		urllib.urlretrieve (url, filename="uniprot_" + uniprot_ac + ".tab")
		with open("uniprot_" + uniprot_ac + ".tab", 'r') as active_file:
			i = 0
			output_row = []
			for line in active_file:
				line = line.rstrip('\n')
				i+=1
				if i==1:
					continue
				else:
					output_row.append(line)
					output_rows.append(output_row)

os.chdir(path_intermediate)

with open('gts_with_evidence.tsv', 'wb') as newfile:
	writer = csv.writer(newfile, delimiter = '\t')
        writer.writerows(output_rows)

#Print completion of current task
print('Done')


