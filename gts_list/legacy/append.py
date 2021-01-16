#This script appends evidence columns to compiled master table.

import sys
import urllib
import os
import gt_list_config
import bs4 as bs
import csv
import re
import pandas as pd

#Assign variable for datetime
from datetime import datetime
today = datetime.now()
now = today.strftime('%Y_%m_%d')

###############################################
###                                         ###
### Set working directory to downloads path ###
###                                         ###
###############################################

working_path = gt_list_config.working_path

path_downloads = gt_list_config.path_downloads

path_intermediate = gt_list_config.path_intermediate

path_generated = gt_list_config.path_generated

os.chdir(path_intermediate)

##########################
###                    ###
### Print current task ###
###                    ###
##########################

print('Appending additional annotations')

#################################
###                           ###
### Appending annotations     ###
###                           ###
#################################

#os.system('cp gts_master_final.tsv gts_evidence_append.tsv')

#os.system('cp go_reordered3.tsv ' + path_intermediate + '/go_reordered_copy.tsv')

#os.system("sed 's/GENE\ PRODUCT\ ID/entry/' go_reordered_copy.tsv > go_reordered_dictionary.tsv")

#os.system("sed 's/GO\ TERM/go/' go_reordered_copy.tsv > go_reordered_dictionary.tsv")

print('GO')

with open('go_canonical.tsv', "r") as go_file_handle:
        go_file_csv = csv.reader(go_file_handle, delimiter="\t")
        next(go_file_csv) # Skip the header
        # Lets record all of the go to map info as a dictionary of list entries
        go_data = {}
        for row in go_file_csv:
            if row[0] in go_data:
                # Already there, additional GO term to add
                go_data[row[0]].append(row[3])
            else:
                # Not already entered, lets create it!
                go_data[row[0]] = []
                go_data[row[0]].append(row[3])
	

#	for row in csv.DictReader(go_file, delimiter='\t'):
#		ac = row['GENE PRODUCT ID']
#		go[ac] = row['GO TERM']
		
#print(go_data)

with open('gts_append_go.tsv', 'w') as output_file_handle:
        # Lets make it a CSV file by using the built in CSV WRITER functionality instead of doing this by hand
        output_file = csv.writer(output_file_handle, delimiter="\t")
        # Open the GTS file.  We will read line by line and append any of the GO TERMS found for
        # that Accession number as a new last column of the file.
        with open('gts_master_final.tsv', 'r') as gts_file_handle:
            gts_file_csv = csv.reader(gts_file_handle, delimiter="\t")
            # We will want to get the header to re-write back into the output file
            need_header = True
            for row in gts_file_csv:
                if need_header:
                    # Haven't printed the header yet, so lets do that
                    need_header = False
                    # Write out header to the output file.
                    row.append("GO_TERM")
                    output_file.writerow(row)
                    continue
                # Map to the UniprotKB Accession number
                if row[0] in go_data:
                    # The AC is found in the go data!
                    # Get terms in semi-colon separated form
                    go_terms = ";".join(go_data[row[0]])
                    # Add to row as new column
                    row.append(go_terms)
                    # Write out to file
                    output_file.writerow(row)
                else:
                    # No GO Term data, add empty column and move on.
                    row.append("")
                    output_file.writerow(row)

#df=pd.read_csv('gts_evidence_append.tsv', sep='\t', skiprows=1, names=["entry","gene","protein","ec","cazy"])
#print(df["entry"])
#df["go"]=df.entry.apply(lambda x:go[x])
#df['go']=df['entry'].map(go)
#df.to_csv('gts_master_append_go.csv', sep='\t', index=False)
#print(df)

#Print completion of current task
print('Done')

