#This program downloads the table of human UniProtKB accessions with "Glycosyltransferase" [KW-0328] annotations coming from a reviewed complete proteome. Columns retrieved include UniProtKB AC, (Preferred) Gene name, Protein names.

import sys
import urllib
import os
import gt_list_config

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

#keyword_path = gt_list_config.keyword_path

keyword_now = gt_list_config.keyword_now

os.chdir(keyword_now)

##########################
###		       ###
### Print current task ###
###		       ###
##########################

print('Downloading UniProtKB table for proteins with "Glycosyltransferase" [KW-0328] annotations to ' + keyword_now)

##############################
###			   ###
### Retrieve file from URL ###
###			   ###
##############################

#NOTE: deprecated URL as of late 2019 due to obsolete kw complete proteome - urllib.urlretrieve('https://www.uniprot.org/uniprot/?query=organism:"Homo+sapiens+(Human)+[9606]"+keyword:"Complete+proteome+[KW-0181]"+reviewed:yes+keyword:"Glycosyltransferase+[KW-0328]"&columns=id,genes(PREFERRED),protein names&format=tab', filename="uniprot_kw_0328" + today.strftime('%Y_%m_%d') + ".tab")

urllib.urlretrieve ('https://www.uniprot.org/uniprot/?query=organism:9606+AND+reviewed:yes+AND+keyword:"Glycosyltransferase [KW-0328]"+AND+proteome:UP000005640&format=tab', filename="uniprot_kw_0328_proteome.tab")

#Print completion of current task
print('Done')


