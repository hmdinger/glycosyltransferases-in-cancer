#This program downloads the table of human UniProtKB accessions with CAZy glycosyltransferase annotations coming from a reviewed complete proteome. Columns retrieved include UniProtKB AC, (Preferred) Gene name, Protein names, CAZy Family name

import sys
import urllib
import os
import gt_list_config

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

cazy_path = gt_list_config.cazy_path

cazy_now = gt_list_config.cazy_now

#cazy_path = os.path.join(path_downloads,'cazy')
#if not (os.path.exists(cazy_path)) :
#        os.mkdir(cazy_path)

#cazy_now = os.path.join(cazy_path,now)
#if not (os.path.exists(cazy_now)) :
#        os.mkdir(cazy_now)

os.chdir(cazy_now)

##########################
###                    ###
### Print current task ###
###                    ###
##########################

print('Downloading UniProtKB table for proteins with CAZy GT annotations to' + cazy_now)

##############################
###                        ###
### Retrieve file from URL ###
###                        ###
##############################

#NOTE - deprecated due to obsolete keyword complete proteome - urllib.urlretrieve('https://www.uniprot.org/uniprot/?query=organism:"Homo+sapiens+(Human)+[9606]"+keyword:"Complete+proteome+[KW-0181]"+reviewed:yes+database:(type:cazy glycosyltransferase)&columns=id,genes(PREFERRED),protein names,database(CAZy)&format=tab', filename='/home/hmhamilt/phd_dissertation/downloads/uniprot_cazy_' + today.strftime('%Y_%m_%d') + '.tab')

urllib.urlretrieve ('https://www.uniprot.org/uniprot/?query=organism:9606+reviewed:yes+proteome:UP000005640+database:(type:cazy glycosyltransferase)&columns=id,genes(PREFERRED),protein names,database(CAZy)&format=tab', filename="uniprot_cazy_proteome.tab")


#Print completion of current task
print('Done')


