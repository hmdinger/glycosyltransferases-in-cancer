#This program downloads the table of human GO annotations for the parent "GO:0016757 transferase activity, transferring glycosyl groups" and relevant child terms. Columns retrieved are TBD.

import sys
import urllib
import os
import gt_list_config
import bs4 as bs
import csv
import re
import requests

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

go_path = gt_list_config.go_path

go_now = gt_list_config.go_now

#go_path = os.path.join(path_downloads,'go')
#if not (os.path.exists(go_path)) :
#        os.mkdir(go_path)

#go_now = os.path.join(go_path,now)
#if not (os.path.exists(go_now)) :
#        os.mkdir(go_now)

os.chdir(go_now)

##########################
###                    ###
### Print current task ###
###                    ###
##########################

print('Downloading GO table for slim built from "GO:0016757 transferase activity, transferring glycosyl groups"')

##############################
###                        ###
### Retrieve file from URL ###
###                        ###
##############################

#Example URL from web interface search: https://www.ebi.ac.uk/QuickGO/annotations?goId=GO:0016758,GO:0016763,GO:0008373,GO:0097363&taxonId=9606&taxonUsage=exact&geneProductSubset=Swiss-Prot&qualifier=enables&proteome=gcrpCan,complete&geneProductType=protein

#Base API URL https://www.ebi.ac.uk/QuickGO/api/index.html

#NOTE - python command was generated from API interface

#urllib.urlretrieve ('"https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?geneProductType=protein&geneProductSubset=Swiss-Prot&goId=GO%3A0016758%2CGO%3A0016763%2CGO%3A0008373%2CGO%3A0097363&taxonId=9606&qualifier=enables"')

requestURL = "https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?geneProductType=protein&geneProductSubset=Swiss-Prot&goId=GO%3A0016758%2CGO%3A0016763%2CGO%3A0008373%2CGO%3A0097363&taxonId=9606&qualifier=enables"

r = requests.get(requestURL, headers={ "Accept" : "text/tsv"})

if not r.ok:
  r.raise_for_status()
  sys.exit()

responseBody = r.text

#output_rows = []

#for line in responseBody :
#    data = line.split('\t')
#    output_rows.append(data)


#print(responseBody)
with open('go.tsv', 'wb') as f :
    f.write(responseBody)

#Print completion of current task
print('Done')


