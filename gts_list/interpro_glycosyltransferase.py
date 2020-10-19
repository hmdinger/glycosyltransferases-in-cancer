#This program downloads the table of InterPro entries using search term "glycosyltransferase."

import sys
import urllib
import os
import gt_list_config
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

interpro_path = gt_list_config.interpro_path

interpro_now = gt_list_config.interpro_now

#interpro_path = os.path.join(path_downloads,'interpro')
#if not (os.path.exists(interpro_path)) :
#        os.mkdir(interpro_path)

#interpro_now = os.path.join(interpro_path,now)
#if not (os.path.exists(interpro_now)) :
#	os.mkdir(interpro_now)

os.chdir(interpro_now)

##########################
###		       ###
### Print current task ###
###		       ###
##########################

print('Downloading InterPro table for entries with "glycosyltransferase" search hits to ' + interpro_now)

##############################
###			   ###
### Retrieve file from URL ###
###			   ###
##############################

#NOTE - a few options for searches using EBI API and InterPro API, none ideal so switching to two parts: 1) retrieve id list for query, 2) retrieve proteins from id list. List will likely contain other species, so filter for human in step 2. 
#urllib.urlretrieve ('https://www.ebi.ac.uk/api/entry/interpro/search/text/glycosyltransferase', filename="interpro_search.tab")
#urllib.urlretrieve ('https://www.ebi.ac.uk/ebisearch/ws/rest/interpro7?query=glycosyltransferase%20AND%20(source_database:interpro%5E2%20OR%20*:*)&format=tsv&fields=id,description,name,source_database&start=0&size=100', filename="interpro_search.tab")
#urllib.urlretrieve ('https://www.ebi.ac.uk/ebisearch_k8s/ws/rest/interpro?query=glycosyltransferase&format=tsv&fields=id,entry,description,name,source_database', filename="interpro_search.tab")

#EBI search for retrieval of ID list from interpro for search term "glycosyltransferase"
urllib.urlretrieve ('https://www.ebi.ac.uk/ebisearch_k8s/ws/rest/interpro?query=glycosyltransferase&size=1000&format=idlist', filename="interpro_search.txt")

print ('Retrieving proteins for hit InterPro IDs to ' + interpro_now)

#InterPro retrieval of human proteins generated from list of IDs retrieved in previous step

#Loop through IDs

good_lines = []

with open('interpro_search.txt') as gt_ids :
	for gt_id in gt_ids :
		maybe_lines = []
		gt_id = gt_id.rstrip("\n")
		print (gt_id)		
		url = "https://www.ebi.ac.uk/interpro/api/protein/reviewed/entry/interpro/" + gt_id + "/taxonomy/uniprot/9606?format=tsv"
		print (url)
		urllib.urlretrieve (url, filename="interpro_interim.tsv")

		#for temp_file in file_list :
        	filesize = os.path.getsize("interpro_interim.tsv")
	
		if filesize == 0:
			continue

		else: 
        		with open('interpro_interim.tsv', 'r') as test_file:
                		lines = test_file.readlines()
                        	for line in lines:
                                	data = []
					line = line.rstrip("\n")
					line = line.rstrip("\r")
					data = line.split("\t")
					data.append(gt_id)
					print (data)
					maybe_lines = []
					for datum in data:		
						maybe_lines.append(datum)
					good_lines.append(maybe_lines)

with open('interpro_proteins.tsv', 'wb') as good_file:
        writer = csv.writer(good_file, delimiter = '\t')
        writer.writerows(good_lines)

#Print completion of current task
print('Done')


