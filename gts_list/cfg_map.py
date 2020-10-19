#This script parses uniprot accessions from the reference details of a glycosyltransferase molecule page at functional glycomics.org.

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

cfg_path = gt_list_config.cfg_path

cfg_now = gt_list_config.cfg_now

#cfg_path = os.path.join(path_downloads,'cfg')
#if not (os.path.exists(cfg_path)) :
#        os.mkdir(cfg_path)

#cfg_now = os.path.join(cfg_path,now)
#if not (os.path.exists(cfg_now)) :
#        os.mkdir(cfg_now)

os.chdir(cfg_now)

##########################
###                    ###
### Print current task ###
###                    ###
##########################

print('Retrieving protein accessions from functionalglycomics.org')

#################################
###                           ###
### Parsing html from website ###
###                           ###
#################################

with open('cfg.tsv', 'r') as csv_file:
    lines = csv_file.readlines()

cfg_names = []
output_rows = []
header = []

i = 0
for line in lines :
    i += 1
    if i == 1 :
	header.append('uniprotkb_ac')
        header.append('cfg_id')
        output_rows.append(header)
    else : 
	data = line.split('\t')
        cfg_names.append(data[0])

cfg_names = sorted(set(cfg_names))
print(cfg_names)

for cfg_name in cfg_names :
    print(cfg_name)
    url = "http://www.functionalglycomics.org/glycomics/molecule/jsp/glycoEnzyme/viewGlycoEnzyme.jsp?sideMenu=true&gbpId={}&pageType=reference".format(cfg_name)
    source = urllib.urlopen(url).read()
    soup = bs.BeautifulSoup(source,'lxml')

    output_row = []
    output_row1 = []
    output_row2 = []

    def niceprot(href):
       return href and re.compile("niceprot").search(href)
    
    link = soup.find_all(href=niceprot)
    # ex: <a href="http://us.expasy.org/cgi-bin/niceprot.pl?Q9NPC4" target="blank">Q9NPC4</a>

    link_string = str(link).strip('[]')

    id_tmp = re.sub('<.*?>', '', link_string)
    
    if "," in id_tmp :
	id1 = id_tmp.split(', ')[0]
        id2 = id_tmp.split(', ')[1]
	output_row1.append(id1)
        output_row1.append(cfg_name)
        output_row2.append(id2)
        output_row2.append(cfg_name)
        output_rows.append(output_row1)
        output_rows.append(output_row2)
    
    elif id_tmp == "" :
	continue   

    else: 
	id = id_tmp

    	output_row.append(id)
    	output_row.append(cfg_name)
        output_rows.append(output_row)

os.chdir(path_intermediate)

with open('cfg_map.tsv', 'wb') as csvfile :
    writer = csv.writer(csvfile, delimiter = '\t')
    writer.writerows(output_rows)

#Print completion of current task
print('Done')

