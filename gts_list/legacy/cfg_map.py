#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This script parses uniprot accessions from the reference details of a
glycosyltransferase molecule page at functional glycomics.org.
"""
import re
import os
import csv
import urllib.request
from datetime import datetime
import bs4 as bs
import gt_list_config

#Assign variable for datetime
today = datetime.now()
now = today.strftime('%Y_%m_%d')

# Set working directory to downloads path
path_downloads = gt_list_config.path_downloads
path_intermediate = gt_list_config.path_intermediate
cfg_path = gt_list_config.cfg_path
cfg_now = gt_list_config.cfg_now

print(cfg_now)

os.chdir(cfg_now)

# Print current task
print('Retrieving protein accessions from functionalglycomics.org')

# Parsing html from website
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
    url = str('http://www.functionalglycomics.org/glycomics/molecule/' \
                    + 'jsp/glycoEnzyme/viewGlycoEnzyme.jsp?sideMenu=true&' \
                    + 'gbpId={}&pageType=reference').format(cfg_name)
    print('\n'+url)
    source = urllib.request.urlopen(url).read()
    soup = bs.BeautifulSoup(source,'lxml')

    output_row = []
    output_row1 = []
    output_row2 = []

    def niceprot(href):
        """
        NiceProt is a tool that provides a user-friendly tabular view of
        SWISS-PROT entries. The 'NiceProt View of SWISS-PROT' is accessible
        from the top and  bottom of  each  SWISS-PROT entry on ExPASy. You can
        use this tool  to link  to any  SWISS-PROT by using the following style
        of URL: http://www.expasy.ch/cgi-bin/niceprot.pl?P01585 (where the last
        part of the URL is a valid primary accession number).
        """

        return href and re.compile("niceprot").search(href)

    link = soup.find_all(href=niceprot)
    # ex: <a href="http://us.expasy.org/cgi-bin/niceprot.pl?Q9NPC4" target="blank">Q9NPC4</a>

    link_string = str(link).strip('[]')

    id_tmp = re.sub('<.*?>', '', link_string)

    if "," in id_tmp:
        id1 = id_tmp.split(', ')[0]
        id2 = id_tmp.split(', ')[1]
        output_row1.append(id1)
        output_row1.append(cfg_name)
        output_row2.append(id2)
        output_row2.append(cfg_name)
        output_rows.append(output_row1)
        output_rows.append(output_row2)

    elif id_tmp == "":
        continue

    else:
        id = id_tmp
        output_row.append(id)
        output_row.append(cfg_name)
        output_rows.append(output_row)

# os.chdir(path_intermediate)

with open(path_intermediate + '/cfg_map.tsv', 'wb') as csvfile :
    writer = csv.writer(csvfile, delimiter = '\t')
    writer.writerows(output_rows)

#Print completion of current task
print('Done')
