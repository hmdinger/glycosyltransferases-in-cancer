#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

This script scrapes the table of cfg glycosyltransferases from the results page
at functional glycomics.org.
"""

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

cfg_path = gt_list_config.cfg_path

cfg_now = gt_list_config.cfg_now

#cfg_path = os.path.join(path_downloads,'cfg')
#if not (os.path.exists(cfg_path)) :
#        os.mkdir(cfg_path)

#cfg_now = os.path.join(cfg_path,now)
#if not (os.path.exists(cfg_now)) :
#        os.mkdir(cfg_now)

os.chdir(cfg_now)


# Print current task
print('Scraping table from functionalglycomics.org')

# Retrieve file from URL
functionalglycomics = 'http://www.functionalglycomics.org/glycomics/' \
            + 'molecule/jsp/glycoEnzyme/gtdb.jsp?species=' \
            + 'Homo+sapiens&linkage_attaching=%3F&linkage_anomeric' \
            + '=%3F&linkage_position=%3F&linkage_base=%3F&pgname=&' \
            + 'from=multiple&title=Multiple+Criteria+Search+Results&' \
            + 'slideNumber=multipleQuery'

with urllib.request.urlopen(functionalglycomics) as response:
    source = response.read()

# Parsing table from website
soup = bs.BeautifulSoup(source,'lxml')
table = soup.find('table', { 'class' : 'glycoEnzymeTablePresentation2' })

output_rows = []
table_rows = table.find_all('tr')
i = 1

for tr in table_rows :
    i +=1
    if i < 5 :
        continue
    columns = tr.find_all('td')
    if len(columns) == 6 :
        output_row = []
        for column in columns:
            output_row.append(column.text)
        output_rows.append(output_row)
    else :
        output_row = []
        output_row.append(output_rows[-1][0])
        output_row.append(output_rows[-1][1])
        for column in columns:
            output_row.append(column.text)
        output_row.append(output_rows[-1][-1])
        output_rows.append(output_row)

with open('cfg.tsv', 'w') as csvfile :
    writer = csv.writer(csvfile, delimiter = '\t')
    writer.writerows(output_rows)

#Print completion of current task
print('Done')