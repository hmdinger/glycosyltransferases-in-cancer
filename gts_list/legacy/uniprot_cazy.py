#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This program downloads the table of human UniProtKB accessions with CAZy
glycosyltransferase annotations coming from a reviewed complete proteome.
Columns retrieved include UniProtKB AC, (Preferred) Gene name, Protein names,
CAZy Family name
"""

import os
import urllib.request
from datetime import datetime
import gt_list_config

#Assign variable for datetime
today = datetime.now()
now = today.strftime('%Y_%m_%d')

# Set working directory to downloads path
path_downloads = gt_list_config.path_downloads

cazy_path = gt_list_config.cazy_path

cazy_now = gt_list_config.cazy_now

os.chdir(cazy_now)

# Print current task
print('Downloading UniProtKB table for proteins with CAZy GT annotations to' + cazy_now)


# Retrieve file from URL
uniprot = 'https://www.uniprot.org/uniprot/?query=' \
                + 'organism:9606+reviewed:yes+proteome:UP000005640' \
                + '+database:(type:cazy%20glycosyltransferase)&columns=id' \
                + ',genes(PREFERRED),protein%20names,database(CAZy)&format=tab'\

with urllib.request.urlopen (uniprot)  as response:
    with open('uniprot_cazy_proteome.tsv', 'wb') as file:
        req = response.read()
        file.write(req)

#Print completion of current task
print('Done')
