#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This program downloads the table of human UniProtKB accessions with
"Glycosyltransferase" [KW-0328] annotations coming from a reviewed complete
proteome. Columns retrieved include UniProtKB AC, (Preferred) Gene name,
Protein names.
"""

import urllib.request
import os
from datetime import datetime
import gt_list_config

#Assign variable for datetime
today = datetime.now()
now = today.strftime('%Y_%m_%d')

# Set working directory to downloads path ###
path_downloads = gt_list_config.path_downloads

# keyword_path = gt_list_config.keyword_path
keyword_now = gt_list_config.keyword_now

os.chdir(keyword_now)

# Print current task ###
print('Downloading UniProtKB table annotations to ' + keyword_now)

# Retrieve file from URL ###
uniprot = 'https://www.uniprot.org/uniprot/?query=' \
                + 'organism:9606+AND+reviewed:yes+AND+' \
                + 'keyword:%22Glycosyltransferase%20[KW-0328]%22+AND+' \
                + 'proteome:UP000005640&format=tab'


with urllib.request.urlopen (uniprot)  as response:
    with open('uniprot_kw_0328_proteome.tsv', 'wb') as file:
        req = response.read()
        file.write(req)

#Print completion of current task
print('Done')
