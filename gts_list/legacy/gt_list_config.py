#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Configure and create working and download directories.
PROJECT_DIRECTORY is assigned by the root directory and then each directory is 
built from that. 
This will create a set of working directories in `/working_dir_`
"""

import os
from pathlib import Path
from datetime import datetime

# Assign variable for datetime
NOW = datetime.now().strftime('%Y_%m_%d')

# Build directory strings
PROJECT_DIRECTORY = os.path.abspath(os.path.join(
                    os.path.dirname(__file__), os.path.pardir))
WORKING_PATH = PROJECT_DIRECTORY + '/working_dir/'
DOWNLOADS = os.path.join(WORKING_PATH,'downloads')
UNIPROT = os.path.join(DOWNLOADS,'uniprot')
KEYWORD = os.path.join(UNIPROT,'kw_0328')
keyword_NOW = os.path.join(KEYWORD,NOW)
EVERYTHING = os.path.join(UNIPROT,'all')
all_NOW = os.path.join(EVERYTHING,NOW)
CAZY = os.path.join(DOWNLOADS,'cazy')
cazy_NOW = os.path.join(CAZY,NOW)
CFG = os.path.join(DOWNLOADS,'cfg')
cfg_NOW = os.path.join(CFG,NOW)
GO = os.path.join(DOWNLOADS,'go')
go_NOW = os.path.join(GO,NOW)
INTERPRO = os.path.join(DOWNLOADS,'interpro')
interpro_NOW = os.path.join(INTERPRO,NOW)
path_intermediate = os.path.join(WORKING_PATH, 'intermediate')
path_generated = os.path.join(WORKING_PATH, 'generated')

directory_list = [WORKING_PATH, DOWNLOADS, UNIPROT, KEYWORD, keyword_NOW,
                EVERYTHING,all_NOW, CAZY, cazy_NOW, CFG, cfg_NOW, GO, go_NOW,
                INTERPRO, interpro_NOW, path_intermediate, path_generated]

# Create directories from list if they don't exist
for directory in directory_list:
    Path(directory).mkdir(parents=True, exist_ok=True)
