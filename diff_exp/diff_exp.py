#!/usr/bin/env python3 
# -*- coding: utf-8 -*-

################################################################################
                        ##Differnetial Expression##
"""
The following steps are each accomplished by a single script, which are joined 
by a wrapper to execute all. More details for each step can be found commented 
in the script for that step.

NOTE: *** denotes some manual curation process was applied to that step
"""
################################################################################

import os
import sys
import urllib.request
import os
import ssl
from os.path import basename
import bs4 as bs
import csv
from pathlib import Path
from datetime import datetime

# Assign variable for datetime
NOW = datetime.now().strftime('%Y_%m_%d')

# Build directory strings
PROJECT_DIRECTORY = os.path.abspath(os.path.join(
                    os.path.dirname(__file__), os.path.pardir))

WORKING_PATH = PROJECT_DIRECTORY + '/working_dir/diff_exp'
DOWNLOADS = os.path.join(WORKING_PATH, NOW, 'downloads')
# Make output directories

#______________________________________________________________________________#
def configure_paths():
    """
    Configure and create working and download directories.
    PROJECT_DIRECTORY is assigned by the root directory and then each
    directory is built from that.
    This will create a set of working directories in `/working_dir_`
    """

    directory_list = [WORKING_PATH, DOWNLOADS]

    # Create directories from list if they don't exist
    for directory in directory_list:
        Path(directory).mkdir(parents=True, exist_ok=True)
#______________________________________________________________________________#
def bioXpress():
    """
    This program downloads the BioXpress v3 and v4 results for differential 
    expression in cancer and the list of gts from GlyGen.
    """

    print('Downloading BioXpress tables and Glycosyltransferase list')

    v2 = 'https://hive.biochemistry.gwu.edu/' \
        +'beta/bioxpress/content/BioXpress_interface_overall_final_v2.0.csv'
    v4_study = 'https://data.oncomx.org/ln2wwwdata' \
        + '/reviewed/human_cancer_mRNA_expression_per_study.csv'
    v4_tissue = 'https://data.oncomx.org/ln2wwwdata' \
        + '/reviewed/human_cancer_mRNA_expression_per_tissue.csv'
    gts = 'https://data.glygen.org/ln2wwwdata/reviewed/' \
        + 'human_protein_glycosyltransferase.csv'

    bioX = os.path.join(DOWNLOADS, 
                    'BioXpress_interface_overall_final_v2.0.csv')
    study = os.path.join(DOWNLOADS, 
                    'human_cancer_mRNA_expression_per_study.csv')
    tissue = os.path.join(DOWNLOADS, 
                    'human_cancer_mRNA_expression_per_tissue.csv')
    hp_gts = os.path.join(DOWNLOADS, 
                    'human_protein_glycosyltransferase.csv')

    # Use tuples to link name(url[1]) and URL(url[0])
    urls = [(v2,bioX), (v4_study, study), (v4_tissue, tissue), (gts, hp_gts)]

    # NOTE - this is a temporary fix and should check on the SSL certificate
    # issues
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE

    for url in urls :
        with urllib.request.urlopen (url[0], context=ctx)  as response:
            with open(url[1], 'wb') as file:
                req = response.read()
                file.write(req)
        print(basename(url[0]), ' downloaded')

    #Print completion of current task
    print('Done')

#______________________________________________________________________________#
def define_vars():
    """
    Define lists for GTs and cancers
    """

    # path_downloads = diff_exp_config.path_downloads
    #
    # path_intermediate = diff_exp_config.path_intermediate
    #
    # bioxpress_now = diff_exp_config.bioxpress_now
    #
    # glygen_now = diff_exp_config.glygen_now
    #
    # os.chdir(path_downloads)
    # Print current task
    print('Parsing tables')

    # Extracting gts from tables
    print('Defining list of GTs')

    hp_gts = os.path.join(DOWNLOADS, 
                    'human_protein_glycosyltransferase.csv')
    bioX = os.path.join(DOWNLOADS, 
                    'BioXpress_interface_overall_final_v2.0.csv')
    gts_tsv = os.path.join(DOWNLOADS, 'gts.tsv')
    cancers_tsv = os.path.join(DOWNLOADS, 'cancers.tsv')

    output_rows = []
    with open(hp_gts, "r") as gt_list:
        gt_csv = csv.reader(gt_list)
        next(gt_csv)
        gts = []
        for row in gt_csv :
            if row[0] in gts :
                continue
            else:
                gts.append(row[0].split('-', 1)[0])
        print(len(gts))

    output_rows = []

    with open(bioX, "r") as input:
        input_csv = csv.reader(input)
        next(input_csv)
        cancers = []
        for row in input_csv :
            if row[8] in cancers:
                continue
            else:
                cancers.append(row[8])
        print(cancers)

    for cancer in cancers :
            output_row = []
            output_row.append(cancer)
            output_rows.append(output_row)

    with open(cancers_tsv, 'w') as csvfile :
            writer = csv.writer(csvfile, delimiter = '\t')
            writer.writerows(output_rows)

    #Print completion of current task
    print('Done')

#______________________________________________________________________________#
def parse_expression():
"""
    
"""
#______________________________________________________________________________#
def main():
    """
    Main
    """

    # configure_paths()
    # bioXpress()
    define_vars()

#______________________________________________________________________________#
if __name__ == "__main__":
    main()