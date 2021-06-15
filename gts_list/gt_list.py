#!/usr/bin/env python3
# -*- coding: utf-8 -*-

################################################################################
                        ##GT List##
"""
    generation of the list of human glycosyltransferases (GTs)
"""
################################################################################

import re
import os
import sys
import csv
import json
import urllib.request

from pathlib import Path
from datetime import datetime

import bs4 as bs

# Assign variable for datetime
NOW = datetime.now().strftime('%Y_%m_%d')

# Build directory strings
PROJECT_DIRECTORY = os.path.abspath(os.path.join(
                    os.path.dirname(__file__), os.path.pardir))
WORKING_PATH = PROJECT_DIRECTORY + '/working_dir/gt_lists'
DOWNLOADS = os.path.join(WORKING_PATH, NOW, 'downloads')

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
def uniprot_kw_0328():
    """
    This program downloads the table of human UniProtKB accessions with
    "Glycosyltransferase" [KW-0328] annotations coming from a reviewed complete
    proteome. Columns retrieved include UniProtKB AC, (Preferred) Gene name,
    Protein names.
    """

    # Print current task ###
    print('Downloading UniProtKB table annotations to ' + DOWNLOADS)

    # Retrieve file from URL ###
    uniprot = 'https://www.uniprot.org/uniprot/?query=' \
                    + 'organism:9606+AND+reviewed:yes+AND+' \
                    + 'keyword:%22Glycosyltransferase%20[KW-0328]%22+AND+' \
                    + 'proteome:UP000005640&format=tab'

    kw_0328 = os.path.join(DOWNLOADS, 'uniprot_kw_0328_proteome.tsv')
    with urllib.request.urlopen (uniprot)  as response:
        with open(kw_0328, 'wb') as file:
            req = response.read()
            file.write(req)

    #Print completion of current task
    print("uniprot_kw_0328_proteome.tsv saved")
#______________________________________________________________________________#
def uniprot_cazy():
    """
    This program downloads the table of human UniProtKB accessions with CAZy
    glycosyltransferase annotations coming from a reviewed complete proteome.
    Columns retrieved include UniProtKB AC, (Preferred) Gene name,
    Protein names, CAZy Family name
    """

    # Print current task
    print('Downloading UniProtKB table for proteins with CAZy GT annotations to'
        + DOWNLOADS)

    # Retrieve file from URL
    uniprot = 'https://www.uniprot.org/uniprot/?query=organism:9606+' \
                    + 'reviewed:yes+proteome:UP000005640+database:' \
                    + '(type:cazy%20glycosyltransferase)&columns=id,genes' \
                    + '(PREFERRED),protein%20names,database(CAZy)&format=tab'\

    cazy_proteome = os.path.join(DOWNLOADS,'uniprot_cazy_proteome.tsv')

    with urllib.request.urlopen (uniprot)  as response:
        with open(cazy_proteome, 'wb') as file:
            req = response.read()
            file.write(req)

    #Print completion of current task
    print('Done')
#______________________________________________________________________________#
def cfg_table_parse():
    """
    This script scrapes the table of cfg glycosyltransferases from the results
    page at functional glycomics.org.
    """

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
    count = 1

    for row in table_rows:
        count +=1
        if count < 5 :
            continue
        output_row = []
        columns = row.find_all('td')
        if len(columns) == 6:
            for column in columns:
                output_row.append(column.text)
            output_rows.append(output_row)
        else :
            output_row.append(output_rows[-1][0])
            output_row.append(output_rows[-1][1])
            for column in columns:
                output_row.append(column.text)
            output_row.append(output_rows[-1][-1])
            output_rows.append(output_row)

    cfg_tsv = os.path.join(DOWNLOADS, 'cfg.tsv')

    with open(cfg_tsv, 'w') as csvfile :
        writer = csv.writer(csvfile, delimiter = '\t')
        writer.writerows(output_rows)

    #Print completion of current task
    print('Done')
#______________________________________________________________________________#
def niceprot(href):
    """
    NiceProt is a tool that provides a user-friendly tabular view of
    SWISS-PROT entries. The 'NiceProt View of SWISS-PROT' is accessible
    from the top and  bottom of  each  SWISS-PROT entry on ExPASy. You can
    use this tool  to link  to any  SWISS-PROT by using the following style
    of URL: http://www.expasy.ch/cgi-bin/niceprot.pl?P01585 (where the last
    part of the URL is a valid primary accession number).

    ex: <a href="http://us.expasy.org/cgi-bin/niceprot.pl?Q9NPC4"
    target="blank">Q9NPC4</a>
    """

    return href and re.compile("niceprot").search(href)
#______________________________________________________________________________#
def cfg_map():
    """
    This script parses uniprot accessions from the reference details of a
    glycosyltransferase molecule page at functional glycomics.org.
    """

    # Print current task
    print('Retrieving protein accessions from functionalglycomics.org')

    cfg = os.path.join(DOWNLOADS, 'cfg.tsv')

    cfg_names = []
    output_rows = [['uniprotkb_ac', 'cfg_id']]

    # reading from CFG_NOW/cfg.tsv
    with open(cfg, 'r') as csv_file:
        next(csv_file) # skip header
        for line in csv_file.readlines() :
            cfg_names.append(line.split('\t')[0])

    cfg_names = sorted(set(cfg_names))

    # Parsing html from website
    for cfg_name in cfg_names:
        url = str('http://www.functionalglycomics.org/glycomics/molecule/' \
                        + 'jsp/glycoEnzyme/viewGlycoEnzyme.jsp?sideMenu=true&' \
                        + 'gbpId={}&pageType=reference').format(cfg_name)
        source = urllib.request.urlopen(url).read()
        soup = bs.BeautifulSoup(source,'lxml')
        link_string = str(soup.find_all(href=niceprot)).strip('[]')
        id_tmp = re.sub('<.*?>', '', link_string)

        if id_tmp == "":
            continue

        if "," in id_tmp:
            output_rows.extend(([id_tmp.split(', ')[0], cfg_name],
                                [id_tmp.split(', ')[1], cfg_name]))
        else:
            output_rows.append([id_tmp, cfg_name])

    cfg = os.path.join(DOWNLOADS, 'cfg_map.tsv')
    with open(cfg, 'w') as csvfile :
        writer = csv.writer(csvfile, delimiter = '\t')
        writer.writerows(output_rows)

    #Print completion of current task
    print('Done')
#______________________________________________________________________________#
def go_custom_slim_0016757():
    """
    This program downloads the table of human GO annotations for the parent
    "GO:0016757 transferase activity, transferring glycosyl groups" and
    relevant child terms. Columns retrieved are TBD.
    """

    # Print current task
    print('Downloading GO table for slim built from "GO:0016757 transferase' \
                                 + 'activity, transferring glycosyl groups"')

    url = 'https://www.ebi.ac.uk/QuickGO/services/annotation/' \
                + 'downloadSearch?geneProductType=protein&geneProductSubset' \
                + '=Swiss-Prot&goId=GO%3A0016758%2CGO%3A0016763%2CGO%3' \
                + 'A0008373%2CGO%3A0097363&taxonId=9606&qualifier=enables'

    req = urllib.request.Request(url)
    req.add_header('Accept', "text/tsv")
    response = urllib.request.urlopen(req)

    print(response.status)
    if response.status != 200:
        response.raise_for_status()
        sys.exit('Response status != 200: OK')

    response_body = response.read()

    go_terms = os.path.join(DOWNLOADS, 'go.tsv')
    with open(go_terms, 'wb') as file:
        file.write(response_body)

    #Print completion of current task
    print('Done')
#______________________________________________________________________________#
def go_reorder():
    """
    This is a simple script to remove the first column of the go slim
    output such that UniProt accession is in the first column and use only
    canonical accessions.
    """

    go_terms = os.path.join(DOWNLOADS, 'go.tsv')
    go_reordered = os.path.join(DOWNLOADS, 'go_canonical.tsv')

    with open(go_terms, 'r') as read:
        csv_reader = csv.reader(read, delimiter = '\t')
        with open(go_reordered, 'w') as write:
            csv_writer = csv.writer(write, delimiter = '\t')
            for row in csv_reader:
                # Sub for (cut -d$'\t' -f2- go.tsv > go_reordered.tsv)
                new_row = row[1:]
                # Sub for (awk -F'\t' '{sub(/-.*$/,"",$1)}1' OFS='\t')
                new_row[0] = new_row[0].split('-')[0]
                csv_writer.writerow(new_row)
#______________________________________________________________________________#
def interpro_glycosyltransferase():
    """
    This program downloads the table of InterPro entries using search term
    "glycosyltransferase."
    """

    # Print current task ###
    print('Downloading InterPro table for entries with ' \
            + '"glycosyltransferase" search hits to ' + DOWNLOADS + '\n')

    # EBI search for ID list from interpro for search term "glycosyltransferase"
    gt_ids_url = 'https://www.ebi.ac.uk/ebisearch/ws/rest/interpro?' \
            + 'query=glycosyltransferase&size=1000&format=idlist'

    gt_ids=[]
    with urllib.request.urlopen(gt_ids_url) as reader:
        # with open(interpro_search, 'wb') as write:
        #     copyfileobj(reader, write)
        for read in reader.readlines():
            gt_ids.append(read.decode('utf-8').rstrip('\n'))
    print ('Retrieving proteins for hit InterPro IDs to ' + DOWNLOADS + '\n')

    # InterPro retrieval of human proteins generated from list of IDs retrieved
    # in previous step

    # Loop through IDs
    gt_id = ''
    url = "https://www.ebi.ac.uk/interpro/api/protein/reviewed/entry/interpro/"\
            + gt_id + "/taxonomy/uniprot/9606?format=json"

    good_lines = []
    for gt_id in gt_ids:
        print (gt_id, '\n')

        with urllib.request.urlopen(url) as query:
            data = query.read()
        encoding = query.info().get_content_charset('utf-8')
        interpro_json = json.loads(data.decode(encoding))
        for result in interpro_json['results']:
            result['metadata'].update({'gt_id':gt_id})
            good_lines.append(result['metadata'])

    interpro_proteins = os.path.join(DOWNLOADS,'interpro_proteins.tsv')
    with open(interpro_proteins, 'w') as good_file:
        writer = csv.DictWriter(good_file, good_lines[0].keys(), delimiter='\t')
        writer.writeheader()
        writer.writerows(good_lines)

  #   #Print completion of current task
    print('Done')
#______________________________________________________________________________#
def compile_unique():
    """
    This script generates a unique list of UniProt accessions from across all
    sources.
    """

    # Print current task
    print('Compiling list of unique protein accessions from all sources')

    keyword = os.path.join(DOWNLOADS, 'uniprot_kw_0328_proteome.tsv')
    cazy = os.path.join(DOWNLOADS, 'uniprot_cazy_proteome.tsv')
    cfg = os.path.join(DOWNLOADS, 'cfg_map.tsv')
    go_terms = os.path.join(DOWNLOADS, 'go_canonical.tsv')
    interpro = os.path.join(DOWNLOADS, 'interpro_proteins.tsv')
    unique_gts = os.path.join(DOWNLOADS, 'unique_gts.tsv')

    annotation_files = [keyword, cazy, cfg, go_terms, interpro]

    total_gts = []

    for afile in annotation_files:
        print('Digesting: \n', afile)
        with open(afile, 'r') as tabfile:
            next(tabfile)
            lines = tabfile.readlines()
            for line in lines:
                total_gts.append(line.split('\t')[0])
    print('Total GTs identified: ', len(total_gts))
    total_gts = sorted(set(total_gts))
    print('Total unique GTs identified: ', len(total_gts))

    with open(unique_gts, 'w') as csvfile :
        writer = csv.writer(csvfile, delimiter = '\t')
        for gts in total_gts:
            writer.writerow([gts])

    print('Done')
#______________________________________________________________________________#
def uniprot_all_gts():
    """
    This program downloads the table of human UniProtKB accessions for all
    proteins with evidence of glycosyltransferase activity compiled across
    resources
    """

    print('Downloading UniProtKB table for all proteins with' \
            + ' glycosyltransferase evidence across resources')

    uniprot_acs = []

    cazy_rows = []
    unique_gts = os.path.join(DOWNLOADS, 'unique_gts.tsv')
    with open(unique_gts, 'r') as file:
        uniprot = csv.reader(file)
        for accession in uniprot:
            uniprot_acs.append(accession[0])

    for uniprot_ac in uniprot_acs:
        print(uniprot_ac)
        url = "https://www.uniprot.org/uniprot/?query=organism:9606+AND" \
            + "+reviewed:yes+AND+accession:" + uniprot_ac + "&format=" \
            + "tab&columns=id,genes(PREFERRED),protein%20names,ec," \
            + "database(CAZy)"
        with urllib.request.urlopen(url) as query:
            data = query.read()
        encoding = query.info().get_content_charset('utf-8')
        if len(cazy_rows) < 1:
            cazy_rows.append(data.decode(encoding).split('\n')[0].split('\t'))
        cazy_rows.append(data.decode(encoding).split('\n')[1].split('\t'))

    gts_with_evidence = os.path.join(DOWNLOADS, 'gts_with_evidence.tsv')

    with open(gts_with_evidence, 'w') as newfile:
        writer = csv.writer(newfile, delimiter = '\t')
        writer.writerows(cazy_rows)

    # #Print completion of current task
    print('Done')
#______________________________________________________________________________#
def append():
    """
    #This script appends evidence columns to compiled master table.
    """

    # Print current task
    print('Appending additional annotations')
    go_canonical = os.path.join(DOWNLOADS, 'go_canonical.tsv')
    gts_master_final = os.path.join(DOWNLOADS, 'gts_with_evidence.tsv')
    gts_append_go = os.path.join(WORKING_PATH, (NOW + '_final_gts_go.tsv'))


    print('GO')
    with open(go_canonical, "r") as go_file_handle:
        go_file_csv = csv.reader(go_file_handle, delimiter="\t")
        next(go_file_csv) # Skip the header
        # record all of the go to map info as a dictionary of list entries
        go_data = {}
        for row in go_file_csv:
            if row[0] in go_data:
                # Already there, additional GO term to add
                go_data[row[0]].append(row[3])
            else:
                # Not already entered, lets create it!
                go_data[row[0]] = [row[3]]
        print('go_data digested')

    header = ['entry','gene','protein','ec','cazy','GO_TERM']
    with open(gts_append_go, 'w') as output_file_handle:
        output_file = csv.writer(output_file_handle, delimiter="\t")
        output_file.writerow(header)# \tprotein\tec\tcazy])
        # Open the GTS file.  We will read line by line and append any of the
        # GO TERMS found for that Accession number as a new last column of the
        # file.
        with open(gts_master_final, 'r') as gts_file_handle:
            gts_file_csv = csv.reader(gts_file_handle, delimiter="\t")
            for row in gts_file_csv:
                # Map to the UniprotKB Accession number
                if row[0] in go_data:
                    # The AC is found in the go data!
                    # Get terms in semi-colon separated form
                    go_terms = ";".join(go_data[row[0]])
                    # Add to row as new column
                    row.append(go_terms)
                    # Write out to file
                    output_file.writerow(row)
                else:
                    # No GO Term data, add empty column and move on.
                    row.append("")
                    output_file.writerow(row)

    #Print completion of current task
    print('Done')
#______________________________________________________________________________#
def main():
    """
    Main
    """

    configure_paths()
    # uniprot_kw_0328()
    # uniprot_cazy()
    # cfg_table_parse()
    # cfg_map()
    # go_custom_slim_0016757()
    # go_reorder()
    # interpro_glycosyltransferase()
    # compile_unique()
    # uniprot_all_gts()
    append()

#______________________________________________________________________________#
if __name__ == "__main__":
    main()
