#This program downloads miRNA target genes from PMID: 30384176 and OncoMX, scRNA-seq summaries from OncoMX, a custom literature mining dataset and the corresponding updated data from OncoMX, and glycan GlycoTree data from GlyGen. 

import sys
import urllib
import os
import multi_config
import ssl
from os.path import basename
import pandas as pd

#Assign variable for datetime
from datetime import datetime
today = datetime.now()
now = today.strftime('%Y_%m_%d')

###############################################
###			  		    ###
### Set working directory to downloads path ###
###			  		    ###
###############################################

path_downloads = multi_config.path_downloads

pmc_now = multi_config.pmc_now

glygen_now = multi_config.glygen_now

oncomx_now = multi_config.oncomx_now

hgnc_now = multi_config.hgnc_now

path_intermediate = multi_config.path_intermediate

os.chdir(path_downloads)

##########################
###		       ###
### Print current task ###
###		       ###
##########################

print('Downloading data')

##############################
###			   ###
### Retrieve file from URL ###
###			   ###
##############################

#List links in format variable = 'url'
sdemirna_single_cancer = 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6279243/bin/NIHMS1511243-supplement-9.xls'
sdemirna_eight = 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6279243/bin/NIHMS1511243-supplement-11.xls'
sdemirna_targets = 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6279243/bin/NIHMS1511243-supplement-12.xls'
bioxpress_mirna_v2 = 'https://data.oncomx.org/ln2wwwdata/reviewed/human_cancer_mRNA_expression_per_tissue.csv'
bioxpress_mirna_v3 = 'https://data.oncomx.org/ln2wwwdata/reviewed/human_cancer_miRNA_expression.csv'
bioxpress_mirna_v4 = 'https://data.oncomx.org/ln2wwwdata/reviewed/human_cancer_mRNA_expression_per_study.csv'
sc = 'https://data.oncomx.org/ln2wwwdata/reviewed/human_cancer_scRNA_preferential_expression.csv'
dexter = 'https://data.oncomx.org/ln2wwwdata/reviewed/human_cancer_expression_literature.csv'
#other literature
glycans ='https://data.glygen.org/ln2wwwdata/reviewed/GLY_000284.csv'
mirna_map = 'https://www.genenames.org/cgi-bin/genegroup/download?id=476&type=branch'

#FDA
#EDRN
#neoepitope
#clinical
#glycotyper

urls = [bioxpress_mirna_v2, bioxpress_mirna_v3, bioxpress_mirna_v4, sc, dexter, glycans, mirna_map]
# sdemirna_single_cancer, sdemirna_eight, sdemirna_targets]

#NOTE - this is a temporary fix and should check on the SSL certificate issues
#ctx = ssl.create_default_context()
#ctx.check_hostname = False
#ctx.verify_mode = ssl.CERT_NONE
#ctx = ssl.SSLContext()
ssl._create_default_https_context = ssl._create_unverified_context


for url in urls :
	a = basename(url)
	substring = "NIH"
	print(url)
	print(a)
	if url == glycans :
		os.chdir(glygen_now)
		print('Retrieving data from ' + url)
                #urllib.urlretrieve (url, context=ctx, filename=a)
		urllib.urlretrieve (url, filename=a)
	elif url == mirna_map :
                os.chdir(hgnc_now)
                print('Retrieving data from ' + url)
                #urllib.urlretrieve (url, context=ctx, filename='mirna_accessions.tsv')
                urllib.urlretrieve (url, filename='mirna_accessions.tsv')
	elif substring in a :
                os.chdir(pmc_now)
                print('Retrieving data from ' + url)
                #urllib.urlretrieve (url, context=ctx, filename=a)
                urllib.urlretrieve (url, filename=a)
                read_file = pd.read_excel (a)
                a = a.rstrip('.xls')
                os.chdir(path_intermediate)
                read_file.to_csv (a + '.csv', index = None, header = True)
	else:
                os.chdir(oncomx_now)
                print('Retrieving data from ' + url)
                #urllib.urlretrieve (url, context=ctx, filename=a)
                urllib.urlretrieve (url, filename=a)

os.system('cp /home/hmhamilt/phd_dissertation/diff_exp/intermediate/gts_clean.tsv ' + path_intermediate + '/')

#Print completion of current task
print('Done')


