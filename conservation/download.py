#This program downloads the BioXpress v3 and v4 results for differential expression in cancer and the list of gts from GlyGen.

import sys
import urllib
import os
import cons_config
import ssl
from os.path import basename

#Assign variable for datetime
from datetime import datetime
today = datetime.now()
now = today.strftime('%Y_%m_%d')

###############################################
###			  		    ###
### Set working directory to downloads path ###
###			  		    ###
###############################################

path_downloads = cons_config.path_downloads

oma_now = cons_config.oma_now

glygen_now = cons_config.glygen_now

bgee_now = cons_config.bgee_now

oncomx_now = cons_config.oncomx_now

path_intermediate = cons_config.path_intermediate

os.chdir(glygen_now)

##########################
###		       ###
### Print current task ###
###		       ###
##########################

print('Downloading mouse and ortholog data')

##############################
###			   ###
### Retrieve file from URL ###
###			   ###
##############################

#v2 = 'https://hive.biochemistry.gwu.edu/beta/bioxpress/content/BioXpress_interface_overall_final_v2.0.csv'
#v4_study = 'https://data.oncomx.org/ln2wwwdata/reviewed/human_cancer_mRNA_expression_per_study.csv'
#v4_tissue = 'https://data.oncomx.org/ln2wwwdata/reviewed/human_cancer_mRNA_expression_per_tissue.csv'
mouse_gts = 'https://data.glygen.org/ln2wwwdata/reviewed/GLY_000030.csv'
#oma = 'https://omabrowswer.org/api/pairs/9606/10090/?rel_type=1%3A1'
#This file is HUGE and contains all species pairs of all relations, not a good way to organize the information. oma = 'https://omabrowser.org/All/oma-pairs.txt.gz'
oma = 'https://omabrowser.org/cgi-bin/gateway.pl?f=PairwiseOrthologs&p1=HUMAN&p2=MOUSE&p3=SwissProt_AC'
oncomx_bgee_hs = 'https://data.oncomx.org/ln2wwwdata/reviewed/human_normal_expression.csv'
oncomx_bgee_mm = 'https://data.oncomx.org/ln2wwwdata/reviewed/mouse_normal_expression.csv'
oncomx_bgee_hs_mapped = 'https://data.oncomx.org/ln2wwwdata/reviewed/human_normal_expression_custom.csv'
oncomx_bgee_mm_mapped = 'https://data.oncomx.org/ln2wwwdata/reviewed/mouse_normal_expression_custom.csv'

urls = [mouse_gts,oma,oncomx_bgee_hs,oncomx_bgee_mm,oncomx_bgee_hs_mapped,oncomx_bgee_mm_mapped]

#NOTE - this is a temporary fix and should check on the SSL certificate issues
ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE

for url in urls :
	a = basename(url)
	if url == mouse_gts :
		os.chdir(glygen_now)
		print('Retrieving data from ' + url)
		urllib.urlretrieve (url, context=ctx, filename=a)
	elif url == oma :
		os.chdir(oma_now)
		print('Retrieving data from ' + url)
                urllib.urlretrieve (url, context=ctx, filename='PairwiseOrthologs.txt')
	else:
		os.chdir(oncomx_now)
		print('Retrieving data from ' + url)
		urllib.urlretrieve (url, context=ctx, filename=a)

#gunzip oma-pairs.txt.gz

os.system('cp /home/hmhamilt/phd_dissertation/diff_exp/intermediate/gts_clean.tsv ' + path_intermediate + '/')

#Print completion of current task
print('Done')


