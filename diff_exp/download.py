#This program downloads the BioXpress v3 and v4 results for differential expression in cancer and the list of gts from GlyGen.

import sys
import urllib
import os
import diff_exp_config
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

path_downloads = diff_exp_config.path_downloads

bioxpress_now = diff_exp_config.bioxpress_now

glygen_now = diff_exp_config.glygen_now

os.chdir(bioxpress_now)

##########################
###		       ###
### Print current task ###
###		       ###
##########################

print('Downloading BioXpress tables and Glycosyltransferase list')

##############################
###			   ###
### Retrieve file from URL ###
###			   ###
##############################

v2 = 'https://hive.biochemistry.gwu.edu/beta/bioxpress/content/BioXpress_interface_overall_final_v2.0.csv'
v4_study = 'https://data.oncomx.org/ln2wwwdata/reviewed/human_cancer_mRNA_expression_per_study.csv'
v4_tissue = 'https://data.oncomx.org/ln2wwwdata/reviewed/human_cancer_mRNA_expression_per_tissue.csv'
gts = 'https://data.glygen.org/ln2wwwdata/reviewed/GLY_000004.csv'

urls = [v2,v4_study,v4_tissue,gts]

os.system('pwd')

#NOTE - this is a temporary fix and should check on the SSL certificate issues
ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE


for url in urls :
	a = basename(url)
	if url == gts :
		os.chdir(glygen_now)
		print('Retrieving data from ' + url)
		urllib.urlretrieve (url, context=ctx, filename=a)
	else :
		print('Retrieving data from ' + url)
                urllib.urlretrieve (url, context=ctx, filename=a)

#Print completion of current task
print('Done')


