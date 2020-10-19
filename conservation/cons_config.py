#Setup output directory for differential expression downloads and outputs

import os

####################################
###                              ###
### Assign variable for datetime ###
###                              ###
####################################

from datetime import datetime

today = datetime.now()
timestring = today.strftime("%Y_%M_%D")
now = today.strftime('%Y_%m_%d')

###############################
###			    ###
### Make output directories ###
###			    ###
###############################

working_path = "/home/hmhamilt/phd_dissertation/conservation/"
if not (os.path.exists(working_path)) :
	os.mkdir(working_path)

path_downloads = os.path.join(working_path,'downloads')
if not (os.path.exists(path_downloads)) :
	os.mkdir(path_downloads)

oma_path = os.path.join(path_downloads,'oma')
if not (os.path.exists(oma_path)) :
	os.mkdir(oma_path)

oma_now = os.path.join(oma_path,now)
if not (os.path.exists(oma_now)) :
        os.mkdir(oma_now)

glygen_path = os.path.join(path_downloads,'glygen')
if not (os.path.exists(glygen_path)) :
        os.mkdir(glygen_path)

glygen_now = os.path.join(glygen_path,now)
if not (os.path.exists(glygen_now)) :
        os.mkdir(glygen_now)

uniprot_path = os.path.join(path_downloads,'uniprot')
if not (os.path.exists(uniprot_path)) :
        os.mkdir(uniprot_path)

uniprot_now = os.path.join(uniprot_path,now)
if not (os.path.exists(uniprot_now)) :
        os.mkdir(uniprot_now)

bgee_path = os.path.join(path_downloads,'bgee')
if not (os.path.exists(bgee_path)) :
        os.mkdir(bgee_path)

bgee_now = os.path.join(bgee_path,now)
if not (os.path.exists(bgee_now)) :
        os.mkdir(bgee_now)

oncomx_path = os.path.join(path_downloads,'oncomx')
if not (os.path.exists(oncomx_path)) :
        os.mkdir(oncomx_path)

oncomx_now = os.path.join(oncomx_path,now)
if not (os.path.exists(oncomx_now)) :
        os.mkdir(oncomx_now)

path_intermediate = os.path.join(working_path, 'intermediate')
if not (os.path.exists(path_intermediate)) :
	os.mkdir(path_intermediate)

path_generated = os.path.join(working_path, 'generated')
if not (os.path.exists(path_generated)) :
	os.mkdir(path_generated)
