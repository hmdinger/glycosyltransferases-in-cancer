#Setup output directory for GT list generation downloads

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

working_path = "/home/hmhamilt/phd_dissertation/gts_list/"
if not (os.path.exists(working_path)) :
	os.mkdir(working_path)

path_downloads = os.path.join(working_path,'downloads')
if not (os.path.exists(path_downloads)) :
	os.mkdir(path_downloads)

uniprot_path = os.path.join(path_downloads,'uniprot')
if not (os.path.exists(uniprot_path)) :
        os.mkdir(uniprot_path)

keyword_path = os.path.join(uniprot_path,'kw_0328')
if not (os.path.exists(keyword_path)) :
        os.mkdir(keyword_path)

keyword_now = os.path.join(keyword_path,now)
if not (os.path.exists(keyword_now)) :
        os.mkdir(keyword_now)

all_path = os.path.join(uniprot_path,'all')
if not (os.path.exists(all_path)) :
        os.mkdir(all_path)

all_now = os.path.join(all_path,now)
if not (os.path.exists(all_now)) :
        os.mkdir(all_now)

cazy_path = os.path.join(path_downloads,'cazy')
if not (os.path.exists(cazy_path)) :
        os.mkdir(cazy_path)

cazy_now = os.path.join(cazy_path,now)
if not (os.path.exists(cazy_now)) :
        os.mkdir(cazy_now)

cfg_path = os.path.join(path_downloads,'cfg')
if not (os.path.exists(cfg_path)) :
        os.mkdir(cfg_path)

cfg_now = os.path.join(cfg_path,now)
if not (os.path.exists(cfg_now)) :
        os.mkdir(cfg_now)

go_path = os.path.join(path_downloads,'go')
if not (os.path.exists(go_path)) :
        os.mkdir(go_path)

go_now = os.path.join(go_path,now)
if not (os.path.exists(go_now)) :
        os.mkdir(go_now)

interpro_path = os.path.join(path_downloads,'interpro')
if not (os.path.exists(interpro_path)) :
        os.mkdir(interpro_path)

interpro_now = os.path.join(interpro_path,now)
if not (os.path.exists(interpro_now)) :
        os.mkdir(interpro_now)

path_intermediate = os.path.join(working_path, 'intermediate')
if not (os.path.exists(path_intermediate)) :
	os.mkdir(path_intermediate)

path_generated = os.path.join(working_path, 'generated')
if not (os.path.exists(path_generated)) :
	os.mkdir(path_generated)

#out_dir='/home/hmhamilt/phd_dissertation/glycosyltransferase_list/downloads/' + today.strftime('%Y_%m_%d')
#os.mkdir(out_dir)

