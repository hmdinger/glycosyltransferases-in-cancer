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

working_path = "/home/hmhamilt/phd_dissertation/diff_exp/"
if not (os.path.exists(working_path)) :
	os.mkdir(working_path)

path_downloads = os.path.join(working_path,'downloads')
if not (os.path.exists(path_downloads)) :
	os.mkdir(path_downloads)

bioxpress_path = os.path.join(path_downloads,'bioxpress')
if not (os.path.exists(bioxpress_path)) :
	os.mkdir(bioxpress_path)

bioxpress_now = os.path.join(bioxpress_path,now)
if not (os.path.exists(bioxpress_now)) :
        os.mkdir(bioxpress_now)

glygen_path = os.path.join(path_downloads,'glygen')
if not (os.path.exists(glygen_path)) :
        os.mkdir(glygen_path)

glygen_now = os.path.join(glygen_path,now)
if not (os.path.exists(glygen_now)) :
        os.mkdir(glygen_now)

path_intermediate = os.path.join(working_path, 'intermediate')
if not (os.path.exists(path_intermediate)) :
	os.mkdir(path_intermediate)

path_generated = os.path.join(working_path, 'generated')
if not (os.path.exists(path_generated)) :
	os.mkdir(path_generated)
