#This program downloads data from UniProtKB for human and mouse 1:1 glycosyltransferase orthologs with differential expression in cancer

import sys
import urllib
import os
import cons_config
import glob
import csv

#Assign variable for datetime
from datetime import datetime
today = datetime.now()
now = today.strftime('%Y_%m_%d')

###############################################
###			  		    ###
### Set working directory to downloads path ###
###			  		    ###
###############################################


oth_downloads = cons_config.path_downloads

oma_now = cons_config.oma_now

glygen_now = cons_config.glygen_now

uniprot_now = cons_config.uniprot_now

bgee_now = cons_config.bgee_now

oncomx_now = cons_config.oncomx_now

path_intermediate = cons_config.path_intermediate

path_generated = cons_config.path_generated

path_diff_gen = '/home/hmhamilt/phd_dissertation/diff_exp/generated'

os.chdir(path_intermediate)

##########################
###		       ###
### Print current task ###
###		       ###
##########################

print('Downloading UniProtKB gene info for human and mouse and creating ortholog table')

##############################
###			   ###
### Retrieve file from URL ###
###			   ###
##############################

#Define mouse gts
print('Defining lists of human and mouse GTs')

mm_gt_details = glygen_now + '/GLY_000030.csv'

output_rows = []
with open(mm_gt_details, "r") as mm_gt_list:
        mm_gt_csv = csv.reader(mm_gt_list)
        next(mm_gt_csv)
        mm_glygen_gts = []
        for row in mm_gt_csv :
                if row[0] in mm_glygen_gts :
                        continue
                else:
			row[0] = row[0].split('-', 1)[0]
                        mm_glygen_gts.append(row[0])
        print(len(mm_glygen_gts))

for glygen_gt in mm_glygen_gts :
        output_row = []
        output_row.append(glygen_gt)
        output_rows.append(output_row)

with open('mouse_glygen_gts.tsv', 'wb') as csvfile :
        writer = csv.writer(csvfile, delimiter = '\t')
        writer.writerows(output_rows)

#Make master ortholog table from OMA and GlyGen data
print('Printing ortholog table')

output_rows = []
header = ['human_uniprotkb_ac', 'human_gene', 'mouse_uniprotkb_ac', 'mouse_gene', 'is_glygen']
output_rows.append(header)

with open('gt_human_mouse_orthologs.tsv', 'r') as orth_file:
	orthologs = csv.reader(orth_file, delimiter="\t")
	next(orthologs)
	for ortholog in orthologs:
		uniprot_ac = ortholog[0]
		mouse_ac = ortholog[1]
		print(uniprot_ac)
		#uniprot_ac = uniprot_ac.rstrip('\n')
		os.chdir(uniprot_now)
		url = "https://www.uniprot.org/uniprot/?query=organism:9606+AND+reviewed:yes+AND+accession:" + uniprot_ac + "&format=tab&columns=id,genes(PREFERRED)"
		urllib.urlretrieve (url, filename="uniprot_" + uniprot_ac + ".tab")
		url_mouse = "https://www.uniprot.org/uniprot/?query=organism:10090+AND+reviewed:yes+AND+accession:" + mouse_ac + "&format=tab&columns=id,genes(PREFERRED)"
                urllib.urlretrieve (url_mouse, filename="uniprot_" + mouse_ac + ".tab")

		output_row = []
		with open('uniprot_' + uniprot_ac + '.tab', "r") as human_file:
			human = csv.reader(human_file, delimiter="\t")
			i = 0
			#next(human) - next is breaking the loop in this case, so using index instead
			for hum_gt in human:
				#hum_gt = hum_gt.rstrip('\n')
				hum_ac = hum_gt[0]
				hum_gene = hum_gt[1]
				i+=1
				if i==1:
					continue
				elif hum_ac == "":
					continue
				else:
					if hum_ac == uniprot_ac:
						print(hum_gt)
						#QC check to make sure rigth file is being opened"
						output_row.append(hum_ac)
						output_row.append(hum_gene)
					else:
						continue
		
		with open("uniprot_" + mouse_ac + ".tab", 'r') as mouse_file:
                        mouse = csv.reader(mouse_file, delimiter="\t")
                        i = 0
                        #next(mouse)
			for mm_gt in mouse:
                                #hum_gt = hum_gt.rstrip('\n')
                                mm_ac = mm_gt[0]
                                mm_gene = mm_gt[1]
                                i+=1
                                if i==1:
                                       continue
				elif mm_ac == "":
					continue
                                else:
                                	if mm_ac == mouse_ac:
						print(mm_gt)
                                        	#QC check to make sure rigth file is being opened and right mouse protein appended"
	                                        output_row.append(mm_ac)
        	                                output_row.append(mm_gene)
                	                else:
                        	                continue

		os.chdir(path_intermediate)
		
		with open('mouse_glygen_gts.tsv', "r") as glygen_file:
			glygen_lines = csv.reader(glygen_file, delimiter="\t")		
			for line in glygen_lines:
				glygen_ac = line[0]	
				if glygen_ac == mouse_ac:
					output_row.append("yes")
				else: 
					continue 
		output_rows.append(output_row)

with open('orthologs_with_info.tsv', 'wb') as newfile:
	writer = csv.writer(newfile, delimiter = '\t')
        writer.writerows(output_rows)

#Filter ortholog file for only those with mouse ortholog in OMA
cmd = """awk -F'\t' '$3 != \"\"' orthologs_with_info.tsv > """ + path_generated + """/orthologs_info_filtered.tsv"""
print(cmd)
os.system(cmd)

#Print completion of current task
print('Done')


