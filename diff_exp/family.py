#This script summarizes family counts for GTs in cancer.

import sys
import urllib
import os
import diff_exp_config
import csv
import re

#Assign variable for datetime
from datetime import datetime
today = datetime.now()
now = today.strftime('%Y_%m_%d')

###############################################
###                                         ###
### Set working directory to downloads path ###
###                                         ###
###############################################

path_downloads = diff_exp_config.path_downloads

path_intermediate = diff_exp_config.path_intermediate

path_generated = diff_exp_config.path_generated

bioxpress_now = diff_exp_config.bioxpress_now

glygen_now = diff_exp_config.glygen_now

os.chdir(path_generated)

##########################
###                    ###
### Print current task ###
###                    ###
##########################

print('Defining variables')

##################################################
###			       			###
### Extracting cancers and families from tables ###
###		 			   	###
###################################################


#EXAMPLE ROW:
#0UniProtKB_AC   1 RefSeq  2Gene    3log2FoldChange  4p_value 5adjusted p_value        6Significant    7 Trend   8TCGA Cancer    9 Cancer Ontology 10PMID   11 UBERON_ID      12 PtsTrend       13 PtsUp  14 PtsDown15 TotPts 16 Percent
#A0PJZ3  NP_001073862.1  GXYLT2  -1.0    1.23e-02        3.19e-02        Yes     Down    BLCA    DOID:11054 / Urinary bladder cancer [UBC]       -       UBERON:0001255 12       7       12      19      63.16
#Define inputs

file = 'de_gts_master_filtered.tsv'
#fc_filt = 'de_gts_master_fc_filtered.tsv'
#fc_high = 'de_gts_master_fc_filtered_high.tsv'
#fc_higher = 'de_gts_master_fc_filtered_higher.tsv'
#fc_highest = 'de_gts_master_fc_filtered_highest.tsv'

#filtered = [filt, fc_filt, fc_high, fc_higher, fc_highest]

#Define variables

os.chdir(path_intermediate)
	
cancers_clean = []
with open('cancers_no_manual.tsv','r') as csv_file:
    cancers = csv_file.readlines()
    print(cancers)
    for cancer in cancers:
	cancers_clean.append(cancer.strip())
    print(cancers_clean)

#gts_clean = []
#with open(path_generated + '/temp.tsv','r') as csv_file:
#    gts = csv_file.readlines()
#    for gt in gts:
#	gts_clean.append(gt.strip())

cancer_total = len(cancers_clean)

print('Cancers ' + str(cancer_total))
#print('GTs ' + str(len(gts_clean)))

os.system('cp /home/hmhamilt/phd_dissertation/gts_list/intermediate/gts_append_go.tsv ' + path_intermediate + '/gts_append_go.tsv')
cmd = """awk -F'\t' 'NR!=1{print $5}' gts_append_go.tsv | sort -u > families.tsv"""
os.system(cmd)
os.system('head families.tsv')

families = []
families_clean = []
families_sep = []
with open('families.tsv', 'r') as fam_file:
	families = fam_file.readlines()
	for family in families:
		families_clean.append(family.strip(";\n"))
	for family in families_clean:
		if ";" in family:
			family = family.split(";")
			i = 0
			#family = "\n".join(family)
			while i < len(family):
				families_sep.append(family[i])
				i += 1
		elif family == "":
			continue
		else: 
			families_sep.append(family)
families_sort = sorted(set(families_sep))
families_final = []
for family in families_sort:
	if "GT" in family:
		families_final.append(family)
	else:
		continue
print(families_final)

cancer_total = len(cancers_clean)
families_total = len(families_final)			
print('Cancers: ' + str(cancer_total))
print('Families: ' + str(families_total))

output_rows = []
header = ['cancer', 'family', 'total_in_family', 'total_uniprot', 'freq', 'total_cancer_degts', 'expected_family_cancer_degts', 'actual_family_cancer_degts', 'difference']

output_rows.append(header)

#Counting family GTs across all cancers

for family in families_final:
	output_row = []
	total_family_gts = 0
	family_gts = []
	with open(glygen_now + '/GLY_000004.csv', 'r') as csv_file:
		#lines = csv_file.readlines()
		lines = csv.reader(csv_file, delimiter=",")
		for line in lines:
			if "|" in line[7]:
               	        	test = line[7].split("|")
                       		i = 0
                       		#family = "\n".join(family)
	                        while i < len(test):
       		                        if family == test[i]:
						total_family_gts +=1
						family_gts.append(line[0])
               		                i += 1
			elif family == line[7]:
				total_family_gts +=1
				family_gts.append(line[0])
			else:
				continue
	
	print(family + ': ' + str(len(sorted(set(family_gts)))))
	fam_gts_clean = []
	for gt in family_gts:
		gt.strip("'")
		gt = gt.split("-")
		gt = gt[0]
		fam_gts_clean.append(gt)
	fam_gts_clean = sorted(set(fam_gts_clean))
	print(family + ': ' + str(len(fam_gts_clean)))
	total_uniprot = 20359
	freq = float(total_family_gts) / total_uniprot
	total_all_cancer_sig_degs = 0
	all_sig_degs = []
	with open('BioXpress_all_sig.tsv', 'r') as biox_file:
		lines = csv.reader(biox_file, delimiter="\t")
		next(lines)
		for line in lines:
			all_sig_degs.append(line[0])

	all_sig_degs = sorted(set(all_sig_degs))
	total_all_cancer_sig_degs = len(all_sig_degs)
	print("All significant DEGs: " + str(total_all_cancer_sig_degs))

	expected_family_all_cancer_sig_degts = total_all_cancer_sig_degs*freq

	actual_family_all_cancer_sig_degts = 0

        for deg in all_sig_degs:
		if deg in fam_gts_clean:
                	actual_family_all_cancer_sig_degts += 1

        difference = actual_family_all_cancer_sig_degts - expected_family_all_cancer_sig_degts

        output_row.append("All")
        output_row.append(family)
        output_row.append(total_family_gts)
        output_row.append(total_uniprot)
        output_row.append(freq)
        output_row.append(total_all_cancer_sig_degs)
        output_row.append(expected_family_all_cancer_sig_degts)
        output_row.append(actual_family_all_cancer_sig_degts)
        output_row.append(difference)
        output_rows.append(output_row)

#Counting family GTs for individual cancers

for cancer in cancers_clean:
	for family in families_final:
		output_row = []
		total_family_gts = 0
	        family_gts = []
        	with open(glygen_now + '/GLY_000004.csv', 'r') as csv_file:
                	lines = csv.reader(csv_file, delimiter=",")
	                for line in lines:
        	                if "|" in line[7]:
                	                test = line[7].split("|")
                        	        i = 0
                                	while i < len(test):
                                        	if family == test[i]:
                                                	total_family_gts +=1
	                                                family_gts.append(line[0])
        	                                i += 1
                	        elif family == line[7]:
                        	        total_family_gts +=1
                                	family_gts.append(line[0])
	                        else:
        	                        continue

	        print(family + ': ' + str(len(sorted(set(family_gts)))))
        	fam_gts_clean = []
	        for gt in family_gts:
        	        gt.strip("'")
                	gt = gt.split("-")
	                gt = gt[0]
        	        fam_gts_clean.append(gt)
	        fam_gts_clean = sorted(set(fam_gts_clean))
        	print(family + ': ' + str(len(fam_gts_clean)))
	        total_uniprot = 20359
        	freq = float(total_family_gts) / total_uniprot
		total_cancer_sig_degs = 0
	        cancer_sig_degs = []
        	with open('BioXpress_all_sig.tsv', 'r') as biox_file:
                	lines = csv.reader(biox_file, delimiter="\t")
	                next(lines)
        	        for line in lines:
				if line[8] == cancer:
	                	        cancer_sig_degs.append(line[0])
					#if line[0] in fam_gts_clean:
                                        #        actual_family_cancer_degts += 1

	        cancer_sig_degs = sorted(set(cancer_sig_degs))
	        total_cancer_sig_degs = len(cancer_sig_degs)
	        print(cancer + " significant DEGs: " + str(total_cancer_sig_degs))

	        expected_family_cancer_sig_degts = total_cancer_sig_degs*freq
		
		actual_family_cancer_sig_degts = 0		

	        for deg in cancer_sig_degs:
        	        if deg in fam_gts_clean:
                	        actual_family_cancer_sig_degts += 1

	        difference = actual_family_cancer_sig_degts - expected_family_cancer_sig_degts

		#total_cancer_sig_degs = 0
		#actual_family_cancer_degts = 0
		#with open(path_generated + '/de_gts_master_reordered.tsv', 'r') as sig_file:
		#	lines = csv.reader(sig_file, delimiter="\t")
		#	for line in lines:
		#		if line[8] == cancer:
		#			total_cancer_degts += 1
		#				if line[0] in fam_gts_clean:
		#					actual_family_cancer_degts += 1
		output_row.append(cancer)
		output_row.append(family)
		output_row.append(total_family_gts)
		output_row.append(total_uniprot)
		output_row.append(freq)
		output_row.append(total_cancer_sig_degs)
		output_row.append(expected_family_cancer_sig_degts)
		output_row.append(actual_family_cancer_sig_degts)
		output_row.append(difference)					
		output_rows.append(output_row)


os.chdir(path_generated)

with open("enrichment_by_family.tsv", "w") as out_file:
	writer = csv.writer(out_file, delimiter='\t')
	writer.writerows(output_rows)

