#This script parses the BioXpress v2 table..

import sys
import urllib
import os
import diff_exp_config
import bs4 as bs
import csv
import pandas as pd

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

bioxpress_now = diff_exp_config.bioxpress_now

glygen_now = diff_exp_config.glygen_now

os.chdir(path_downloads)

##########################
###                    ###
### Print current task ###
###                    ###
##########################

print('Parsing tables')

##################################
###			       ###
### Extracting gts from tables ###
###			       ###
##################################

#Define gts
print('Defining list of GTs')

gt_details = glygen_now + '/GLY_000004.csv'
output_rows = []

with open(gt_details, "r") as gt_list:
	gt_csv = csv.reader(gt_list)
	next(gt_csv)
	gts = []
	for row in gt_csv :
		if row[0] in gts :
			continue
		else:
			row[0] = row[0].split('-', 1)[0]
	 		gts.append(row[0])
	print(len(gts))

for gt in gts :
	output_row = []
	output_row.append(gt)
	output_rows.append(output_row)

os.chdir(path_intermediate)

with open('gts.tsv', 'wb') as csvfile :
        writer = csv.writer(csvfile, delimiter = '\t')
        writer.writerows(output_rows)

#Parse file

v2 = bioxpress_now + '/BioXpress_interface_overall_final_v2.0.csv'
#v4_study = bioxpress_now + '/human_cancer_mRNA_expression_per_study.csv'
#v4_tissue = bioxpress_now + '/human_cancer_mRNA_expression_per_tissue.csv'

#tables = [v2,v4_study,v4_tissue]

table = v2
output_rows = []

#for table in tables :
print(table)
#if table == v2 :
with open(table, "r") as input:
	input_csv = csv.reader(input)
	next(input_csv)
	cancers = []
	for row in input_csv :
		#print(row)
		#print(row[8])
		if row[8] in cancers:
			continue
		else:
			cancers.append(row[8])
	#print(len(cancers))
	#Shouldn't be necessary since we're saying to skip it if it's already there
	cancers = sorted(set(cancers))
	print(cancers)

for cancer in cancers :
        output_row = []
        output_row.append(cancer)
        output_rows.append(output_row)

with open('cancers.tsv', 'wb') as csvfile :
        writer = csv.writer(csvfile, delimiter = '\t')
        writer.writerows(output_rows)

	#Make copy of table so we don't accidentally overwrite
	#os.system('cp ' + table + ' ' + path_intermediate + '/v2_copy.csv')
	#os.system('ls ' + path_intermediate)
	#os.chdir(path_intermediate)
	#os.system('pwd')

			#table_copy = 'v2_copy.csv'

			#Example header and row of table_copy
			#UniProtKB_AC,RefSeq,Gene,log2FoldChange,p_value,adjusted p_value,Significant,Trend,TCGA Cancer,Cancer Ontology,#Patients,Data Source,PMID,UBERON_ID
			#Q6UXB8,NP_001186088.1; NP_699201.2; XP_005248974.1; XP_011512677.1,PI16,-7.72,4.13e-46,8.02e-42,Yes,Down,BLCA,DOID:11054 / Urinary bladder cancer [UBC],19/19(100.0),RNASeqV2,-,UBERON:0001255

			#Split patient counts into separate columns
			#df = pd.read_csv(table_copy)
			#df.dropna(inplace = True)
			#new = df["#Patients"].str.split("/", n = 1, expand = True)
			#df["No. with Trend"]= new[0]
			#df[#Placeholder"]= new[1]
			#df.drop(columns = ["#Patients"], inplace = True)
			#with open("temp.csv", "w") as file:
			#	df.to_csv(file, index=False) 
			#print(data)

	#else : 
	#	continue


#for gt in total_gts :
#        output_row = []
#        output_row.append(gt)
#        output_rows.append(output_row)



#with open('cfg.tsv', 'wb') as csvfile :
#    writer = csv.writer(csvfile, delimiter = '\t')
#    writer.writerows(output_rows)

#Print completion of current task
print('Done')

