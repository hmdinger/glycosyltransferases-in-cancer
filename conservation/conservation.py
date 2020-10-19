#This program downloads data from UniProtKB for human and mouse 1:1 glycosyltransferase orthologs with differential expression in cancer

import sys
import urllib
import os
import cons_config
import glob
import csv
import pandas as pd
import numpy as np

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

uniprot_now = cons_config.uniprot_now

bgee_now = cons_config.bgee_now

oncomx_now = cons_config.oncomx_now

path_intermediate = cons_config.path_intermediate

path_generated = cons_config.path_generated

path_diff_gen = '/home/hmhamilt/phd_dissertation/diff_exp/generated'

path_diff_int = '/home/hmhamilt/phd_dissertation/diff_exp/intermediate'

os.chdir(path_generated)

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

#Define gts from orthologs table
print('Defining GTs')

output_rows = []
with open('orthologs_info_filtered.tsv', "r") as gt_file:
        gt_rows = csv.reader(gt_file, delimiter="\t")
        next(gt_rows)
        ortho_gts = []
	gt_mouse_genes = []
        for row in gt_rows:
                if row[0] in ortho_gts :
                        continue
                else:
                        ortho_gts.append(row[0])
			gt_mouse_genes.append(row[3])

print('There are ' + str(len(ortho_gts)) + ' gts with ortholog and expression information.')

for gt in ortho_gts :
        output_row = []
        output_row.append(gt)
        output_rows.append(output_row)

with open('ortholog_gts.tsv', 'wb') as csvfile :
        writer = csv.writer(csvfile, delimiter = '\t')
        writer.writerows(output_rows)

#Make human normal expression file just for the gts
print('Filtering human normal expression data for gts with orthologs')
head = """head -n 1 """ + oncomx_now + """/human_normal_expression_custom.csv > """ + path_intermediate + """/human_normal_custom_filtered.csv"""
os.system(head)
for gt in ortho_gts :
	cmd = """awk -F',' '$2 == \"\\\"""" + gt + """\\\"\"' """ + oncomx_now + """/human_normal_expression_custom.csv >> """ + path_intermediate + """/human_normal_custom_filtered.csv"""
	os.system(cmd)

#Make mouse normal expression file just for the gts
print('Filtering mouse normal expression data for gts with orthologs')
head = """head -n 1 """ + oncomx_now + """/mouse_normal_expression_custom.csv > """ + path_intermediate + """/mouse_normal_custom_filtered.csv"""
os.system(head)
for gene in gt_mouse_genes :
        cmd = """awk -F',' '$1 == \"\\\"""" + gene + """\\\"\"' """ + oncomx_now + """/mouse_normal_expression_custom.csv >> """ + path_intermediate + """/mouse_normal_custom_filtered.csv"""
        os.system(cmd)

#Define cancers
#cancers_clean = []
#with open(path_diff_int + '/cancers_no_manual.tsv','r') as csv_file:
#            cancers = csv_file.readlines()
#            print(cancers)
#            for cancer in cancers:
#                cancers_clean.append(cancer.strip())
#            print(cancers_clean)

#Define cancer-tissue map from unique tissues in bgee file
print('Making a dictionary for cancer terms and their corresponding tissue')
map = {
	"BLCA":	"urinary bladder",
	"BRCA":	"breast",
	"CESC":	"uterine cervix",
	"ESCA":	"esophagus",
	"HNSC": "no match",
	"KICH":	"kidney",
	"KIRC":	"kidney",
	"KIRP":	"kidney",
	"LIHC": "liver",
	"LUAD":	"lung",
	"LUSC":	"lung",
	"PAAD":	"pancreas",
	"PRAD":	"prostate gland",
	"STAD":	"stomach",
	"THCA":	"thyroid gland",
	"UCEC":	"uterus"
}
print(map)
print(map["BLCA"])
#print(map.values())

#Print tissues
#tissues_list = []
output_rows = []
#for tissue in set(map.values()):
#	tissues = []
#	if tissue in tissues:
#		continue
#	else: 
#		tissues.append(tissue)
#	tissues_list.append(tissues)
#	output_rows.append(tissues_list)
	
#with open('tissue_list.csv', "w") as tissue_file:
#	writer = csv.writer(tissue_file)
#	writer.writerow(tissues_list)

for tissue in set(map.values()) :
	print(tissue)
        output_row = []
        output_row.append(tissue)
       	output_rows.append(output_row)

with open('ortholog_tissues.tsv', 'wb') as csvfile :
        writer = csv.writer(csvfile, delimiter = '\t')
        writer.writerows(output_rows)



#Make master expression conservation table
print('Generating expression conservation table')

orth_file = pd.read_csv('orthologs_info_filtered.tsv', sep="\t")
de_file = pd.read_csv(path_diff_gen + '/de_gts_master_filtered.tsv', sep="\t")
hs_norm = pd.read_csv(path_intermediate + '/human_normal_custom_filtered.csv', sep=",")
mm_norm = pd.read_csv(path_intermediate + '/mouse_normal_custom_filtered.csv', sep=",")

print(orth_file.head())
print(de_file.head())
print(hs_norm.head())
print(mm_norm.head())

common1 = pd.merge(orth_file, de_file, how='left', left_on=['human_uniprotkb_ac'], right_on=['UniProtKB_AC'])
print(common1.head())
common1['tissue'] = common1['TCGA Cancer'].map(map)
print(common1.head())

common2 = pd.merge(common1, hs_norm, how='left', left_on=['human_uniprotkb_ac', 'tissue'], right_on=['uniprotkb_ac', 'uberon_anatomical_name_short'])
common2_fixed = common2.dropna(subset=['uberon_anatomical_name_short'])
print(common2_fixed.head())

common3 = pd.merge(common2_fixed, mm_norm, how='left', left_on=['mouse_gene', 'tissue'], right_on=['gene_symbol', 'uberon_anatomical_name_short'])
common3_fixed = common3.dropna(subset=['uberon_anatomical_name_short_y'])
print(common3_fixed.head())

for col in common3_fixed.columns:
	print(col)

common3_better_headers = common3_fixed.rename(columns = {'TCGA Cancer': 'cancer', 'Trend': 'cancer_trend', 'adjusted p_value': 'adjusted_p_value', 'log2FoldChange': 'log2fc', 'uberon_anatomical_id_x': 'uberon_id', 'uberon_anatomical_name_short_x': 'uberon_name', 'expression_level_anatomical_relative_x': 'human_healthy_expression', 'expression_score_x': 'human_score', 'expression_level_anatomical_relative_y': 'mouse_healthy_expression', 'expression_score_y': 'mouse_score'})

score = {
        "HIGH": 1,
        "MEDIUM": 0,
        "LOW": -1,
        "ABSENT": -10,
	"Up": -2,
	"Down": 2
}

print(score)
print(score["ABSENT"])

common3_better_headers['conservation_score'] = (common3_better_headers['human_healthy_expression'].map(score)) + (common3_better_headers['mouse_healthy_expression'].map(score))

common3_better_headers['impact_score'] = (common3_better_headers['conservation_score']) + (common3_better_headers['cancer_trend'].map(score))

print(common3_better_headers.head())

#common3.to_csv('pandas_test3_fixed.tsv', sep='\t', columns=['human_uniprotkb_ac', 'human_gene', 'mouse_uniprotkb_ac', 'mouse_gene', 'TCGA Cancer', 'Trend', 'adjusted p_value', 'log2FoldChange', 'uberon_anatomical_id_x', 'uberon_anatomical_name_short_x', 'expression_level_anatomical_relative_x', 'expression_score_x', 'expression_level_anatomical_relative_y', 'expression_score_y'], index=False)

common3_better_headers.to_csv('conservation_master.tsv', sep='\t', columns=['human_uniprotkb_ac', 'human_gene', 'mouse_uniprotkb_ac', 'mouse_gene', 'cancer', 'cancer_trend', 'adjusted_p_value', 'log2fc', 'uberon_id', 'uberon_name', 'human_healthy_expression', 'human_score', 'mouse_healthy_expression', 'mouse_score', 'conservation_score', 'impact_score'], index=False)

print('Generating table filtered for absent calls')
output_rows = []
with open('conservation_master.tsv', "r") as master_file:
        rows = csv.reader(master_file, delimiter="\t")
        i = 0
        for row in rows:
                if i == 0:
                        output_rows.append(row)
                else:
                        if int(row[15]) < -4:
                                continue
                        else:
                                output_rows.append(row)
                i += 1

with open('no_absent_conservation.tsv', "w") as out_file:
        writer = csv.writer(out_file, delimiter = "\t")
        writer.writerows(output_rows)

print('Generating subset of high value GTs based on conservation of expression profile')
output_rows = []
with open('conservation_master.tsv', "r") as master_file:
	rows = csv.reader(master_file, delimiter="\t")
	i = 0
	for row in rows:
		if i == 0:
			output_rows.append(row)
		else:
			if (row[15] == '4') or (row[15] == '-4'):
				output_rows.append(row)
			else:
				continue
		i += 1

with open('high_value_conservation.tsv', "w") as out_file:
	writer = csv.writer(out_file, delimiter = "\t") 
	writer.writerows(output_rows)

#Filter tables for higher log2FC
with open("de_gts_master_filtered.tsv", "w") as out_file:
        writer = csv.writer(out_file, delimiter='\t')
        writer.writerows(output_rows)

output_rows=[]
print('Filtering for increased foldchange')
with open('high_value_conservation.tsv', "r") as cons_file:
        cons_file_csv = csv.reader(cons_file, delimiter="\t")
        needheader=True
        for row in cons_file_csv:
                if needheader:
                        needheader=False
                        output_rows.append(row)
                        continue
                elif (((float(row[7])) >= 1.0) == True):
                        output_rows.append(row)
                elif (((float(row[7])) <= -1.0) == True):
                        output_rows.append(row)
                else:
                        continue

with open("high_value_conservation_high_threshold.tsv", "w") as out_file:
        writer = csv.writer(out_file, delimiter='\t')
        writer.writerows(output_rows)

print('Done')


