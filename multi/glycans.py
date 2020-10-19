#This program generates a table for glycans impacted by DEGTs in cancer and counts affected glycans/linkages per cancer

import sys
import urllib
import os
import multi_config
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

path_downloads = multi_config.path_downloads

pmc_now = multi_config.pmc_now

glygen_now = multi_config.glygen_now

oncomx_now = multi_config.oncomx_now

hgnc_now = multi_config.hgnc_now

manual_path = multi_config.manual_path

path_intermediate = multi_config.path_intermediate

path_diff_gen = '/home/hmhamilt/phd_dissertation/diff_exp/generated'

path_diff_down_glygen = '/home/hmhamilt/phd_dissertation/diff_exp/downloads/glygen/2020_10_12'

path_generated = multi_config.path_generated

os.chdir(path_intermediate)

##########################
###		       ###
### Print current task ###
###		       ###
##########################

print('Downloading UniProtKB gene info for human and mouse and creating ortholog table')

##############################
###			   ###
### Define variables       ###
###			   ###
##############################

#Define gts from gts_list
print('Defining GTs')

copy = 'cp ' + path_diff_down_glygen + '/GLY_000004.csv glygen_gts.csv'
os.system(copy)

output_rows = []
with open('glygen_gts.csv', "r") as gt_file:
        gt_rows = csv.reader(gt_file)
        next(gt_rows)
        gts = []
        for row in gt_rows:
                if row[0] in gts :
                        continue
                else:
			stripped = row[0].replace('"','')
			split = stripped.split('-', 1)[0]
                        gts.append(split)

print('There are ' + str(len(gts)) + 'gts.')

print(gts)

#Making dictionary for gts with genes
print('Making a dictionary for gt UniProtKB Acs and corresponding gene symbols')
check = []
gt_dict = {}
with open('glygen_gts.csv', "r") as gt_file:
        gt_rows = csv.DictReader(gt_file)
        next(gt_rows)
        for d in gt_rows:
                if d in check:
                        continue
                else:
                        check.append(d)
                        gt_dict.setdefault(d['uniprotkb_canonical_ac'],d['gene_symbol'])

#print(gt_dict)
gts = gt_dict.keys()
genes = gt_dict.values()

#Remove quotes and hyphens"
for gt in gts:
	gt_stripped = gt.replace('"','')
	gt_split = gt.split('-', 1)[0]
	gt_dict[gt_split] = gt_dict.pop(gt)

#print(gt_dict)
gts = gt_dict.keys()
genes = gt_dict.values()


###############################################
### Make master table for glycans and DEGTs ###
###############################################

#Filter glycan table for human gt glycans
print('Filtering input for human glycans')
###NOTE - using file from manual path provided directly by Rahi
#copy = 'cp ' + glygen_now + '/GLY_000284.csv glycans.csv'
copy = 'cp ' + manual_path + '/annotated_glycans.csv glycans.csv'
os.system(copy)

output_rows = []
for gt in gts:
	with open('glycans.csv', "r") as gly_file:
		rows = csv.reader(gly_file)
		for row in rows:
			if output_rows == []:
				output_rows.append(row)
			else:
				enzyme = row[3]
                        #	enzyme = enzyme.replace('"','')
                        #	enzyme = enzyme.split('-', 1)[0]
				if enzyme == gt:
					output_rows.append(row)
				else:
					continue
	
with open('glycans_human.csv', "w") as out_file:
	writer = csv.writer(out_file)
	writer.writerows(output_rows)

#Merge with cancer table
print('Joining glycans with DEGTs')

de_file = pd.read_csv(path_diff_gen + '/de_gts_master_filtered.tsv', sep="\t")
gly_file = pd.read_csv(path_intermediate + '/glycans_human.csv', sep=",")

print(de_file.head())
print('Differential expression input file has ' + str(len(de_file.index)) + ' rows.')
print(gly_file.head())
print('Glycan input file has ' + str(len(gly_file.index)) + ' rows.')

#NOTE - following steps needed for reformatting Acs of GLY_000284 but not needed for annotated_glycans.csv. 
#Additional commented out steps in this section are pertinent to GLY_000284.

#gly_file['uniprotkb_canonical_ac'] = gly_file['uniprotkb_canonical_ac'].str.split('-', 1).str[0]
#print(gly_file['uniprotkb_canonical_ac'].head())
#gly_file.to_csv('test.tsv', sep='\t')

#common = pd.merge(de_file, gly_file, how='left', left_on=['UniProtKB_AC'], right_on=['uniprotkb_canonical_ac'])
common = pd.merge(de_file, gly_file, how='left', left_on=['UniProtKB_AC'], right_on=['uniprot'])
print(common.head())
print('Merged DEGT-glycan file has ' + str(len(common.index)) + ' rows.')

#common_full = common.dropna(subset=['uniprotkb_canonical_ac'])
common_full = common.dropna(subset=['uniprot'])
print(common_full.head())
print('After removal of duplicates, merged file has ' + str(len(common_full.index)) + ' rows.')

print('Renaming headers')
better_headers = common_full.rename(columns = {'UniProtKB_AC': 'uniprotkb_ac', 'Gene': 'gene', 'TCGA Cancer': 'cancer', 'Trend': 'cancer_trend', 'adjusted p_value': 'adjusted_p_value', 'log2FoldChange': 'log2fc'})
print(better_headers.head())

print('Writing merged file to tsv')
better_headers.to_csv(path_generated + '/gt_glycans_with_glyc.tsv', sep='\t', columns=['uniprotkb_ac', 'gene', 'cancer', 'cancer_trend', 'adjusted_p_value', 'log2fc', 'glytoucan_ac', 'residue_name', 'residue_id','parent_residue_id'], index=False)

print('Removing glycan column from merged file')
better_headers.to_csv(path_generated + '/gt_glycans_no_glyc.tsv', sep='\t', columns=['uniprotkb_ac', 'gene', 'cancer', 'cancer_trend', 'adjusted_p_value', 'log2fc', 'residue_name', 'residue_id', 'parent_residue_id'], index=False)

##########################################
### Make high value table for new file ###
##########################################
print('Generating high value tables based on glycans affected by DE and printing dataframes for images')

thresholds = [1.0,2.0,3.0]
for threshold in thresholds:
	in_files = ['gt_glycans_with_glyc.tsv', 'gt_glycans_no_glyc.tsv']
	for file in in_files:
		file_split = file.split(".",1)[0]
		print(file, threshold, file_split)
		output_rows = []
		with open(path_generated + '/' + file, "r") as gtgly_file:
			rows = csv.reader(gtgly_file, delimiter="\t")
			needheader=True
			for row in rows:
				if needheader:
					needheader=False
					output_rows.append(row)
				elif (((float(row[5])) >= threshold) == True) or (((float(row[5])) <= -threshold) == True) :
					output_rows.append(row)
				else:
					continue

		with open(path_generated + '/' + file_split + '_' + str(threshold) + '.tsv', "w") as outfile:
			writer = csv.writer(outfile, delimiter = "\t")
			writer.writerows(output_rows)

		df_file = pd.read_csv(path_generated + '/' + file_split + '_' + str(threshold) + '.tsv', sep="\t")
		print('There are ' + str(len(df_file.index)) + ' rows with |log2fc| >= ' + str(threshold))
		if file == 'gt_glycans_no_glyc.tsv':
			print(file)
			df_short = df_file[['cancer', 'cancer_trend', 'gene', 'residue_id']]
		else:
			print(file)
			df_short = df_file[['cancer', 'cancer_trend', 'gene', 'glytoucan_ac']]
		print(df_short.head())
		print('There are ' + str(len(df_short.index)) + ' rows in the resulting data frame file.')
		df_short_unique = df_short.drop_duplicates()
		print(df_short_unique.head())
		print('There are ' + str(len(df_short_unique.index)) + ' rows remaining after duplicate filtration.')
                #df_short_unique.to_csv(path_generated + '/' + file + 'for_sankey' + str(threshold) + '.tsv', sep='\t', columns=['cancer', 'cancer_trend', 'gene', 'glytoucan_ac'], index=False)
		df_short_unique.to_csv(path_generated + '/' + file_split + '_for_sankey_' + str(threshold) + '.tsv', sep='\t', index=False)

		
###################################
### Summarize counts by gene    ###
###################################

#Summarize glycan counts by gene for |log2fc| >= 1.0
print('Summarize counts of glycans by gene')
output_rows = []
header = ['uniprotkb_ac','gene','disease_count','glycan_count','position_count']
output_rows.append(header)

summary = pd.read_csv(path_generated + '/gt_glycans_filt_1.0.tsv', sep="\t")
print(summary.head())

for gt in gts:
        gene = gt_dict[gt]
        if gt in summary.values:
		output_row = []
                summary_gt = summary.loc[summary['uniprotkb_ac'] == gt]
                disease_list = summary_gt['cancer'].unique()
                disease_count = len(disease_list)
		glycan_list = summary_gt['glytoucan_ac'].unique()
		glycan_count = len(glycan_list)
		position_list = summary_gt['residue_id'].unique()
		position_count = len(position_list)
		output_row.append(gt)
		output_row.append(gene)
		output_row.append(disease_count)
		output_row.append(glycan_count)
		output_row.append(position_count)
		output_rows.append(output_row)
	else:
		continue

with open(path_generated + "/glycan_counts_by_gene.tsv", "w") as out_file:
	writer = csv.writer(out_file, delimiter = "\t")
	writer.writerows(output_rows)

###################################
### Summarize by glycan         ###
###################################

output_rows = []
header = ['glytoucan_ac','total_gts','gts_up_count','gts_down_count','position_count']
output_rows.append(header)

summary = pd.read_csv(path_generated + '/gt_glycans_filt_1.0.tsv', sep="\t")
print(summary.head())

glycan_list = summary['glytoucan_ac'].unique()

for glycan in glycan_list:
	output_row = []
        summary_glycan = summary.loc[summary['glytoucan_ac'] == glycan]
        uniprot_list = summary_glycan['uniprotkb_ac'].unique()
        genes_count = len(uniprot_list)
        genes_down = summary_glycan.loc[summary['cancer_trend'] == "DOWN"]
        genes_down_count = len(genes_down)
        genes_up = summary_glycan.loc[summary['cancer_trend'] == "UP"]
        genes_up_count = len(genes_up)
        gts_count = 0
        gts_up_count = 0
        gts_down_count = 0
        #for gt in gts:
         #       if gt in uniprot_list:
          #              gts_count +=1
           #     if gt in genes_down:
            #            gts_down_count +=1
             #   if gt in genes_up:
              #          gts_up_count +=1
        position = summary_glycan['residue_id'].unique()
        position_count = len(position_list)
	output_row.append(glycan)
        output_row.append(genes_count)
	output_row.append(genes_up_count)
	output_row.append(genes_down_count)
        output_row.append(position_count)

with open(path_generated + "/glycan_counts_by_glycan.tsv", "w") as out_file:
        writer = csv.writer(out_file, delimiter = "\t")
        writer.writerows(output_rows)

##################################
### Summarize counts by cancer ###
##################################

output_rows = []
#header = ['cancer','total_genes','down_genes','up_genes','glycan','position']
#output_rows.append(header)

summary = pd.read_csv(path_generated + '/gt_glycans_filt_1.0.tsv', sep="\t")
print(summary.head())

disease_list = summary['cancer'].unique()

for disease in disease_list:
	output_row = []
	summary_disease = summary.loc[summary['cancer'] == disease]
	uniprot_list = summary_disease['uniprotkb_ac'].unique()
	genes_count = len(uniprot_list)
        genes_down = summary_disease.loc[summary['cancer_trend'] == "DOWN"]
	genes_down_count = len(genes_down)
	genes_up = summary_disease.loc[summary['cancer_trend'] == "UP"]
	genes_up_count = len(genes_up)
	glycan_list = summary_disease['glytoucan_ac'].unique()
	glycan_count = len(glycan_list)
	position = summary_disease['residue_id'].unique()
	position_count = len(position)
	output_row.append(disease)
	output_row.append(genes_count)
	output_row.append(genes_down_count)
	output_row.append(genes_up_count)
        output_row.append(glycan_count)
        output_row.append(position_count)
        output_rows.append(output_row)

with open(path_generated + '/glycan_counts_by_cancer.tsv', "w") as outfile:
	writer = csv.writer(outfile, delimiter = "\t")
	writer.writerows(output_rows)

print('Done')

