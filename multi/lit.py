#This program generates a table for DEGTs literature mining evidence

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
			#ac_stripped = row[0].replace('"','')
                        #ac_split = stripped.split('-', 1)[0]
                        #gene_stripped = row[3].replace('"','')
                        #gene_split = stripped.split('-', 1)[0]
                        gt_dict.setdefault(d['uniprotkb_canonical_ac'],d['gene_symbol'])


print(gt_dict)
gts = gt_dict.keys()
genes = gt_dict.values()

for gt in gts:
	gt_stripped = gt.replace('"','')
	gt_split = gt.split('-', 1)[0]
	gt_dict[gt_split] = gt_dict.pop(gt)

print(gt_dict)
gts = gt_dict.keys()
genes = gt_dict.values()

#Make dictionary for disease mapping for disease mentions
copy = 'cp ' + manual_path + '/cdo_copy.txt cdo_slim.tsv'
os.system(copy)

#check = []
disease = {}
with open('cdo_slim.tsv', "r") as disease_file:
        rows = csv.DictReader(disease_file, delimiter = "\t")
        next(rows)
        for d in rows:
 #               if d in check:
  #                      continue
   #             else:
    #                    check.append(d)
		disease.setdefault(d['Child_Organ_system'],d['Slim'])


print(disease)
child_terms = disease.keys()
slim_terms = disease.values()

for child in child_terms:
	child_split = child.split('/',1)[0]
	child_strip = child_split.replace(' ','')
	disease[child_strip] = disease.pop(child)
#	print(map["LUAD"])

print(disease)

for child, slim in disease.items():
	#slim_replace = slim.replace('/','/\ ')
	slim_split = slim.split('/',1)[0]
	slim_strip = slim_split.replace(' ','')
	disease[child] = slim_strip

print(disease)
children = disease.keys()
slims = disease.values()


#######################################
### Make master table for all genes ###
#######################################

print('Printing master lit mining file')
output_rows = []
header = ['Entry','GeneName','PMID','DOID','DOID_Name','DiseaseMention','DiseaseDetectedFrom','SentenceType','Expression_OR_Activation','Sentence']
output_rows.append(header)
for gene in genes:
	gene_path = '/home/hmhamilt/phd_dissertation/multi/intermediate/dexter/279GT_CSV/' + gene + '.csv'
	if os.path.exists(gene_path):
		with open(gene_path, "r") as gene_file:
			rows = csv.reader(gene_file)
			next(rows)
			for row in rows:
				if row == "":
					continue
				else:
					output_rows.append(row)
	else:
		continue

with open('all_genes_lit.tsv', "w") as out_file:
	writer = csv.writer(out_file, delimiter = "\t")
	writer.writerows(output_rows)


###########################################
### Print combined table using old file ###
###########################################

#Make master DEGT with literature table
print('Generating expression literature table')

de_file = pd.read_csv(path_diff_gen + '/de_gts_master_filtered.tsv', sep="\t")
lit_file = pd.read_csv(path_intermediate + '/all_genes_lit.tsv', sep="\t")
corr_file = pd.read_csv(path_intermediate + '/do_tcga.tsv', sep="\t")

print(de_file.head())
print(lit_file.head())

#Define cancer-cancer map from mentions in literature
print('Making a dictionary for cancer terms and their literature mentions')

common1 = pd.merge(de_file, corr_file, how='left', left_on=['TCGA Cancer'], right_on=['tcga'])
print(common1.head())

print(len(de_file.index))
print(len(common1.index))

#common1_null = common1.loc[common1['mirna'].isnull()]
#print(common1_null.head())
#print(len(common1_null.index))

common1_full = common1.dropna(subset=['tcga'])
print(common1_full.head())
print(len(common1_full.index))

common2 = pd.merge(common1, lit_file, how='left', left_on=['UniProtKB_AC','do'], right_on=['Entry','DOID_Name'])
print(common2.head())
print(len(common2.index))

common2_full = common2.dropna(subset=['Entry'])
print(common2_full.head())
print(len(common2_full.index))

for col in common2_full.columns:
	print(col)

better_headers = common2_full.rename(columns = {'UniProtKB_AC': 'uniprotkb_ac', 'Gene': 'gene', 'TCGA Cancer': 'cancer', 'Trend': 'cancer_trend', 'adjusted p_value': 'adjusted_p_value', 'log2FoldChange': 'log2fc', 'PMID_y': 'pmid', 'DiseaseMention': 'disease_mention', 'DiseaseDetectedFrom': 'detected_from', 'SentenceType': 'sentence_type', 'Expression_OR_Activation': 'action', 'Sentence': 'sentence'})

print(better_headers.head())

better_headers.to_csv(path_generated + '/gt_lit_old.tsv', sep='\t', columns=['uniprotkb_ac', 'gene', 'cancer', 'cancer_trend', 'adjusted_p_value', 'log2fc', 'pmid', 'disease_mention', 'detected_from', 'sentence_type', 'action', 'sentence'], index=False)

###########################################
### Print combined table using new file ###
###########################################

copy = 'cp ' + oncomx_now + '/human_cancer_expression_literature.csv oncomx_literature.csv'
os.system(copy)

#Make master DEGT with literature table for all genes all cancers
print('Generating expression literature table from newer dataset')

de_file = pd.read_csv(path_diff_gen + '/de_gts_master_filtered.tsv', sep="\t")
lit_file = pd.read_csv(path_intermediate + '/oncomx_literature.csv', sep=",")
print(de_file.head())
print(lit_file.head())
print(len(lit_file.index))

de_file['doid'] = de_file['Cancer Ontology'].str.split(' /',1).str[0]
print(de_file['doid'].head())
de_file.to_csv('test1_newcol.tsv', sep='\t')

print(lit_file['doid'].head())
lit_file['slim'] = lit_file['doid'].map(disease)
print(lit_file['slim'].head())
lit_file.to_csv('test2_domap.tsv', sep='\t')
print(len(lit_file.index))

hits = lit_file.dropna(subset=['slim'])
print(hits['slim'])
print(len(hits.index))

common = pd.merge(de_file, lit_file, how='left', left_on=['UniProtKB_AC','doid'], right_on=['uniprotkb_ac','slim'])
print(common.head())

common_full = common.dropna(subset=['uniprotkb_ac','slim'])
print(common_full.head())

better_headers = common_full.rename(columns = {'UniProtKB_AC' : 'uniprotkb_ac','Gene': 'gene', 'TCGA Cancer': 'cancer', 'Trend': 'cancer_trend', 'adjusted p_value': 'adjusted_p_value', 'log2FoldChange': 'log2fc', 'PMID_y': 'pmid', 'uniprotkb_ac' : 'ac'})

better_headers.to_csv(path_generated + '/gt_lit_new.tsv', sep='\t', columns=['uniprotkb_ac', 'gene', 'cancer', 'cancer_trend', 'adjusted_p_value', 'log2fc', 'pmid', 'disease_mention', 'disease_extracted_from', 'expression_change', 'sample_1', 'sample_2', 'is_same_patient', 'sentence_type', 'sentence_text'], index=False)


##########################################
### Make high value table for new file ###
##########################################
print('Generating high value table based on DE and literature trends agreeing')

output_rows = []
with open(path_generated + '/gt_lit_new.tsv', "r") as lit_file:
	rows = csv.reader(lit_file, delimiter="\t")
	needheader=True
	for row in rows:
		if needheader:
			needheader=False
			output_rows.append(row)
		elif ((row[3] == "Down") and (row[9] == "DOWN")) or ((row[3] == "Up") and (row[9] == "UP")):
			output_rows.append(row)
		else:
			continue

with open(path_generated + '/gt_lit_same.tsv', "w") as outfile:
	writer = csv.writer(outfile, delimiter = "\t")
	writer.writerows(output_rows)
		
###################################
### Summarize old table by gene ###
###################################

print('Summarize counts of literature')
output_rows = []
header = ['uniprotkb_ac','gene','disease_count','pmid_count','sentence_count']
output_rows.append(header)

summary = pd.read_csv('all_genes_lit.tsv', sep="\t")
print(summary.head())
#def same(s):
#        if (s['gt_cancer_trend'] == s['mirna_cancer_trend']):
#                return 'same'
#        else:
#                return 'different'

#summary['is_same'] = summary.apply(same, axis=1)

for gt in gts:
        gene = gt_dict[gt]
        if gt in summary.values:
		output_row = []
                summary_gt = summary.loc[summary['Entry'] == gt]
                disease_list = summary_gt['DOID'].unique()
                disease_count = len(disease_list)
		pmid_list = summary_gt['PMID'].unique()
		pmid_count = len(pmid_list)
		sentence_list = summary_gt['Sentence'].unique()
		sentence_count = len(sentence_list)
		output_row.append(gt)
		output_row.append(gene)
		output_row.append(disease_count)
		output_row.append(pmid_count)
		output_row.append(sentence_count)
		output_rows.append(output_row)
	else:
		continue

with open(path_generated + "/counts_old.tsv", "w") as out_file:
	writer = csv.writer(out_file, delimiter = "\t")
	writer.writerows(output_rows)


###################################
### Summarize new table by gene ###
###################################

output_rows = []
header = ['uniprotkb_ac','gene','cancer','pmid_count','sentence_count','sentence_up_count','sentence_down_count']
output_rows.append(header)

summary = pd.read_csv('oncomx_literature.csv', sep=",")
print(summary.head())

#def same(s):
#        if (s['gt_cancer_trend'] == s['mirna_cancer_trend']):
#                return 'same'
#        else:
#                return 'different'

#summary['is_same'] = summary.apply(same, axis=1)

for gt in gts:
        gene = gt_dict[gt]
        if gt in summary.values:
                summary_gt = summary.loc[summary['uniprotkb_ac'] == gt]
		disease_list = summary_gt['doid'].unique()
		for disease in disease_list:
			output_row = []
		        disease_summary = summary_gt.loc[summary['doid'] == disease]
		        pmid_list = disease_summary['pmid'].unique()
			pmid_count = len(pmid_list)
		        sentence_list = disease_summary['sentence_text'].unique()
			sentence_up = disease_summary.loc[summary['expression_change'] == "UP"]
		        sentence_up_count = len(sentence_up)
			sentence_down = disease_summary.loc[summary['expression_change'] == "DOWN"]
			sentence_down_count = len(sentence_down)
		        output_row.append(gt)
                	output_row.append(gene)
		        output_row.append(disease)
                	output_row.append(pmid_count)
			output_row.append(sentence_up_count)
			output_row.append(sentence_down_count)	
                	output_rows.append(output_row)
        else:
                continue

with open(path_generated + "/counts_gene_cancer.tsv", "w") as out_file:
        writer = csv.writer(out_file, delimiter = "\t")
        writer.writerows(output_rows)

#################################
### Summarize counts by caner ###
#################################

output_rows = []
header = ['cancer','total_genes','down_genes','up_genes','total_gts','down_gts','up_gts','pmids','sentences']
output_rows.append(header)

summary = pd.read_csv('oncomx_literature.csv', sep=",")
print(summary.head())

disease_list = summary['doid'].unique()

for disease in disease_list:
	output_row = []
	summary_disease = summary.loc[summary['doid'] == disease]
	uniprot_list = summary_disease['uniprotkb_ac'].unique()
	genes_count = len(uniprot_list)
        genes_down = summary_disease.loc[summary['expression_change'] == "DOWN"]
	genes_down_count = len(genes_down)
	genes_up = summary_disease.loc[summary['expression_change'] == "UP"]
	genes_up_count = len(genes_up)
	gts_count = 0
	gts_up_count = 0
	gts_down_count = 0
	for gt in gts:
		if gt in uniprot_list:
			gts_count +=1
		if gt in genes_down:
			gts_down_count +=1
		if gt in genes_up:
			gts_up_count +=1
	pmid_list = summary_disease['pmid'].unique()
	pmid_count = len(pmid_list)
	sentences = summary_disease['sentence_text'].unique()
	sentence_count = len(sentences)
	output_row.append(disease)
	output_row.append(genes_count)
	output_row.append(genes_down_count)
	output_row.append(genes_up_count)
	output_row.append(gts_count)
	output_row.append(gts_down_count)
	output_row.append(gts_up_count)
        output_row.append(pmid_count)
        output_row.append(sentence_count)
        output_rows.append(output_row)

with open(path_generated + '/counts_by_cancer.tsv', "w") as outfile:
	writer = csv.writer(outfile, delimiter = "\t")
	writer.writerows(output_rows)

quit()

print('Filtering for increased foldchange to identify high-value pairs')
thresholds = [1.0, 2.0, 3.0]
for threshold in thresholds:
	output_rows = []
	with open('similar_gt_scrna.tsv', "r") as similar_file:
        	sims = csv.reader(similar_file, delimiter="\t")
	        needheader=True
   		for row in sims:
                	if needheader:
                        	needheader=False
	                        output_rows.append(row)
        	                continue
                	elif (((float(row[5])) >= threshold) == True) or (((float(row[5])) <= -threshold) == True) :
                        	output_rows.append(row)
	                else:
        	                continue
	with open("high_value_gt_scrna_similar_pairs_" + str(threshold) + ".tsv", "w") as out_file:
        	writer = csv.writer(out_file, delimiter='\t')
	        writer.writerows(output_rows)

quit()

print('Done')

