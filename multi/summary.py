#This program downloads data from UniProtKB for human and mouse 1:1 glycosyltransferase orthologs with differential expression in cancer

import sys
import urllib
import os
import cons_config
import glob
import csv
import pandas as pd
from array import array

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

gts = []
#with open('ortholog_gts.tsv', "r") as gt_file:

#gts = list(csv.reader(open('ortholog_gts.tsv')))

gt_dict = {}
with open('orthologs_info_filtered.tsv', "r") as gt_file:
        gt_rows = csv.DictReader(gt_file, delimiter="\t")
        next(gt_rows)
        for d in gt_rows:
		#gt_dict.setdefault(d['human_uniprotkb_ac'],[]).append(d['human_gene'])
		gt_dict.setdefault(d['human_uniprotkb_ac'],d['human_gene'])


#                if row[0] in gts :
#                        continue
#                else:
#                        gts.append(row[0])


print(gt_dict)
gts = gt_dict.keys()
genes = gt_dict.values()

print('There are ' + str(len(gts)) + ' gts with ortholog and expression information.')

print('Defining tissues')
output_rows = []
tissues = []
with open('ortholog_tissues.tsv', "r") as tissue_file:
        tissue_list = csv.reader(tissue_file, delimiter = "\t")
	#tissue_list = tissue_file.readlines()
	for row in tissue_list:
		tissues.append(row[0])
print(tissues)
#print(tissue_list)

print('There are ' + str(len(tissues)) + ' tissues with ortholog and expression information.')

#Summarize scores by gene
print('Making gene summary table with absent calls in score')
output_rows = []
header = ['uniprotkb_ac', 'gene', 'gene_score', 'number_tissues', 'average_score']
output_rows.append(header)

for gt in gts:
	output_row = []
	gene = gt_dict[gt]
	print(gt, gene)
	gene_score = 0
	tissue_list = []
	gt_tissue = []
	tissue_count = 0
	average = 0

	with open('conservation_master.tsv', "r") as cons_file:
        	rows = csv.reader(cons_file, delimiter="\t")
		for row in rows:
			cons_gt = row[0]
			gt_tissue = row[9]
			if cons_gt == gt:
				print(row)
				gene_score += int(row[15])
				if gt_tissue in tissue_list:
					continue
				else:
					tissue_list.append(gt_tissue)
			else: 
				continue
		tissue_count = len(tissue_list)
		if tissue_count > 0:
			average = gene_score/float(tissue_count)
			output_row.append(gt)
			output_row.append(gene)
			output_row.append(gene_score)
			output_row.append(tissue_count)
			output_row.append(average)
			output_rows.append(output_row)
		else:
			continue

with open('summary_by_gene.tsv', "w") as outfile:
	writer = csv.writer(outfile, delimiter = "\t")
	writer.writerows(output_rows)

print('Making gene summary table with absent calls in score')
output_rows = []
header = ['uniprotkb_ac', 'gene', 'gene_score', 'number_tissues', 'average_score']
output_rows.append(header)

for gt in gts:
        output_row = []
        gene = gt_dict[gt]
        print(gt, gene)
        gene_score = 0
        tissue_list = []
        gt_tissue = []
        tissue_count = 0
        average = 0

        with open('no_absent_conservation.tsv', "r") as cons_file:
                rows = csv.reader(cons_file, delimiter="\t")
                for row in rows:
                        cons_gt = row[0]
                        gt_tissue = row[9]
                        if cons_gt == gt:
                                print(row)
				gene_score += int(row[15])
                                if gt_tissue in tissue_list:
                                        continue
                                else:
                                        tissue_list.append(gt_tissue)
                        else:
                                continue
                tissue_count = len(tissue_list)
                if tissue_count > 0:
                        average = gene_score/float(tissue_count)
                        output_row.append(gt)
                        output_row.append(gene)
                        output_row.append(gene_score)
                        output_row.append(tissue_count)
                        output_row.append(average)
                        output_rows.append(output_row)
                else:
                        continue

with open('summary_by_gene_no_absent.tsv', "w") as outfile:
        writer = csv.writer(outfile, delimiter = "\t")
        writer.writerows(output_rows)

#Summarize score by tissues
print('Making tissue summary table for all calls')
output_rows = []
header = ['tissue', 'tissue_score', 'number_genes', 'average_score']
output_rows.append(header)

for tissue in tissues:
        output_row = []
	print(tissue)
        tissue_score = 0
        gene_list = []
        gt_gene = []
        gene_count = 0
        average = 0

        with open('conservation_master.tsv', "r") as cons_file:
                rows = csv.reader(cons_file, delimiter="\t")
                for row in rows:
                        cons_gt = row[0]
			gt_gene = row[1]
                        cons_tissue = row[9]
                        if cons_tissue == tissue:
                                print(row)
                                tissue_score += int(row[15])
                                if gt_gene in gene_list:
                                        continue
                                else:
                                        gene_list.append(gt_gene)
                        else:
                                continue
                gene_count = len(gene_list)
                if gene_count > 0:
			average = tissue_score/float(gene_count)
                        output_row.append(tissue)
                        output_row.append(tissue_score)
                        output_row.append(gene_count)
                        output_row.append(average)
                        output_rows.append(output_row)
                else:
                        continue

with open('summary_by_tissue.tsv', "w") as outfile:
        writer = csv.writer(outfile, delimiter = "\t")
        writer.writerows(output_rows)

print('Making tissue summary table with absent calls filtered out')
output_rows = []
header = ['tissue', 'tissue_score', 'number_genes', 'average_score']
output_rows.append(header)

for tissue in tissues:
        output_row = []
        print(tissue)
        tissue_score = 0
        gene_list = []
        gt_gene = []
        gene_count = 0
        average = 0

        with open('no_absent_conservation.tsv', "r") as cons_file:
                rows = csv.reader(cons_file, delimiter="\t")
                for row in rows:
                        cons_gt = row[0]
                        gt_gene = row[1]
                        cons_tissue = row[9]
                        if cons_tissue == tissue:
                                print(row)
                                tissue_score += int(row[15])
                                if gt_gene in gene_list:
                                        continue
                                else:
                                        gene_list.append(gt_gene)
                        else:
                                continue
                gene_count = len(gene_list)
                if gene_count > 0:
			average = tissue_score/float(gene_count)
                        output_row.append(tissue)
                        output_row.append(tissue_score)
                        output_row.append(gene_count)
                        output_row.append(average)
                        output_rows.append(output_row)
                else:
                        continue

with open('summary_by_tissue_no_absent.tsv', "w") as outfile:
        writer = csv.writer(outfile, delimiter = "\t")
	writer.writerows(output_rows)	

print('Done')


