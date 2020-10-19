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
        for row in gt_rows:
                if row[0] in ortho_gts :
                        continue
                else:
                        ortho_gts.append(row[0])

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
for gt in ortho_gts :
        cmd = """awk -F',' '$2 == \"\\\"""" + gt + """\\\"\"' """ + oncomx_now + """/mouse_normal_expression_custom.csv >> """ + path_intermediate + """/mouse_normal_custom_filtered.csv"""
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
	"THCA":	"thyroid",
	"UCEC":	"uterus"
}
print(map)
print(map["BLCA"])

#Make master expression conservation table
print('Generating expression conservation table')

output_rows = []
header = ['human_uniprotkb_ac', 'human_gene', 'mouse_uniprotkb_ac', 'mouse_gene', 'cancer', 'cancer_trend', 'log2FC', 'p_value', 'uberon_id', 'uberon_name', 'human_healthy_expression', 'human_score', 'mouse_healthy_expression', 'mouse_score']
output_rows.append(header)

print('Retrieving gts with mouse orthologs and cancer data')
with open('orthologs_info_filtered.tsv', 'r') as orth_file:
	orthologs = csv.reader(orth_file, delimiter="\t")
	next(orthologs)
	k = 0
	for ortholog in orthologs:
		k += 1
		hs_ac = ortholog[0]
		hs_gene = ortholog[1]
		ms_ac = ortholog[2]
		ms_gene = ortholog[3]
#		output_row.append(hs_ac)
#                output_row.append(hs_gene)
#                output_row.append(ms_ac)
#                output_row.append(ms_gene)

		if k == 1 :

			print('Retrieving differential expression details for ' + hs_ac)
			with open(path_diff_gen + '/de_gts_master_filtered.tsv', "r") as de_file:

			#Example row:
			#UniProtKB_AC    RefSeq  Gene    log2FoldChange  p_value adjusted p_value        Significant     Trend   TCGA Cancer     Cancer Ontology PMID    			#UBERON_ID       PtsTrend        PtsUp   PtsDown TotPts  Percent
			#A0PJZ3  NP_001073862.1  GXYLT2  -1.0    1.23e-02        3.19e-02        Yes     Down    BLCA    DOID:11054 / Urinary bladder cancer [UBC]       			#-       UBERON:0001255 12       7       12      19      63.16		

				lines = csv.reader(de_file, delimiter="\t")
				for line in lines:
					output_row = []
					disease = line[8]
					trend = line[7]
					log2fc = line[3]
					p = line[5]
					j = 0
					if line[0] == hs_ac:
						j +=1
						print(j)
						print('Retrieving human normal expression details for ' + hs_ac + ' in ' + disease + '/' + map[disease])
						if j == 1:
							with open(path_intermediate + '/human_normal_custom_filtered.csv', "r") as human_file:
								hs_rows = csv.reader(human_file, delimiter=',')

						#Example row:
						#"gene_symbol","uniprotkb_ac","uberon_anatomical_id","uberon_anatomical_name_long","uberon_anatomical_name_short","expre						#ssion_level_anatomical_relative","expression_score"
						#"TSPAN6","O43657","UBERON:0002107","liver","liver","HIGH","96.3"
								yes_count = 0
								no_count = 0
								for row in hs_rows:
									output_row = []
									hs_ac_bgee = row[1].strip('\"')
									uberon_id = row[2].strip('\"')
									uberon_name = row[4].strip('\"')
									hs_level = row[5].strip('\"')
		                					hs_score = row[6].strip('\"')
									if ((hs_ac_bgee == hs_ac) and (uberon_name == map[disease])):
										print(hs_ac_bgee)
										print(uberon_name)
										print(row)
										yes_count += 1

										output_row.append(hs_ac)
					                                        output_row.append(hs_gene)
                	                				        output_row.append(ms_ac)
					                                        output_row.append(ms_gene)
                                					        output_row.append(disease)
				        	                                output_row.append(trend)
				                	                        output_row.append(log2fc)
                                					        output_row.append(p)
										output_row.append(uberon_id)
										output_row.append(uberon_name)
										output_row.append(hs_level)
										output_row.append(hs_score)
										output_rows.append(output_row)
									else:
                                                                                no_count += 1

								#else:
									#continue
							print("Yes: " + str(yes_count))
							print("No: " + str(no_count))
							total = yes_count + no_count
							print(total)
						
						else:
							continue

					else:
						continue

		else:
			continue

with open('conservation_master.tsv', "wb") as out_file:
	writer = csv.writer(out_file, delimiter = "\t")
	writer.writerows(output_rows)

#Going to make the modified bgee tables and then join
#                if row[0] in go_data:
                    # The AC is found in the go data!
                    # Get terms in semi-colon separated form
 #                   go_terms = ";".join(go_data[row[0]])
                    # Add to row as new column
  #                  row.append(go_terms)
                    # Write out to file
   #                 output_file.writerow(row)
	
#Print completion of current task
print('Done')


