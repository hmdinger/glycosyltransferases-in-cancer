#!/bin/bash

###### SETUP CONFIGURATION #######
# Set the available cancer types from sort | uniq

#cancers="BLCA BRCA CESC COAD ESCA HNSC KICH KIRC KIRP LIHC LUAD LUSC PAAD PRAD READ SARC STAD THCA UCEC"
#echo -e "Cancer \tTotSigGenes \tTotSigGTGenes \tTotUp05Genes \tUpGT05Genes \tTotDown05Genes \tDownGT05Genes \tTot01Genes \tGT01Genes \tTotUp01Genes \tUpGT01Genes \tTotDown01Genes \tDownGT01Genes \tTotFC1UpGenes \tUpGTFC1Genes \tTotFC1DownGenes \tDownGTFC1Genes \tTotSig75UpGenes \tGTSig75UpGenes \tTotSig75DownGenes \tGTSig75DownGenes \tTotSig90UpGenes \tGTSig90UpGenes \tTotSig90DownGenes \tGTSig90DownGenes" >> total_counts.txt

awk -F'\t' '{print $2}' ensembl_gts_list_mm.txt > mouse_ensembl_gts_list.txt

mouse_ensembl_gts=$(<mouse_ensembl_gts_list.txt)

for gt in $mouse_ensembl_gts; do

 #Print rows from master table
        echo "Retrieving mouse rows"

        awk -F'\t' '$1 == "'"${gt}"'"' Mus_musculus_UBERON:0000113_RNA_SEQ.tsv >> filtered_bgee_mouse.txt
done


