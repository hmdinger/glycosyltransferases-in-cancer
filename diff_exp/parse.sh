#!/bin/bash

path_downloads=$(python -c "import diff_exp_config,os; path_downloads = diff_exp_config.path_downloads; print(path_downloads)")
path_intermediate=$(python -c "import diff_exp_config,os; path_intermediate = diff_exp_config.path_intermediate; print(path_intermediate)")
path_generated=$(python -c "import diff_exp_config,os; path_generated = diff_exp_config.path_generated; print(path_generated)")
bioxpress_now=$(python -c "import diff_exp_config,os; bioxpress_now = diff_exp_config.bioxpress_now; print(bioxpress_now)")

cd ${path_intermediate}

pwd

sed -e 's/[\r\n]$//g' cancers.tsv > cancers_clean.tsv
sed -e 's/[\r\n]$//g' gts.tsv > gts_clean.tsv


cancers=$(<cancers_clean.tsv)
delete=-
echo ${cancers[@]/$delete}
cancers=("${cancers[@]/$delete}")

gts=$(<gts_clean.tsv)

echo """${cancers[*]}"""
echo """${gts[*]}"""


cp ${bioxpress_now}/BioXpress_interface_overall_final_v2.0.csv v2_copy.csv
sed 's/,/\t/g' v2_copy.csv > v2_copy.tsv

echo -e "Cancer \tTotSigGenes \tTotSigGTGenes \tTotUp05Genes \tUpGT05Genes \tTotDown05Genes \tDownGT05Genes \tTot01Genes \tGT01Genes \tTotUp01Genes \tUpGT01Genes \tTotDown01Genes \tDownGT01Genes \tTotFC1UpGenes \tUpGTFC1Genes \tTotFC1DownGenes \tDownGTFC1Genes \tTotSig75UpGenes \tGTSig75UpGenes \tTotSig75DownGenes \tGTSig75DownGenes \tTotSig90UpGenes \tGTSig90UpGenes \tTotSig90DownGenes \tGTSig90DownGenes" > ${path_generated}/total_counts.txt

echo -e "PtsTrend \tTotPts \tPercent" > BioXpress_patients.tsv
awk -F'\t' '{if (NR!=1) {print}}' v2_copy.tsv | awk -F'\t' '{print $11}' | awk -F'[/()]' '{print $1 "\t" $2 "\t" $3}' >> BioXpress_patients.tsv
paste v2_copy.tsv BioXpress_patients.tsv > BioXpress_all_with_patients.tsv

head -n 1 BioXpress_all_with_patients.tsv > ${path_generated}/Differential_expression_gts_master.tsv
awk -F'\t' '$7 == "Yes"' BioXpress_all_with_patients.tsv > BioXpress_all_sig.tsv


for cancer in $cancers; do
	# Split patient counts and generate distinct cancer tables

	#Example of table
	#1UniProtKB_AC   2 RefSeq 3 Gene    4log2FoldChange  5p_value adjusted 6p_value       7 Significant    8 Trend  9 TCGA Cancer    10 Cancer Ontology 11#Patients      12 Data Source   13 PMID     14UBERON_ID
	#Q6UXB8  NP_001186088.1; NP_699201.2; XP_005248974.1; XP_011512677.1     PI16    -7.72   4.13e-46        8.02e-42        Yes     Down    BLCA    DOID:11054 / Urinary bladder cancer [UBC]       19/19(100.0)    RNASeqV2        -       UBERON:0001255

	echo -e "PtsTrend \tTotPts \tPercent" > BioXpress_patients.tsv
        echo "***awk -F'\t' '$9 == '"${cancer}"'' BioXpress_all_with_patients.tsv > ${cancer}_bioxpress.tsv"
        awk -F'\t' '$9 == "'"${cancer}"'"' BioXpress_all_with_patients.tsv > ${cancer}_bioxpress.tsv

        #Filter by significance
	echo "***awk -F'\t' '$7 == "Yes"' ${cancer}_bioxpress.tsv > ${cancer}_bioxpress_sig.tsv"
	awk -F'\t' '$7 == "Yes"' ${cancer}_bioxpress.tsv > ${cancer}_bioxpress_sig.tsv

        #Remove entries without UniProt accessions
        echo "***awk -F'\t' '$1 != "-"' ${cancer}_bioxpress_sig.tsv > ${cancer}_bioxpress_sig_uni.tsv"
	awk -F'\t' '$1 != "-"' ${cancer}_bioxpress_sig.tsv > ${cancer}_bioxpress_sig_uni.tsv

        #Find GTs at p < 0.05
	echo "***${cancer}: mapping GTs"
	for gt in $gts; do

               awk -F'\t' '$1 == "'"${gt}"'"' ${cancer}_bioxpress_sig_uni.tsv >> ${cancer}_gts.tsv

        done

	#Find All Up and Down genes
        awk -F'\t' '$8 == "Up"' ${cancer}_bioxpress_sig_uni.tsv > ${cancer}_all_up_05.tsv
        awk -F'\t' '$8 == "Down"' ${cancer}_bioxpress_sig_uni.tsv > ${cancer}_all_down_05.tsv

	#Find Up and Down .05 GTS
        awk -F'\t' '$8 == "Up"' ${cancer}_gts.tsv > ${cancer}_05_up_gts.tsv
        awk -F'\t' '$8 == "Down"' ${cancer}_gts.tsv > ${cancer}_05_down_gts.tsv

	#Find All Up and Down.01 Genes
	awk -F'\t' -v thresh01=0.01 '$6 < thresh01' ${cancer}_bioxpress_sig_uni.tsv > ${cancer}_all_01.tsv 
	awk -F'\t' '$8 == "Up"' ${cancer}_all_01.tsv > ${cancer}_all_up_01.tsv
        awk -F'\t' '$8 == "Down"' ${cancer}_all_01.tsv > ${cancer}_all_down_01.tsv

	#Find GTs Up and Down  at p < 0.01
	awk -F'\t' -v thresh01=0.01 '$6 < thresh01' ${cancer}_gts.tsv > ${cancer}_01_gts.tsv  
	awk -F'\t' '$8 == "Up"' ${cancer}_01_gts.tsv > ${cancer}_01_up_gts.tsv
	awk -F'\t' '$8 == "Down"' ${cancer}_01_gts.tsv > ${cancer}_01_down_gts.tsv

	#Find All .01 FC > 1 at 75 and 90
	awk -F'\t' -v fcup1=1.0 '$4 >= fcup1' ${cancer}_all_up_01.tsv > ${cancer}_all_up_fc.tsv
        awk -F'\t' -v fcdown1=-1.0 '$4 <= fcdown1' ${cancer}_all_down_01.tsv > ${cancer}_all_down_fc.tsv
        awk -F'\t' -v sigpt75=75.0 '$17 >= sigpt75' ${cancer}_all_up_fc.tsv > ${cancer}_all_up_sig75.tsv
        awk -F'\t' -v sigpt75=75.0 '$17 >= sigpt75' ${cancer}_all_down_fc.tsv > ${cancer}_all_down_sig75.tsv
        awk -F'\t' -v sigpt90=90.0 '$17 >= sigpt90' ${cancer}_all_up_fc.tsv > ${cancer}_all_up_sig90.tsv
        awk -F'\t' -v sigpt90=90.0 '$17 >= sigpt90' ${cancer}_all_down_fc.tsv > ${cancer}_all_down_sig90.tsv

        #Find GTs .01 FC > 1 at 75 and 90
        awk -F'\t' -v fcup1=1.0 '$4 >= fcup1' ${cancer}_01_up_gts.tsv > ${cancer}_gts_up_fc.tsv
        awk -F'\t' -v fcdown1=-1.0 '$4 <= fcdown1' ${cancer}_01_down_gts.tsv > ${cancer}_gts_down_fc.tsv
        awk -F'\t' -v sigpt75=75.0 '$17 >= sigpt75' ${cancer}_gts_up_fc.tsv > ${cancer}_gts_up_sig75.tsv
        awk -F'\t' -v sigpt75=75.0 '$17 >= sigpt75' ${cancer}_gts_down_fc.tsv > ${cancer}_gts_down_sig75.tsv
        awk -F'\t' -v sigpt90=90.0 '$17 >= sigpt90' ${cancer}_gts_up_fc.tsv > ${cancer}_gts_up_sig90.tsv
        awk -F'\t' -v sigpt90=90.0 '$17 >= sigpt90' ${cancer}_gts_down_fc.tsv > ${cancer}_gts_down_sig90.tsv

	#Make individual significant lists
        awk -F'\t' '{print $1}' ${cancer}_05_up_gts.tsv > ${cancer}_05_up_list.tsv
        awk -F'\t' '{print $1}' ${cancer}_05_down_gts.tsv > ${cancer}_05_down_list.tsv
	awk -F'\t' '{print $1}' ${cancer}_01_up_gts.tsv > ${cancer}_01_up_list.tsv
        awk -F'\t' '{print $1}' ${cancer}_01_down_gts.tsv > ${cancer}_01_down_list.tsv
        awk -F'\t' '{print $1}' ${cancer}_gts_up_fc.tsv > ${cancer}_gts_up_fc_list.tsv
        awk -F'\t' '{print $1}' ${cancer}_gts_down_fc.tsv > ${cancer}_gts_down_fc_list.tsv
	awk -F'\t' '{print $1}' ${cancer}_gts_up_sig75.tsv > ${cancer}_gts_up_sig75_list.tsv
        awk -F'\t' '{print $1}' ${cancer}_gts_down_sig75.tsv > ${cancer}_gts_down_sig75_list.tsv
        awk -F'\t' '{print $1}' ${cancer}_gts_up_sig90.tsv > ${cancer}_gts_up_sig90_list.tsv
        awk -F'\t' '{print $1}' ${cancer}_gts_down_sig90.tsv > ${cancer}_gts_down_sig90_list.tsv

        #Count genes in total dataset
        echo "***${cancer}"
        allsig="${cancer}_bioxpress_sig_uni.tsv"
        allsigcounts=$(wc -l ${allsig} | awk '{print $1}')
	allup="${cancer}_all_up_05.tsv"
	allupcounts=$(wc -l ${allup} | awk '{print $1}')
	alldown="${cancer}_all_down_05.tsv"
	alldowncounts=$(wc -l ${alldown} | awk '{print $1}')
	all01="${cancer}_all_01.tsv"
	all01counts=$(wc -l ${all01} | awk '{print $1}')
	all01up="${cancer}_all_up_01.tsv"
	all01upcounts=$(wc -l ${all01up} | awk '{print $1}')
	all01down="${cancer}_all_down_01.tsv"
	all01downcounts=$(wc -l ${all01down} | awk '{print $1}')
	allfcup="${cancer}_all_up_fc.tsv"
	allfcupcounts=$(wc -l ${allfcup} | awk '{print $1}')
        allfcdown="${cancer}_all_down_fc.tsv"
        allfcdowncounts=$(wc -l ${allfcdown} | awk '{print $1}')
        allsig75up="${cancer}_all_up_sig75.tsv"
        allsig75upcounts=$(wc -l ${allsig75up} | awk '{print $1}')
        allsig75down="${cancer}_all_down_sig75.tsv"
        allsig75downcounts=$(wc -l ${allsig75down} | awk '{print $1}')
        allsig90up="${cancer}_all_up_sig90.tsv"
        allsig90upcounts=$(wc -l ${allsig90up} | awk '{print $1}')
        allsig90down="${cancer}_all_down_sig90.tsv"
        allsig90downcounts=$(wc -l ${allsig90down} | awk '{print $1}')
        cangts="${cancer}_gts.tsv"
        cangtcounts=$(wc -l ${cangts} | awk '{print $1}')
	can05upgts="${cancer}_05_up_gts.tsv"
	can05upgtcounts=$(wc -l ${can05upgts} | awk '{print $1}')
	can05downgts="${cancer}_05_down_gts.tsv"
	can05downgtcounts=$(wc -l ${can05downgts} | awk '{print $1}')
        can01gts="${cancer}_01_gts.tsv"
        can01gtcounts=$(wc -l ${can01gts} | awk '{print $1}')
        can01upgts="${cancer}_01_up_gts.tsv"
        can01upgtcounts=$(wc -l ${can01upgts} | awk '{print $1}')
        can01downgts="${cancer}_01_down_gts.tsv"
        can01downgtcounts=$(wc -l ${can01downgts} | awk '{print $1}')
	canfcup="${cancer}_gts_up_fc.tsv"
	canfcupcounts=$(wc -l ${canfcup} | awk '{print $1}')
        canfcdown="${cancer}_gts_down_fc.tsv"
        canfcdowncounts=$(wc -l ${canfcdown} | awk '{print $1}')
	cansig75up="${cancer}_gts_up_sig75.tsv"
	cansig75upcounts=$(wc -l ${cansig75up} | awk '{print $1}')
        cansig75down="${cancer}_gts_down_sig75.tsv"
        cansig75downcounts=$(wc -l ${cansig75down} | awk '{print $1}')
        cansig90up="${cancer}_gts_up_sig90.tsv"
        cansig90upcounts=$(wc -l ${cansig90up} | awk '{print $1}')
        cansig90down="${cancer}_gts_down_sig90.tsv"
        cansig90downcounts=$(wc -l ${cansig90down} | awk '{print $1}')


        echo -e "${cancer} \t${allsigcounts} \t${cangtcounts} \t${allupcounts} \t${can05upgtcounts} \t${alldowncounts} \t${can05downgtcounts} \t${all01counts} \t${can01gtcounts} \t${all01upcounts} \t${can01upgtcounts} \t${all01downcounts} \t${can01downgtcounts} \t${allfcupcounts} \t${canfcupcounts} \t${allfcdowncounts} \t${canfcdowncounts} \t${allsig75upcounts} \t${cansig75upcounts} \t${allsig75downcounts} \t${cansig75downcounts} \t${allsig90upcounts} \t${cansig90upcounts} \t${allsig90downcounts} \t${cansig90downcounts}" >> ${path_generated}/total_counts.txt

        #Print master significant lists
        echo -e "UpGT05Genes \tUpGT01Genes \tUpGTFC1Genes \tUpGTSig75Genes \tUPGTSig90Genes \tDownGT05Genes \tDownGT01Genes \tDownGTFC1Genes \tDownGTSig75Genes \tDownGTSig90Genes" > ${path_generated}/${cancer}_sig_lists.txt
        paste ${cancer}_05_up_list.tsv ${cancer}_01_up_list.tsv ${cancer}_gts_up_fc_list.tsv ${cancer}_gts_up_sig75_list.tsv ${cancer}_gts_up_sig90_list.tsv ${cancer}_05_down_list.tsv ${cancer}_01_down_list.tsv ${cancer}_gts_down_fc_list.tsv ${cancer}_gts_down_sig75_list.tsv ${cancer}_gts_down_sig90_list.tsv | awk -F'\t' '{printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' >> ${path_generated}/${cancer}_sig_lists.txt

	#Print master gts in bioxpress
	cat ${cancer}_gts.tsv >> ${path_generated}/Differential_expression_gts_master.tsv
	
done
