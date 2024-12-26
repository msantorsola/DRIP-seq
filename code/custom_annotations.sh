#!/bin/bash

#SBATCH --time 20:00:00
#SBATCH --job-name=custom_anno  
#SBATCH -p knl_usr_prod
#SBATCH -N 1
#SBATCH -n 20 
#SBATCH --mem=86000
#SBATCH --output=custom_anno .%J.out
#SBATCH --error=custom_anno .%J.err
#SBATCH --account=<ACCOUNT_NAME>

module load profile/global
module load profile/bioinf

module load autoload bedtools2

workDir=$WORK/drip_project/
tmpDir=$WORK/drip_project/tmp
aligDir=$WORK/drip_project/alignment
logDir=$WORK/drip_project/logs
qualDir=$WORK/drip_project/qual
peakDir=$WORK/drip_project/peak_calling
rawDir=$WORK/drip_raw
data=$WORK/data_annotation

#create a logfile 
exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1>${logDir}/custom_annotations.log 2>&1

cd ${workDir}


#bedtools window -a A.bed -b B.bed -w 10000

#A.bed
#$WORK/drip_project/qual/IDR_analysis/regular_noModel-idr
#$WORK/drip_project/qual/IDR_analysis/broad_noModel-idr

#$WORK/drip_project/peak_calling/regular_noModel_intersect_common.bed
#$WORK/drip_project/peak_calling/broad_noModel_intersect_common.bed

#bedtools window -a $WORK/drip_project/qual/IDR_analysis/regular_noModel-idr -b -w 10000

#regular IDR
bedtools window -a $WORK/drip_project/qual/IDR_analysis/regular_noModel-idr -b $WORK/data_annotation/promoter_enhancer/wgEncodeSydhHistoneU2osH3k9me3UcdPk.narrowPeak -w 10000 > $WORK/drip_project/peak_calling/annotations/regular_noModel-idr_vs_U2osH3k9me3.bed 
bedtools window -a $WORK/drip_project/qual/IDR_analysis/regular_noModel-idr -b $WORK/data_annotation/promoter_enhancer/wgEncodeSydhHistoneU2osH3k36me3bUcdPk.narrowPeak -w 10000 > $WORK/drip_project/peak_calling/annotations/regular_noModel-idr_vs_U2osH3k36me3b.bed 
bedtools window -a $WORK/drip_project/qual/IDR_analysis/regular_noModel-idr -b $WORK/data_annotation/promoter_enhancer/wgEncodeSydhHistoneMcf7H3k09me3UcdPk.narrowPeak -w 10000 > $WORK/drip_project/peak_calling/annotations/regular_noModel-idr_vs_Mcf7H3k09me3.bed 
bedtools window -a $WORK/drip_project/qual/IDR_analysis/regular_noModel-idr -b $WORK/data_annotation/promoter_enhancer/wgEncodeSydhHistoneMcf7H3k27acUcdPk.narrowPeak -w 10000 > $WORK/drip_project/peak_calling/annotations/regular_noModel-idr_vs_Mcf7H3k27ac.bed 
bedtools window -a $WORK/drip_project/qual/IDR_analysis/regular_noModel-idr -b $WORK/data_annotation/promoter_enhancer/wgEncodeSydhHistoneMcf7H3k27me3bUcdPk.narrowPeak -w 10000 > $WORK/drip_project/peak_calling/annotations/regular_noModel-idr_vs_Mcf7H3k27me3.bed 
bedtools window -a $WORK/drip_project/qual/IDR_analysis/regular_noModel-idr -b $WORK/data_annotation/promoter_enhancer/wgEncodeSydhHistoneMcf7H3k36me3bUcdPk.narrowPeak -w 10000 > $WORK/drip_project/peak_calling/annotations/regular_noModel-idr_vs_Mcf7H3k36me3.bed 

for anno in $(ls $WORK/data_annotation/g_quadruplex_forming_repeats/*_GQ.tsv)
do
	bedtools window -a $WORK/drip_project/qual/IDR_analysis/regular_noModel-idr -b ${anno} > $WORK/drip_project/peak_calling/annotations/regular_noModel-idr_vs_U2osH3k9me3.bed

#cpgIsland
bedtools window -a $WORK/drip_project/qual/IDR_analysis/regular_noModel-idr -b $WORK/data_annotation/CpGIsland/cpgIslandExt_1.txt -w 10000 > $WORK/drip_project/peak_calling/annotations/regular_noModel-idr_vs_cpgIslandExt.bed

#promoter/enhancer
$WORK/data_annotation/promoter_enhancer/wgEncodeSydhHistoneU2osH3k9me3UcdPk.narrowPeak
$WORK/data_annotation/promoter_enhancer/wgEncodeSydhHistoneU2osH3k36me3bUcdPk.narrowPeak
$WORK/data_annotation/promoter_enhancer/wgEncodeSydhHistoneMcf7H3k09me3UcdPk.narrowPeak
$WORK/data_annotation/promoter_enhancer/wgEncodeSydhHistoneMcf7H3k27acUcdPk.narrowPeak
$WORK/data_annotation/promoter_enhancer/wgEncodeSydhHistoneMcf7H3k27me3bUcdPk.narrowPeak
$WORK/data_annotation/promoter_enhancer/wgEncodeSydhHistoneMcf7H3k36me3bUcdPk.narrowPeak

#nonB dna
bedtools window -a $WORK/drip_project/qual/IDR_analysis/regular_noModel-idr -b $WORK/data_annotation/G4-nonB_DB/GSE63874_Na_K_PDS_minus_hits_intersect.bed > $WORK/drip_project/peak_calling/annotations/regular_noModel-idr_vs_GSE63874_Na_K_PDS_minus.bed
bedtools window -a $WORK/drip_project/qual/IDR_analysis/regular_noModel-idr -b $WORK/data_annotation/G4-nonB_DB/GSE63874_Na_K_PDS_plus_hits_intersect.bed > $WORK/drip_project/peak_calling/annotations/regular_noModel-idr_vs_GSE63874_Na_K_PDS_plus.bed
bedtools window -a $WORK/drip_project/qual/IDR_analysis/regular_noModel-idr -b $WORK/data_annotation/G4-nonB_DB/GSE63874_Na_K_plus_hits_intersect.bed > $WORK/drip_project/peak_calling/annotations/regular_noModel-idr_vs_GSE63874_Na_K_plus.bed
bedtools window -a $WORK/drip_project/qual/IDR_analysis/regular_noModel-idr -b $WORK/data_annotation/G4-nonB_DB/GSE63874_Na_PDS_minus_hits_intersect.bed > $WORK/drip_project/peak_calling/annotations/regular_noModel-idr_vs_GSE63874_Na_PDS_minus.bed
bedtools window -a $WORK/drip_project/qual/IDR_analysis/regular_noModel-idr -b $WORK/data_annotation/G4-nonB_DB/GSE63874_Na_PDS_plus_hits_intersect.bed > $WORK/drip_project/peak_calling/annotations/regular_noModel-idr_vs_GSE63874_PDS_plus_hits.bed
bedtools window -a $WORK/drip_project/qual/IDR_analysis/regular_noModel-idr -b $WORK/data_annotation/G4-nonB_DB/GSE63874_Na_K_minus_hits_intersect.bed > $WORK/drip_project/peak_calling/annotations/regular_noModel-idr_vs_GSE63874_Na_K_minus.bed


