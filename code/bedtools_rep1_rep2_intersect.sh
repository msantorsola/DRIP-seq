#!/bin/bash

#SBATCH --time 20:00:00
#SBATCH --job-name=reg_inter_replicates 
#SBATCH -p knl_usr_prod
#SBATCH -N 1
#SBATCH -n 20 
#SBATCH --mem=86000
#SBATCH --output=reg_inter_replicates.%J.out
#SBATCH --error=reg_inter_replicates.%J.err
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
idrDir=$WORK/drip_project/qual/IDR_analysis


cd ${idrDir}

exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1>${logDir}/intersect_replicates.log 2>&1



bedtools intersect -wao -a ${peakDir}/regular_Model_rep1_sorted_peaks.narrowPeak -b ${peakDir}/regular_Model_rep2_sorted_peaks.narrowPeak > ${peakDir}/regular_Model_intersect.bed
bedtools intersect -wao -a ${peakDir}/regular_noModel_rep1_sorted_peaks.narrowPeak -b ${peakDir}/regular_noModel_rep2_sorted_peaks.narrowPeak > ${peakDir}/regular_noModel_intersect.bed

bedtools intersect -wa -a ${peakDir}/regular_Model_rep1_sorted_peaks.narrowPeak -b ${peakDir}/regular_Model_rep2_sorted_peaks.narrowPeak > ${peakDir}/regular_Model_intersect_rep1.bed
bedtools intersect -wb -a ${peakDir}/regular_Model_rep1_sorted_peaks.narrowPeak -b ${peakDir}/regular_Model_rep2_sorted_peaks.narrowPeak > ${peakDir}/regular_Model_intersect_rep2.bed


bedtools intersect -wao -a ${peakDir}/broad_Model_rep1_sorted_peaks.broadPeak -b ${peakDir}/broad_Model_rep2_sorted_peaks.broadPeak > ${peakDir}/broad_Model_intersect.bed
bedtools intersect -wao -a ${peakDir}/broad_noModel_rep1_sorted_peaks.broadPeak -b ${peakDir}/broad_noModel_rep2_sorted_peaks.broadPeak > ${peakDir}/broad_noModel_intersect.bed


bedtools intersect -wa -a ${peakDir}/broad_Model_rep1_sorted_peaks.broadPeak -b ${peakDir}/broad_Model_rep2_sorted_peaks.broadPeak > ${peakDir}/broad_Model_intersect_rep1.bed
bedtools intersect -wb -a ${peakDir}/broad_noModel_rep1_sorted_peaks.broadPeak -b ${peakDir}/broad_noModel_rep2_sorted_peaks.broadPeak > ${peakDir}/broad_noModel_intersect_rep2.bed

#Irreproducible Discovery Rate (IDR)


