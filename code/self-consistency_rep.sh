#!/bin/bash

#SBATCH --time 20:00:00
#SBATCH --job-name=self_rep1_reg
#SBATCH -p knl_usr_prod
#SBATCH -N 1
#SBATCH -n 20 
#SBATCH --mem=86000
#SBATCH --output=self_rep1_reg.%J.out
#SBATCH --error=self_rep1_reg.%J.err
#SBATCH --account=<ACCOUNT_NAME>


module load profile/global
module load profile/bioinf

module load autoload r
module load autoload python/3.6.4
module load autoload macs
module load autoload samtools


#create a logfile 
exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1>self_rep1_reg.log 2>&1


MACS(){
    /cineca/prod/opt/applications/macs/2.1.0/python--2.7.12/bin/macs2 "$@"
}
export MACS

date 


#file name
ctrl="siCon_S9_6_1_RV_S1_uniq.bam" 
rep1="siPSFM_S9_6_2_RV_S1_uniq.bam" 
rep2="siPSFD_S9_6_3_RV_S2_uniq.bam"



# Make Directories
mkdir -p $WORK/drip_project_test/IDR_analysis/self-consistency
chmod -R 755 $WORK/drip_project_test/IDR_analysis/self-consistency

mkdir -p $WORK/drip_project_test/IDR_analysis/self-consistency/macs
chmod -R 755 $WORK/drip_project_test/IDR_analysis/self-consistency/macs

# Set paths
macsDir=$WORK/drip_project_test/IDR_analysis/self-consistency/macs
outputDir=$WORK/drip_project_test/IDR_analysis/pooled_pseudoreps

workDir=$WORK/drip_project/
tmpDir=$WORK/drip_project/tmp
aligDir=$WORK/drip_project/alignment
logDir=$WORK/drip_project/logs
qualDir=$WORK/drip_project/qual
peakDir=$WORK/drip_project/peak_calling
rawDir=$WORK/drip_raw
idrDir=$WORK/drip_project/qual/IDR_analysis



# Number of reads in the rep1 BAM file
#nlines=$(samtools view ${baseDir}/${rep1} | wc -l ) 
# half that number
#nlines=$(( (nlines + 1) / 2 ))
nlines="24958394"
lk${tmpDir}/rep1.00.bam
cat ${tmpDir}/rep1_header.sam ${tmpDir}/rep101 | samtools view -bS - > ${tmpDir}/rep1.01.bam


#Peak calling on pseudoreplicates
echo "Peak calling on regular rep1 replicates..."

echo "Calling peaks for regular/model rep1 replicate1"
macs2 callpeak -t ${tmpDir}/rep1.00.bam -c ${aligDir}/${ctrl} --format BAMPE -g hs --keep-dup=auto --bdg -q 0.05 --call-summits -n rep1_pr1_reg_Mod --outdir $tmpDir

echo "Calling peaks for regular/model rep1 replicate2"
macs2 callpeak -t ${tmpDir}/rep1.01.bam -c ${aligDir}/${ctrl} --format BAMPE -g hs --keep-dup=auto --bdg -q 0.05 --call-summits -n rep1_pr2_reg_Mod --outdir $tmpDir

echo "Calling peaks for regular/no model rep1 replicate1"
macs2 callpeak -t ${tmpDir}/rep1.00.bam -c ${aligDir}/${ctrl} --format BAMPE -g hs --keep-dup=auto --bdg --nomodel --extsize 200 -q 0.05 --call-summits -n rep1_pr1_reg_noMod --outdir $tmpDir

echo "Calling peaks for regular/no model rep1 replicate2"
macs2 callpeak -t ${tmpDir}/${rep1}.01.bam -c ${aligDir}/${ctrl} --format BAMPE -g hs --keep-dup=auto --bdg --nomodel --extsize 200 -q 0.05 --call-summits -n rep1_pr2_reg_noMod --outdir $tmpDir


echo "Peak calling on broad rep1 replicates..."

echo "Calling peaks for broad/noModel rep1 replicate1"
macs2 callpeak -t ${tmpDir}/${rep1}.00.bam -c ${aligDir}/${ctrl} --format BAMPE --broad -g hs --broad-cutoff 0.1 --keep-dup=auto --bdg --nomodel --extsize 200 -q 0.05 -n rep1_pr1_broad_noMod --outdir $tmpDir

echo "Calling peaks for broad/noModel rep1 replicate2"
macs2 callpeak -t ${tmpDir}/${rep1}.01.bam -c ${aligDir}/${ctrl} --format BAMPE --broad -g hs --broad-cutoff 0.1 --keep-dup=auto --bdg --nomodel --extsize 200 -q 0.05 -n rep1_pr2_broad_noMod --outdir $tmpDir

echo "Calling peaks for broad/model rep1 replicate1"
macs2 callpeak -t ${tmpDir}/${rep1}.00.bam -c ${aligDir}/${ctrl} --format BAMPE --broad -g hs --broad-cutoff 0.1 --keep-dup=auto --bdg -q 0.05 rep1_pr1_broad_Mod --outdir $tmpDir

echo "Calling peaks for broad/model rep1 replicate2"
macs2 callpeak -t ${tmpDir}/${rep1}.01.bam -c ${aligDir}/${ctrl} --format BAMPE --broad -g hs --broad-cutoff 0.1 --keep-dup=auto --bdg -q 0.05rep1_pr2_broad_Mod --outdir $tmpDir


#Sort peak by -log10(p-value)
echo "Sorting peaks..."
sort -k8,8nr ${tmpDir}/rep1_pr1_reg_noMod_peaks.narrowPeak | head -n 100000 > ${peakDir}/rep1_pr1_reg_noMod_sorted.narrowPeak
sort -k8,8nr ${tmpDir}/rep1_pr2_reg_noMod_peaks.narrowPeak | head -n 100000 > ${peakDir}/rep1_pr2_reg_noMod_sorted.narrowPeak
sort -k8,8nr ${tmpDir}/rep1_pr1_reg_Mod_peaks.narrowPeak | head -n 100000 > ${peakDir}/rep1_pr1_reg_Mod_sorted.narrowPeak
sort -k8,8nr ${tmpDir}/rep1_pr2_reg_Mod_peaks.narrowPeak | head -n 100000 > ${peakDir}/rep1_pr2_reg_Mod_sorted.narrowPeak
sort -k8,8nr ${tmpDir}/rep1_pr1_broad_Mod_peaks.broadPeak | head -n 100000 > ${peakDir}/rep1_pr1_broad_Mod_sorted.broadPeak
sort -k8,8nr ${tmpDir}/rep1_pr2_broad_Mod_peaks.broadPeak | head -n 100000 > ${peakDir}/rep1_pr2_broad_Mod_sorted.broadPeak
sort -k8,8nr ${tmpDir}/rep1_pr1_broad_noMod_peaks.broadPeak | head -n 100000 > ${peakDir}/rep1_pr1_broad_noMod_sorted.broadPeak
sort -k8,8nr ${tmpDir}/rep1_pr2_broad_noMod_peaks.broadPeak | head -n 100000 > ${peakDir}/rep1_pr2_broad_noMod_sorted.broadPeak


#Independent replicate IDR
echo "Running IDR on rep1 replicates..."

echo "Running IDR on regular/noModel rep1 replicates..."
idr --samples ${peakDir}/rep1_pr1_reg_noMod_sorted.narrowPeak ${peakDir}/rep1_pr2_reg_noMod_sorted.narrowPeak \
--input-file-type narrowPeak \
--output-file rep1_pr_rnM-idr \
--rank p.value \
--plot \
--log-output-file rep1_pr_rnM.idr.log


echo "Running IDR on regular/model rep1 replicates..."
idr --samples ${peakDir}/rep1_pr1_reg_Mod_sorted.narrowPeak ${peakDir}/rep1_pr2_reg_Mod_sorted.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file rep1_pr_rm-idr \
--plot \
--log-output-file rep1_pr_rm.idr.log



echo "Running IDR on broad/noModel rep1 replicates..."
idr --samples ${peakDir}/rep1_pr1_broad_noMod_sorted.broadPeak ${peakDir}/rep1_pr2_broad_noMod_sorted.broadPeak \
--input-file-type narrowPeak \
--output-file rep1_pr_bnM-idr \
--rank p.value \
--plot \
--log-output-file rep1_pr_bnM.idr.log

echo "Running IDR on broad/model rep1 replicates..."
idr --samples ${peakDir}/rep1_pr1_broad_Mod_sorted.broadPeak ${peakDir}/rep1_pr2_broad_Mod_sorted.broadPeak \
--input-file-type narrowPeak \
--output-file rep1_pr_bm-idr \
--rank p.value \
--plot \
--log-output-file rep1_pr_bm.idr.log


