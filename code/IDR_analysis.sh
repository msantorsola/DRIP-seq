#!/bin/bash

#SBATCH --time 20:00:00
#SBATCH --job-name=IDR_analysis  
#SBATCH -p knl_usr_prod
#SBATCH -N 1
#SBATCH -n 20 
#SBATCH --mem=86000
#SBATCH --output=IDR_analysis.%J.out
#SBATCH --error=IDR_analysis.%J.err
#SBATCH --account=<ACCOUNT_NAME>


module load profile/global
module load profile/bioinf

module load autoload python/3.6.4


#create a logfile 
exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1>IDR_analysis.log 2>&1

workDir=$WORK/drip_project/
tmpDir=$WORK/drip_project/tmp
aligDir=$WORK/drip_project/alignment
logDir=$WORK/drip_project/logs
qualDir=$WORK/drip_project/qual
peakDir=$WORK/drip_project/peak_calling
rawDir=$WORK/drip_raw

cd ${workDir}


mkdir -p $WORK/drip_project/qual/IDR_analysis
chmod -R 755 $WORK/drip_project/qual/IDR_analysis

mkdir -p $WORK/drip_project/qual/IDR_analysis/true_rep
chmod -R 755 $WORK/drip_project_test/IDR_analysis/true_rep


cd $WORK/drip_project_test/IDR_analysis/true_rep


#Sort peak by -log10(p-value)

#regular/model
echo "Sorting MACS regular peaks with model"
sort -k8,8nr $WORK/drip_project_test/MACS_peaks/regular_peaks_Model_rep1_peaks.narrowPeak > $WORK/drip_project_test/MACS_peaks/regular_Model_rep1_sorted_peaks.narrowPeak
sort -k8,8nr $WORK/drip_project_test/MACS_peaks/regular_peaks_Model_rep2_peaks.narrowPeak > $WORK/drip_project_test/MACS_peaks/regular_Model_rep2_sorted_peaks.narrowPeak

#regular/NoModel
echo "Sorting MACS regular peaks without model"
sort -k8,8nr $WORK/drip_project_test/MACS_peaks/regular_peaks_noModel_rep1_peaks.narrowPeak > $WORK/drip_project_test/MACS_peaks/regular_noModel_rep1_sorted_peaks.narrowPeak
sort -k8,8nr $WORK/drip_project_test/MACS_peaks/regular_peaks_noModel_rep2_peaks.narrowPeak > $WORK/drip_project_test/MACS_peaks/regular_noModel_rep2_sorted_peaks.narrowPeak

#broad/model
echo "Sorting MACS broad peaks with model"
sort -k8,8nr $WORK/drip_project_test/MACS_peaks/broad_peaks_Model_rep1_peaks.narrowPeak > $WORK/drip_project_test/MACS_peaks/broad_Model_rep1_sorted_peaks.narrowPeak
sort -k8,8nr $WORK/drip_project_test/MACS_peaks/broad_peaks_Model_rep2_peaks.narrowPeak > $WORK/drip_project_test/MACS_peaks/broad_Model_rep2_sorted_peaks.narrowPeak

#broad/noModel
echo "Sorting MACS broad peaks without model"
sort -k8,8nr $WORK/drip_project_test/MACS_peaks/broad_peaks_noModel_rep1_peaks.narrowPeak > $WORK/drip_project_test/MACS_peaks/broad_noModel_rep1_sorted_peaks.narrowPeak
sort -k8,8nr $WORK/drip_project_test/MACS_peaks/broad_peaks_noModel_rep2_peaks.narrowPeak > $WORK/drip_project_test/MACS_peaks/broad_noModel_rep2_sorted_peaks.narrowPeak



#Peak consistency between true replicates

#regular/model 
echo "Running IDR on regular-model replicates..."
idr --samples $WORK/drip_project_test/MACS_peaks/regular_Model_rep1_sorted_peaks.narrowPeak $WORK/drip_project_test/MACS_peaks/regular_Model_rep2_sorted_peaks.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file regular_model-idr \
--plot \
--log-output-file regular_model.idr.log

#regular/noModel 
echo "Running IDR on regular-noModel replicates..."
idr --samples $WORK/drip_project_test/MACS_peaks/regular_noModel_rep1_sorted_peaks.narrowPeak $WORK/drip_project_test/MACS_peaks/regular_noModel_rep2_sorted_peaks.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file regular_noModel-idr \
--plot \
--log-output-file regular_noModel.idr.log


#broad/model 
echo "Running IDR on broad-model replicates..."
idr --samples $WORK/drip_project_test/MACS_peaks/broad_Model_rep1_sorted_peaks.narrowPeak $WORK/drip_project_test/MACS_peaks/regular_Model_rep2_sorted_peaks.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file broad_model-idr \
--plot \
--log-output-file broad_model.idr.log

#broad/noModel 
echo "Running IDR on broad-noModel replicates..."
idr --samples $WORK/drip_project_test/MACS_peaks/broad_noModel_rep1_sorted_peaks.narrowPeak $WORK/drip_project_test/MACS_peaks/regular_noModel_rep2_sorted_peaks.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file broad_noModel-idr \
--plot \
--log-output-file broad_noModel.idr.log


#How many common peaks are considered for each TF?
echo "Count total common peaks..."
wc -l *-idr > common_peak_number.txt

#To find out how may of those shared regions have an IDR < 0.05, we can take a look at the log files. 
#Alternatively, since we requested all peaks and their IDR value as output we can also filter the file using an awk command.
echo "Count common peaks with IDR<0.05..."
awk '{if($12 > 1.3) print $0}' regular_model-idr | wc -l > regular_model-idr_0.05.txt
awk '{if($12 > 1.3) print $0}' regular_noModel-idr | wc -l > regular_noModel-idr_0.05.txt
awk '{if($12 > 1.3) print $0}' broad_model-idr | wc -l > broad_Model-idr_0.05.txt
awk '{if($12 > 1.3) print $0}' broad_noModel-idr | wc -l > broad_noModel-idr_0.05.txt
