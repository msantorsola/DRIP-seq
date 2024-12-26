s#!/bin/bash

#SBATCH --time 20:00:00
#SBATCH --job-name=peak_calling 
#SBATCH -p knl_usr_prod
#SBATCH -N 1
#SBATCH -n 20 
#SBATCH --mem=86000
#SBATCH --output=peak_calling.%J.out
#SBATCH --error=peak_calling.%J.err
#SBATCH --account=<ACCOUNT_NAME>


module load profile/global
module load profile/bioinf

module load autoload ig_homo_sapiens/hg19
module load autoload r
module load autoload macs
module load autoload bedtools2


Rscript(){
	$R_HOME/bin/R/Rscript 
}
export Rscript


#GEM(){
#	java -jar -Xmx40G $WORK/gem/gem.jar --d $WORK/gem/Read_Distribution_default.txt --g $WORK/gem/hg19.chrom.sizes --genome $IG_HG19_GENOME "$@"
#}
#export GEM


MACS(){
	/cineca/prod/opt/applications/macs/2.1.0/python--2.7.12/bin/macs2 "$@"
}
export MACS


bedtools(){
	$BEDTOOLS2_HOME/bin/bedtools "$@"
}
export bedtools


#create a logfile 
exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1>peak_calling_MACS2.log 2>&1



workDir=$WORK/drip_project/
tmpDir=$WORK/drip_project/tmp
aligDir=$WORK/drip_project/alignment
logDir=$WORK/drip_project/logs
qualDir=$WORK/drip_project/qual
peakDir=$WORK/drip_project/peak_calling
rawDir=$WORK/drip_raw

cd ${workDir}


sampleName="U2OS"

#MACS2 calling

echo "Peak calling with MACS2..."

echo "Calling regular sharp peaks without model..."
#call regular sharp peaks without model
macs2 callpeak -t ${aligDir}/${sampleName}_rep1_uniq.bam -c ${aligDir}/${sampleName}_ctrl_uniq.bam --format BAMPE  --name "regular_noModel_rep1" --outdir "${peakDir}" -g hs --keep-dup=auto --bdg --nomodel --extsize 200 -q 0.05 --call-summits  
macs2 callpeak -t ${aligDir}/${sampleName}_rep2_uniq.bam -c ${aligDir}/${sampleName}_ctrl_uniq.bam --format BAMPE  --name "regular_noModel_rep2" --outdir "${peakDir}" -g hs --keep-dup=auto --bdg --nomodel --extsize 200 -q 0.05 --call-summits  

echo "done"

echo "Calling regular sharp peaks with model..."
#call regular sharp peaks with model
macs2 callpeak -t ${aligDir}/${sampleName}_rep1_uniq.bam -c ${aligDir}/${sampleName}_ctrl_uniq.bam --format BAMPE  --name "regular_Model_rep1" --outdir "${peakDir}" -g hs --keep-dup=auto --bdg -q 0.05 --call-summits 
macs2 callpeak -t ${aligDir}/${sampleName}_rep2_uniq.bam -c ${aligDir}/${sampleName}_ctrl_uniq.bam --format BAMPE  --name "regular_Model_rep2" --outdir "${peakDir}" -g hs --keep-dup=auto --bdg -q 0.05 --call-summits 
echo "done"

echo "Calling broad peaks without model..."
#call broad peaks without model
macs2 callpeak -t ${aligDir}/${sampleName}_rep1_uniq.bam -c ${aligDir}/${sampleName}_ctrl_uniq.bam --format BAMPE  --name "broad_noModel_rep1" --outdir "${peakDir}" --broad -g hs --broad-cutoff 0.1 --keep-dup=auto --bdg --nomodel --extsize 200 -q 0.05 
macs2 callpeak -t ${aligDir}/${sampleName}_rep2_uniq.bam -c ${aligDir}/${sampleName}_ctrl_uniq.bam --format BAMPE  --name "broad_noModel_rep2" --outdir "${peakDir}" --broad -g hs --broad-cutoff 0.1 --keep-dup=auto --bdg --nomodel --extsize 200 -q 0.05 
echo "done"

echo "Calling broad peaks with model..."
#call broad peaks with model
macs2 callpeak -t ${aligDir}/${sampleName}_rep1_uniq.bam -c ${aligDir}/${sampleName}_ctrl_uniq.bam --format BAMPE  --name "broad_Model_rep1" --outdir "${peakDir}" --broad -g hs --broad-cutoff 0.1 --keep-dup=auto --bdg -q 0.05 
macs2 callpeak -t ${aligDir}/${sampleName}_rep2_uniq.bam -c ${aligDir}/${sampleName}_ctrl_uniq.bam --format BAMPE  --name "broad_Model_rep2" --outdir "${peakDir}" --broad -g hs --broad-cutoff 0.1 --keep-dup=auto --bdg -q 0.05 
echo "done"


#GEM
#peak calling with gem
#java -jar -Xmx40G $WORK/gem/gem.jar --d $WORK/gem/Read_Distribution_default.txt --g $WORK/gem/hg19.chrom.sizes --genome $IG_HG19_GENOME --expt OVCAR_P3_S1_uniq.bed --ctrl OVCAR_CTRL_S2_uniq.bed --f BED --out $WORK/drip_project_test/GEM_peaks --k_min 6 --k_max 20 --k_seqs 1000 --outNP --outMEME --outBED --k_neg_dinu_shuffle --t 20 

#SPP
#Rscript $WORK/run_spp.R -c=OVCAR_CTRL_S2_uniq.bam -i=OVCAR_P3_S1_uniq.bam -fdr=0.05 -odir=$WORK/drip_project_test/SPP_peaks -savr -savp -savd -rf


### call regular sharp peaks with model
#cat sample_names.txt | parallel --max-procs=12 'macs2 callpeak -t {}-A-NC.sorted.bam \
# -c {}-G-NC.sorted.bam -g hs -n {}-A-NC-sharp-model -q 0.01 --outdir {}-A-NC-sharp-model-peaks 2> {}-A-NC-sharp-model.stderr'

### call regular sharp peaks without model
#cat sample_names.txt | parallel --max-procs=12 'macs2 callpeak -t {}-A-NC.sorted.bam \
#-c {}-G-NC.sorted.bam -g hs -n {}-A-NC-sharp-nomodel -q 0.01 \
# --nomodel --extsize 146 --outdir {}-A-NC-sharp-nomodel-peaks 2> {}-A-NC-sharp-nomodel.stderr'

### call broad peaks with model
#cat sample_names.txt | parallel --max-procs=12 'macs2 callpeak -t {}-A-NC.sorted.bam \
# -c {}-G-NC.sorted.bam --broad -g hs --broad-cutoff 0.1 -n {}-A-NC-broad-model -q 0.01 \
# --outdir {}-A-NC-broad-model-peaks 2> {}-A-NC-broad-model.stderr'

### call broad peaks without model
#cat sample_names.txt | parallel --max-procs=12 'macs2 callpeak -t {}-A-NC.sorted.bam \
# -c {}-G-NC.sorted.bam --broad -g hs --broad-cutoff 0.1 -n {}-A-NC-broad-nomodel -q 0.01 \
# --nomodel --extsize 146 --outdir {}-A-NC-broad-nomodel-peaks 2> {}-A-NC-broad-nomodel.stderr'


