#!/bin/bash

#SBATCH --time 20:00:00
#SBATCH --job-name=align_stats  
#SBATCH -p knl_usr_prod
#SBATCH -N 1
#SBATCH -n 20 
#SBATCH --mem=86000
#SBATCH --output=align_stats.%J.out
#SBATCH --error=align_stats.%J.err
#SBATCH --account=<ACCOUNT_NAME>


module load profile/global
module load profile/bioinf

module load autoload samtools
module load autoload r



Rscript(){
	$R_HOME/bin/R/Rscript 
}
export Rscript


samtools(){
     $SAMTOOLS_HOME/bin/samtools  "$@"
}
export samtools


workDir=$WORK/drip_project/
tmpDir=$WORK/drip_project/tmp
aligDir=$WORK/drip_project/alignment
logDir=$WORK/drip_project/logs
qualDir=$WORK/drip_project/qual
peakDir=$WORK/drip_project/peak_calling
rawDir=$WORK/drip_raw

cd ${workDir}


#create a logfile 
exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1>${logDir}/alignment_stats.log 2>&1


#declare *_uniq.bam files

for bam in $(ls ${aligDir}/*_uniq.bam)

do

	#for i in *_R1.fastq.gz; do   SAMPLE=$(echo ${i} | sed "s/_R1\.fastq\.gz//")
	echo "Running flagstat..."
	samtools flagstat ${aligDir}/${bam} > ${qualDir}/$(echo ${bam} | sed "s/_uniq\.bam//").flagstats

	#samtools flagstat $WORK/drip_project_test/OVCAR_P3_S1_uniq.bam > OVCAR_P3_S1.flagstats

	echo "Counting read number in the BAM files..."
	# Number of reads in the BAM file
	samtools view ${aligDir}/${bam} | wc -l > ${qualDir}/$(echo ${bam} | sed "s/_uniq\.bam//").read_num.txt
	#samtools view $WORK/drip_project_test/OVCAR_P3_S1_uniq.bam | wc -l > exp_read_num.txt

done

echo "Removing tmp files..."
# Remove the tmp directory
rm -r ${tmpDir}/*

echo "done"

