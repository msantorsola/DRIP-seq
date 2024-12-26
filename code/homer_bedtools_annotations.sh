#!/bin/bash

#SBATCH --time 20:00:00
#SBATCH --job-name=homer_bedtools
#SBATCH -p knl_usr_prod
#SBATCH -N 1
#SBATCH -n 20 
#SBATCH --mem=86000
#SBATCH --output=homer_bedtools.%J.out
#SBATCH --error=homer_bedtools.%J.err
#SBATCH --account=<ACCOUNT_NAME>


module load profile/global
module load profile/bioinf


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
exec 1>${logDir}/homer_bedtools_annotations.log 2>&1


for bed in $(ls ${peakDir}/*_intersect.bed)

do

        bed_name=`basename $bed .bed`

        annotatePeaks.pl ${bed} hg19 -annStats ${qualDir}/${bed_name}.stats  > ${peakDir}/${bed_name}_homer_annotations.txt

done


