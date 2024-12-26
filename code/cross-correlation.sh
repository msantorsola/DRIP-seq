#!/bin/bash

#SBATCH --time 20:00:00
#SBATCH --job-name=stats_cross  
#SBATCH -p knl_usr_prod
#SBATCH -N 1
#SBATCH -n 20 
#SBATCH --mem=86000
#SBATCH --output=stats_cross.%J.out
#SBATCH --error=stats_cross.%J.err
#SBATCH --account=<ACCOUNT_NAME>


module load profile/global
module load profile/bioinf

module load autoload r
module load autoload python


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
exec 1>${logDir}/cross-correlation.log 2>&1


#Cross-correlation with SPP
echo "Cross-correlation analysis..."

for bam in $(ls ${aligDir}/*_uniq.bam)

	do 
	
	bam2=`basename $bam _uniq.bam`
	
	Rscript $WORK/tools/spp/run_spp.R -c=${bam} -savp -out=${qualDir}/${bam2}.qual -tmpdir="${tmpDir}" -fdr=0.05 -rf -x=-500:85 2> ${logDir}/${bam2}.Rout

done


echo "done"
#The for loop generates three output files. 
#The quality metrics are written in a tab-delimited text file, and the log files contains the standard output text. 
#A third file is created in the same directory as the BAM files. 
#These are pdf files that contain the cross-correlation plot for each sample.

#strand cross-correlation analysis
#CRTL
#time Rscript $WORK/tools/spp/run_spp.R -c=OVCAR_CTRL_S2_uniq.bam -savp -out=$WORK/drip_project_test/align_stats/CRTL-cross-correlation.txt -tmpdir="${tmpDir}" -fdr=0.05 -rf -x=-500:85
#time Rscript $WORK/tools/spp/run_spp.R -c=OVCAR_P3_S1_uniq.bam -savp -out=$WORK/drip_project_test/align_stats/EXP-cross-correlation.txt -tmpdir="${tmpDir}" -fdr=0.05 -rf -x=-500:85
