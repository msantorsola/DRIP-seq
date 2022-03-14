#!/bin/bash


#SBATCH --time 20:00:00
#SBATCH --job-name=231_peaks
#SBATCH --partition=g100_usr_prod
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --mem=86000
#SBATCH --output=231_peaks.%J.out
#SBATCH --error=231_peaks.%J.err
#SBATCH --account=***

#SBATCH --mail-user=***@gmail.com
#SBATCH --mail-type=ALL,TIME_LIMIT_90



module load profile/bioinf
module load autoload fastqc/0.11.5
module load autoload trimmomatic
module load autoload ig_UCSC_Homo_sapiens/hg19
module load autoload bwa
module load autoload picardtools/2.3.0
module load autoload macs/2.1.2

workDir=***/drip_2021_231
tmpDir=***/drip_2021_231/tmp
aligDir=***/drip_2021_231/alignment
logDir=***/drip_2021_231/logs
qualDir=***/drip_2021_231/qual
varDir=***/drip_2021_231/peaks



#../alignment/Sample_2_S2_uniq.bam
#../alignment/Sample_1_S1_uniq.bam
#../alignment/Sample_4_S4_uniq.bam
#../alignment/Sample_3_S3_uniq.bam

#create a logfile 
exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1>231_peaks.log 2>&1

#Peak calling, MACS2
#narrow
$MACS_HOME/bin/macs2 callpeak -t ${aligDir}/Sample_2_S2_uniq.bam -c ${aligDir}/Sample_1_S1_uniq.bam -f BAMPE --keep-dup auto -g hs --outdir ${varDir} -q 1 --name Sample_2_S2_narrow 
$MACS_HOME/bin/macs2 callpeak -t ${aligDir}/Sample_3_S3_uniq.bam -c ${aligDir}/Sample_1_S1_uniq.bam -f BAMPE --keep-dup auto -g hs --outdir ${varDir} -q 1 --name Sample_3_S3_narrow 

#Ctrl_H vs ctrl
$MACS_HOME/bin/macs2 callpeak -t ${aligDir}/Sample_4_S4_uniq.bam -c ${aligDir}/Sample_1_S1_uniq.bam -f BAMPE --keep-dup auto -g hs --outdir ${varDir} -q 1 --name CRTL_rnaseH_narrow 

#broad
$MACS_HOME/bin/macs2 callpeak -t ${aligDir}/Sample_2_S2_uniq.bam  -c ${aligDir}/Sample_1_S1_uniq.bam -f BAMPE --keep-dup auto -g hs --outdir ${varDir} -q 1 --broad --name Sample_2_S2_broad 
$MACS_HOME/bin/macs2 callpeak -t ${aligDir}/Sample_3_S3_uniq.bam  -c ${aligDir}/Sample_1_S1_uniq.bam -f BAMPE --keep-dup auto -g hs --outdir ${varDir} -q 1 --broad --name Sample_3_S3_broad 

#Ctrl_H vs ctrl
$MACS_HOME/bin/macs2 callpeak -t ${aligDir}/Sample_1_S1_uniq.bam -c ${aligDir}/Sample_4_S4_uniq.bam -f BAMPE --keep-dup auto -g hs --outdir ${varDir} -q 1 --broad --name rnaseH_vs_CTRL_broad 

wait

