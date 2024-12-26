#!/bin/bash

#SBATCH --time 20:00:00
#SBATCH --job-name=fastqc_trimming   
#SBATCH -p knl_usr_prod
#SBATCH -N 1
#SBATCH -n 20 
#SBATCH --mem=86000
#SBATCH --output=fastqc_trimming.%J.out
#SBATCH --error=fastqc_trimming.%J.err
#SBATCH --account=<ACCOUNT_NAME>


module load profile/global
module load profile/bioinf

module load autoload trimmomatic
module load fastqc


#create a logfile 
exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1>fastqc.log 2>&1

trimmomatic(){
      java -jar $TRIMMOMATIC_HOME/trimmomatic-0.36.jar "$@"
}
export trimmomatic


#create and go to the working directory
mkdir -p $WORK/drip_project 
chmod -R 755 $WORK/drip_project

mkdir -p $WORK/drip_project/logs $WORK/drip_project/qual $WORK/drip_project/alignment $WORK/drip_project/peak_calling
#mkdir -p $WORK/drip_project_test/fastqc

chmod -R 755 $WORK/drip_project/logs $WORK/drip_project/qual $WORK/drip_project/alignment $WORK/drip_project/peak_calling $WORK/drip_project/tmp


workDir=$WORK/drip_project/
tmpDir=$WORK/drip_project/tmp
aligDir=$WORK/drip_project/alignment
logDir=$WORK/drip_project/logs
qualDir=$WORK/drip_project/qual
peakDir=$WORK/drip_project/peak_calling
rawDir=$WORK/drip_raw

cd ${workDir}

#link the bam locally
#ln -s $WORK/drip_seq.samples/*fastq.gz .
ln -s ${rawDir}/*fastq.gz .


#QC and cleaning 

#declare sample name
declare -a arr=()

for sample in "${arr[@]}"     
      do

      echo "Running FastQC on raw data..."
      #check quality of raw data with fastqc
      fastqc -t 2 ${sample}_R1_001.fastq.gz
      fastqc -t 2 ${sample}_R2_001.fastq.gz

      
      #SAMPLE=$(echo ${i} | sed "s/_R1_001\.fastq\.gz//")
      
      echo "Running trimmomatic on raw data..."
      #with trimmotatic
      trimmomatic PE -threads 20 ${sample}_R1_001.fastq.gz ${sample}_R2_001.fastq.gz ${rawDir}/${sample}_R1.cleaned.fastq.gz ${rawDir}/${sample}_R1.discarded.fastq.gz ${rawDir}/${sample}_R2.cleaned.fastq.gz ${rawDir}/${sample}_R2.discarded.fastq.gz ILLUMINACLIP:$TRIMM_ADAPTERS/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:18 MINLEN:50

      echo "Running FastQC on cleaned data..."
      #check quality of clean data with fastqc
      fastqc -t 2 ${rawDir}/${sample}_R1.cleaned.fastq.gz
      fastqc -t 2 ${rawDir}/${sample}_R2.cleaned.fastq.gz

done

echo "Moving FastQC output..."
mv *fastqc* ${qualDir}

echo "Removing discarded reads..."
rm ${rawDir}/*.discarded.fastq.gz

echo "done"