#!/bin/bash


#SBATCH --time 20:00:00
#SBATCH --job-name=rep2_231
#SBATCH --partition=gll_usr_prod
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --mem=86000
#SBATCH --output=rep2_231.%J.out
#SBATCH --error=rep2_231.%J.err
#SBATCH --qos=gll_qos_shared
#SBATCH --account=***

#SBATCH --mail-user=***@gmail.com
#SBATCH --mail-type=ALL,TIME_LIMIT_90



module load profile/bioinf
module load autoload fastqc/0.11.5
module load autoload trimmomatic
module load autoload ig_UCSC_Homo_sapiens/hg19
module load autoload bwa
module load autoload picardtools/2.3.0
module load autoload samtools
module load autoload sambamba


workDir=***/drip_2021_231
tmpDir=***/drip_2021_231/tmp
aligDir=***/drip_2021_231/alignment
logDir=***/drip_2021_231/logs
qualDir=***/drip_2021_231/qual
varDir=***/drip_2021_231/peaks


#create a logfile 
exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1>rep2_231.log 2>&1


cd ${workDir}

ctrl='Sample_1_S1'
ctrl_H='Sample_4_S4'
rep1='Sample_2_S2'
rep2='Sample_3_S3'


#check quality of your data with fastqc
declare -a fastqs=("Sample_3_S3_L001_R1_001.fastq.gz " "Sample_3_S3_L001_R2_001.fastq.gz" "Sample_3_S3_L002_R1_001.fastq.gz" "Sample_3_S3_L002_R2_001.fastq.gz" "Sample_3_S3_L003_R1_001.fastq.gz" "Sample_3_S3_L003_R2_001.fastq.gz" "Sample_3_S3_L004_R1_001.fastq.gz" "Sample_3_S3_L004_R2_001.fastq.gz")

for fastq in "${fastqs[@]}"

do

fastqc -t 10 ${fastq}

done



#with trimmomatic
java -jar -Xmx40g -XX:ParallelGCThreads=20 -Djava.io.tmpdir=${tmpDir} $TRIMMOMATIC_HOME/trimmomatic-0.36.jar PE -threads 20  ${workDir}/${rep2}_L001_R1_001.fastq.gz ${workDir}/${rep2}_L001_R2_001.fastq.gz ${aligDir}/${rep2}_L001_1_cleaned.fastq.gz ${aligDir}/${rep2}_L001_1_discarded.fastq.gz ${aligDir}/${rep2}_L001_2_cleaned.fastq.gz ${aligDir}/${rep}_L001_2_discarded.fastq.gz ILLUMINACLIP:$TRIMM_ADAPTERS/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:18 MINLEN:25
java -jar -Xmx40g -XX:ParallelGCThreads=20 -Djava.io.tmpdir=${tmpDir} $TRIMMOMATIC_HOME/trimmomatic-0.36.jar PE -threads 20  ${workDir}/${rep2}_L002_R1_001.fastq.gz ${workDir}/${rep2}_L002_R2_001.fastq.gz ${aligDir}/${rep2}_L002_1_cleaned.fastq.gz ${aligDir}/${rep2}_L002_1_discarded.fastq.gz ${aligDir}/${rep2}_L002_2_cleaned.fastq.gz ${aligDir}/${rep}_L002_2_discarded.fastq.gz ILLUMINACLIP:$TRIMM_ADAPTERS/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:18 MINLEN:25
java -jar -Xmx40g -XX:ParallelGCThreads=20 -Djava.io.tmpdir=${tmpDir} $TRIMMOMATIC_HOME/trimmomatic-0.36.jar PE -threads 20  ${workDir}/${rep2}_L003_R1_001.fastq.gz ${workDir}/${rep2}_L003_R2_001.fastq.gz ${aligDir}/${rep2}_L003_1_cleaned.fastq.gz ${aligDir}/${rep2}_L003_1_discarded.fastq.gz ${aligDir}/${rep2}_L003_2_cleaned.fastq.gz ${aligDir}/${rep}_L003_2_discarded.fastq.gz ILLUMINACLIP:$TRIMM_ADAPTERS/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:18 MINLEN:25
java -jar -Xmx40g -XX:ParallelGCThreads=20 -Djava.io.tmpdir=${tmpDir} $TRIMMOMATIC_HOME/trimmomatic-0.36.jar PE -threads 20  ${workDir}/${rep2}_L004_R1_001.fastq.gz ${workDir}/${rep2}_L004_R2_001.fastq.gz ${aligDir}/${rep2}_L004_1_cleaned.fastq.gz ${aligDir}/${rep2}_L004_1_discarded.fastq.gz ${aligDir}/${rep2}_L004_2_cleaned.fastq.gz ${aligDir}/${rep}_L004_2_discarded.fastq.gz ILLUMINACLIP:$TRIMM_ADAPTERS/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:18 MINLEN:25


#align reads
$BWA_HOME/bin/bwa mem -t 40 -T 0 -M -R '@RG\tID:${rep2}_L001\tLB:${rep2}\tSM:${rep2}\tPL:ILLUMINA' $BWA_INDEX/genome.fa ${aligDir}/${rep2}_L001_1_cleaned.fastq.gz  ${aligDir}/${rep2}_L001_2_cleaned.fastq.gz > ${aligDir}/${rep2}_L001.sam
samtools view -hbu ${aligDir}/${rep2}_L001.sam  > ${aligDir}/${rep2}_L001.bam
java -jar -Xmx40g -XX:ParallelGCThreads=20 -Djava.io.tmpdir=${tmpDir} $PICARDTOOLS_HOME/bin/picard.jar SortSam CREATE_INDEX=true INPUT=${aligDir}/${rep2}_L001.bam OUTPUT=${aligDir}/${rep2}_L001_sorted.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=STRICT

$BWA_HOME/bin/bwa mem -t 40 -T 0 -M -R '@RG\tID:${rep2}_L002\tLB:${rep2}\tSM:${rep2}\tPL:ILLUMINA' $BWA_INDEX/genome.fa ${aligDir}/${rep2}_L002_1_cleaned.fastq.gz  ${aligDir}/${rep2}_L002_2_cleaned.fastq.gz > ${aligDir}/${rep2}_L002.sam
samtools view -hbu ${aligDir}/${rep2}_L002.sam  > ${aligDir}/${rep2}_L002.bam
java -jar -Xmx40g -XX:ParallelGCThreads=20 -Djava.io.tmpdir=${tmpDir} $PICARDTOOLS_HOME/bin/picard.jar SortSam CREATE_INDEX=true INPUT=${aligDir}/${rep2}_L002.bam OUTPUT=${aligDir}/${rep2}_L002_sorted.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=STRICT

$BWA_HOME/bin/bwa mem -t 40 -T 0 -M -R '@RG\tID:${rep2}_L003\tLB:${rep2}\tSM:${rep2}\tPL:ILLUMINA' $BWA_INDEX/genome.fa ${aligDir}/${rep2}_L003_1_cleaned.fastq.gz  ${aligDir}/${rep2}_L003_2_cleaned.fastq.gz > ${aligDir}/${rep2}_L003.sam
samtools view -hbu ${aligDir}/${rep2}_L003.sam  > ${aligDir}/${rep2}_L003.bam
java -jar -Xmx40g -XX:ParallelGCThreads=20 -Djava.io.tmpdir=${tmpDir} $PICARDTOOLS_HOME/bin/picard.jar SortSam CREATE_INDEX=true INPUT=${aligDir}/${rep2}_L003.bam OUTPUT=${aligDir}/${rep2}_L003_sorted.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=STRICT

$BWA_HOME/bin/bwa mem -t 40 -T 0 -M -R '@RG\tID:${rep2}_L004\tLB:${rep2}\\tSM:${rep2}\tPL:ILLUMINA' $BWA_INDEX/genome.fa ${aligDir}/${rep2}_L004_1_cleaned.fastq.gz  ${aligDir}/${rep2}_L004_2_cleaned.fastq.gz > ${aligDir}/${rep2}_L004.sam
samtools view -hbu ${aligDir}/${rep2}_L004.sam  > ${aligDir}/${rep2}_L004.bam
java -jar -Xmx40g -XX:ParallelGCThreads=20 -Djava.io.tmpdir=${tmpDir} $PICARDTOOLS_HOME/bin/picard.jar SortSam CREATE_INDEX=true INPUT=${aligDir}/${rep2}_L004.bam OUTPUT=${aligDir}/${rep2}_L004_sorted.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=STRICT

#merge bam files
java -Xmx40g -Djava.io.tmpdir=${tmpDir} -jar $PICARDTOOLS_HOME/bin/picard.jar MergeSamFiles I=${aligDir}/${rep2}_L001_sorted.bam I=${aligDir}/${rep2}_L002_sorted.bam I=${aligDir}/${rep2}_L003_sorted.bam I=${aligDir}/${rep2}_L004_sorted.bam O=${aligDir}/${rep2}_merged.bam SO=coordinate CREATE_INDEX=TRUE


# Marking Duplicates
java -jar -Xmx40g -XX:ParallelGCThreads=20 -Djava.io.tmpdir=\${tmpDir} $PICARDTOOLS_HOME/bin/picard.jar MarkDuplicates CREATE_INDEX=true I=${aligDir}/${rep2}_merged.bam O=${aligDir}/${rep2}_dedup.bam M=${aligDir}/${rep2}_dedup_metrics.txt VALIDATION_STRINGENCY=LENIENT TMP_DIR=${tmpDir} REMOVE_DUPLICATES=true


#Filtering uniquely mapped and properly paired reads
$SAMBAMBA_HOME/bin/sambamba-0.7.0-linux-static  view -t 12 -h -f bam -F "mapping_quality >= 1 and not (unmapped or secondary_alignment) and not ([XA] != null or [SA] != null)" ${aligDir}/${rep2}_dedup.bam -o ${aligDir}/${rep2}_uniq.bam


wait



