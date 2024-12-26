#!/bin/bash

#SBATCH --time 20:00:00
#SBATCH --job-name=align_filter  
#SBATCH -p knl_usr_prod
#SBATCH -N 1
#SBATCH -n 20 
#SBATCH --mem=86000
#SBATCH --output=align_filter.%J.out
#SBATCH --error=align_filter.%J.err
#SBATCH --account=<ACCOUNT_NAME>


module load profile/global
module load profile/bioinf

module load autoload ig_homo_sapiens/hg19
module load autoload bwa
module load autoload samtools
module load autoload picard
module load autoload gatk
module load autoload sambamba/0.6.6




bwa(){
      $BWA_HOME/bin/bwa mem -t 40 -M -R "@RG\tID:${sample}\tLB:${sample}\tSM:${sample}\tPL:ILLUMINA" $IG_HG19_BWA_INDEX/genome.fa "$@"
}
export bwa


samtools(){
     $SAMTOOLS_HOME/bin/samtools  "$@"
}
export samtools


markdupl(){
      java -jar -Xmx40g -XX:ParallelGCThreads=20 $PICARD_HOME/picard.jar MarkDuplicates "$@"
}
export markdupl


mergebam(){
      java -Xmx40g -Djava.io.tmpdir=$WORK/tmp_dir -jar $PICARD_HOME/picard.jar MergeSamFiles "$@"
}
export mergebam


sambamba(){
      $SAMBAMBA_HOME/bin/sambamba_v0.6.6 "$@" 
}
export sambamba


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
exec 1>${logDir}/align_filter.log 2>&1

#declare sample names
declare -a arr=( "siCon_S9_6_1_RV_S1_L001" "siCon_S9_6_1_RV_S1_L002" "siCon_S9_6_1_RV_S1_L003" "siCon_S9_6_1_RV_S1_L004" "siPSFD_S9_6_3_RV_S2_L001" "siPSFD_S9_6_3_RV_S2_L002" "siPSFD_S9_6_3_RV_S2_L003" "siPSFD_S9_6_3_RV_S2_L004" "siPSFM_S9_6_2_RV_S1_L001" "siPSFM_S9_6_2_RV_S1_L002" "siPSFM_S9_6_2_RV_S1_L003" "siPSFM_S9_6_2_RV_S1_L004")

#sort and index bam
for sample in "${arr[@]}"

      do
      
      echo "Mapping with bwa..."
      #align reads
      bwa ${rawDir}/${sample}_R1.cleaned.fastq.gz ${rawDir}/${sample}_R2.cleaned.fastq.gz >  ${rawDir}/${sample}.sam
      
      # Remove unmapped reads
      #samtools view -hS -F 4 /lustrehome/m.santorsola/ChIP-Seq/data/OVCAR_CTRL/OVCAR_CTRL_S2_L001.sam | /lustre/browser/bin/samtools1.3.1/bin/samtools view  -b  -S  >  /lustrehome/m.santorsola/ChIP-Seq/data/OVCAR_CTRL/OVCAR_CTRL_S2_L001_best.mapped_only.bam 

      echo "sam 2 bam & sorting bam..."
      # Sort BAM file  
      samtools view -hbu ${rawDir}/${sample}.sam  | samtools sort -T FIX -o ${aligDir}/${sample}.sorted.bam
      #chmod -R 755 ${sample}.sorted.bam
      
      echo "Indexing bam files by genomic coordinates..."
      # index the bam file  
      samtools index -b ${aligDir}/${sample}.sorted.bam
      #chmod -R 755 ${sample}.sorted.bam.bai

done


#merge bam files
for i in "${arr[@]}"

      do

      echo "Merging lanes per sample..."
      mergebam \
            I=${aligDir}/${i}_L001.sorted.bam \
            I=${aligDir}/${i}_L002.sorted.bam \
            I=${aligDir}/${i}_L003.sorted.bam \
            I=${aligDir}/${i}_L004.sorted.bam \
            O=${aligDir}/${i}_merged.bam \
            SO=coordinate \
            CREATE_INDEX=TRUE

      echo "Marking Duplicates..."
      # Marking Duplicates
      markdupl I=${aligDir}/${i}_merged.bam O=${aligDir}/${i}_dedup.bam M=${aligDir}/${i}_dedup_metrics.txt CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=tmpDir 
      
      echo "Filtering uniquely mapped reads..."
      # Filter uniquely mapped reads
      #samtools view -h ${i}_dedup.bam | grep -v XA:Z | grep -v SA:Z | /lustre/browser/bin/samtools1.3.1/bin/samtools view -b  > /lustrehome/m.santorsola/ChIP-Seq/data/OVCAR_P3_S1_uniquely.mapped.bam 
      sambamba view -t 12 -h -f bam -F "mapping_quality >= 1 and not (unmapped or secondary_alignment) and not ([XA] != null or [SA] != null)" ${aligDir}/${i}_dedup.bam -o ${aligDir}/${i}_uniq.bam
     
done

echo "Removing sam files..."
rm ${rawDir}/*.sam

echo "done"


