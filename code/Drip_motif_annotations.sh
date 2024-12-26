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
module load autoload bedtools2
module load autoload meme


#create a logfile 
exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1>motif_annotations.log 2>&1


bedtools(){
	$BEDTOOLS2_HOME/bin/bedtools "$@"
}
export bedtools


mkdir -p $WORK/drip_project_test/MACS_peaks/homer_anno
chmod -R 755 $WORK/drip_project_test/MACS_peaks/homer_anno


mkdir -p $WORK/drip_project_test/MACS_peaks/meme_motif
chmod -R 755 $WORK/drip_project_test/MACS_peaks/meme_motif

workDir=$WORK/drip_project_test/MACS_peaks/
memeDir=$WORK/drip_project_test/MACS_peaks/meme_motif
homerDir=$WORK/drip_project_test/MACS_peaks/homer_anno

cd ${workDir}


# bed to fasta

declare -a peaks=( "regular_peaks_Model_peaks" "regular_peaks_noModel_peaks" "broad_peaks_Model_peaks" "broad_peaks_Model_peaks")

for peak in "${peaks[@]}"

	do

	bedtools getfasta -fi $IG_HG19_GENOME -bed ${peak}.narrowPeak -fo ${peak}.narrowPeak.fa

	bedtools getfasta -fi $IG_HG19_GENOME -bed ${peak}.gappedPeak -fo ${peak}.gappedPeak.fa

	bedtools getfasta -fi $IG_HG19_GENOME -bed ${peak}.broadPeak -fo ${peak}.broadPeak.fa


	meme ${peak}.narrowPeak.fa -dna -mod anr -nmotifs 10 -oc ${memeDir}

	meme ${peak}.gappedPeak.fa -dna -mod anr -nmotifs 10 -oc ${memeDir}

	meme ${peak}.broadPeak.fa -dna -mod anr -nmotifs 10 -oc ${memeDir}

done





# Annotate regions with HOMER


#MACS
annotatePeaks.pl MACS_peakCalling_noModel_peaks.narrowPeak hg19   > MACS_HOMER_peakAnnotation.txt 
#GEM
annotatePeaks.pl /lustrehome/m.santorsola/ChIP-Seq/data/GEM_peakCalling/GEM_peakCalling.GEM_events.narrowPeak hg19   > /lustrehome/m.santorsola/ChIP-Seq/data/OVCAR_P3/GEM_HOMER_peakAnnotation.txt 
#GPS
annotatePeaks.pl /lustrehome/m.santorsola/ChIP-Seq/data/GEM_peakCalling/GEM_peakCalling.GPS_events.narrowPeak  hg19 > /lustrehome/m.santorsola/ChIP-Seq/data/OVCAR_P3/GPS_HOMER_peakAnnotation.txt 
#SPP
annotatePeaks.pl /lustrehome/m.santorsola/ChIP-Seq/data/OVCAR_P3/OVCAR_P3_S1_uniquely.mapped_VS_OVCAR_CTRL_S2_uniquely.mapped.regionPeak hg19   > /lustrehome/m.santorsola/ChIP-Seq/data/OVCAR_P3/SPP_HOMER_peakAnnotation.txt 

