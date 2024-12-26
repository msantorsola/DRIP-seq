#!/usr/bin/env python

import deeptools.countReadsPerBin as crpb

#The FRiP score as a measure of ChIP-seq quality
#FRiP is very useful for comparing results obtained with the same antibody across cell lines or 
#with different antibodies against the same factor. 
#only for point-source data sets, we calculate the fraction of all mapped reads that fall into peak regions identified by a peak-calling algorithm 


#bed file
#.narrowPeak format  is  a  extension  of  BED  format,  you  will  need  to  convert  the  0-base  start  coordinates
#encoded  in  the  BED  format  to  1-base  start  coordinates  (just  sum  +1  to  every  start  coordinate


# not recommended for gappedPeak and broadPeak!


bed_files = ["peaks.bed"]

#bam_file1 = 
#bam_file2 = 

cr = countReadsPerBin.CountReadsPerBin([bam_file1, bam_file2],
                                        bedFile=bed_files,
                                        numberOfProcessors=10)
reads_at_peaks = cr.run()
#print reads_at_peaks
#The result is a numpy array with a row for each peak region and a column for each BAM file.

reads_at_peaks.shape

#(6295, 2)

#Now, the *total number of reads per peaks* per bam file is computed:

total = reads_at_peaks.sum(axis=0)

#Next, we need to find the *total number of mapped reads* in each of the bam files. For this we use the pysam module.

import pysam
bam1 = pysam.AlignmentFile(bam_file1)
bam2 = pysam.AlignmentFile(bam_file2)

#Now, bam1.mapped and bam2.mapped contain the total number of mapped reads in each of the bam files, respectively.
#Finally, we can compute the FRiP score:

frip1 = float(total[0]) / bam1.mapped
frip2 = float(total[1]) / bam2.mapped
print frip1, frip2

#0.170030741997, 0.216740390353
