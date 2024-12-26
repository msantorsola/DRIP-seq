#!/usr/bin/env Rscript



library(EnsDb.Hsapiens.v75) ##(hg19)
edb <- EnsDb.Hsapiens.v75
organism(edb)
#supportedFilters(edb)

library(ChIPpeakAnno)
library(GenomicFeatures)
library(biomaRt)
library(org.Hs.eg.db)


annoData <- toGRanges(edb,feature="gene")
seqlevelsStyle(annoData) <- "Ensembl"

idr_reg <- read.table('regular_noModel-idr', sep='\t')
colnames(idr_reg) = c('Chr','Start','End','Name','scaledIDR', 'Strand','signalValue', 'pvalue','qvalue','peak','globalIDR', 'localIDR', 'rep1_Start', 'rep1_End', 'rep1_signalValue', 'rep1_summit', 'rep2_Start', 'rep2_End', 'rep2_signalValue', 'rep2_summit')

#myPeakList= GRange object
idr_peaks <- toGRanges(idr_reg, format="BED")

seqlevelsStyle(idr_peaks) <- seqlevelsStyle(annoData)

#macs.anno <- annotatePeakInBatch(macsOutput, AnnotationData=TSS.human.GRCh38)
idr_anno <- annotatePeakInBatch(idr_peaks, AnnedbotationData=annoData, output="both", FeatureLocForDistance="geneEnd", bindingType=c("startSite", "endSite", "fullRange"), 
               maxgap=10000L, select="all")

head(idr_anno)


idr_anno <- as.data.frame(idr_anno)

idr_anno$gene_name <- annoData$gene_name[match(idr_anno$feature,names(annoData))]

write.table(idr_anno, 'regular_gene_anno.csv', sep='\t',row.names=F, quote=F)
idr_anno <- makeGRangesFromDataFrame(idr_anno, keep.extra.columns=TRUE)

png("DistributionAroundTSS.pdf")
binOverFeature(idr_anno, annotationData=annoData,
               radius=5000, nbins=20, FUN=length, errFun=0,
               ylab="count", 
               main="Distribution of aggregated peak numbers around TSS")
dev.off()

head(idr_anno)

over <- getEnrichedGO(idr_anno, orgAnn="org.Hs.eg.db", maxP=0.1, minGOterm=10, multiAdjMethod="BH",condense=FALSE)

write.table(over[["bp"]], 'regular_gene_BP_enrichment.csv', sep='\t',row.names=F, quote=F)
write.table(over[["cc"]], 'regular_gene_CC_enrichment.csv', sep='\t',row.names=F, quote=F)
write.table(over[["mf"]], 'regular_gene_MF_enrichment.csv', sep='\t',row.names=F, quote=F)


library(KEGG.db)
over_path <- getEnrichedPATH(idr_anno, orgAnn="org.Hs.eg.db", 
                    pathAnn="KEGG.db", maxP=0.05, minPATHterm=10, 
                    multiAdjMethod=NULL)

write.table(over_path, 'regular_gene_KEGG_enrichment.csv', sep='\t',row.names=F, quote=F)

#TSS annotations
data(TSS.human.GRCh37)

annotatedPeak <- annotatePeakInBatch(idr_peaks, AnnotationData=TSS.human.GRCh37)

annotatedPeak$gene_name <- annoData$gene_name[match(annotatedPeak$feature,names(annoData))]

write.table(annotatedPeak, 'regular_TSS_anno.csv', sep='\t',row.names=F, quote=F)

#enrichment
over_tss <- getEnrichedGO(annotatedPeak, orgAnn="org.Hs.eg.db", maxP=0.05, minGOterm=10, multiAdjMethod="BH", condense=FALSE)

write.table(over_tss[["bp"]], 'regular_TSS_BP_enrichment.csv', sep='\t',row.names=F, quote=F)
write.table(over_tss[["cc"]], 'regular_TSS_CC_enrichment.csv', sep='\t',row.names=F, quote=F)
write.table(over_tss[["mf"]], 'regular_TSS_MF_enrichment.csv', sep='\t',row.names=F, quote=F)


enrich_relaxed[["bp"]]$gene_name <- mapIds(org.Hs.eg.db, as.vector(enrich_relaxed[["bp"]]$EntrezID), 'SYMBOL', 'ENTREZID')
enrich_relaxed[["cc"]]$gene_name <- mapIds(org.Hs.eg.db, as.vector(enrich_relaxed[["cc"]]$EntrezID), 'SYMBOL', 'ENTREZID')
enrich_relaxed[["mf"]]$gene_name <- mapIds(org.Hs.eg.db, as.vector(enrich_relaxed[["mf"]]$EntrezID), 'SYMBOL', 'ENTREZID')


write.table(enrich_relaxed[['bp']], 'ALL_RELAXED_regular_BP_enrichment.csv', sep='\t',row.names=F, quote=F)
write.table(enrich_relaxed[['cc']], 'ALL_RELAXED_regular_CC_enrichment.csv', sep='\t',row.names=F, quote=F)
write.table(enrich_relaxed[['mf']], 'ALL_RELAXED_regular_MF_enrichment.csv', sep='\t',row.names=F, quote=F)


library(KEGG.db)
over_path_tss <- getEnrichedPATH(annotatedPeak, orgAnn="org.Hs.eg.db", pathAnn="KEGG.db", maxP=0.05, minPATHterm=10, multiAdjMethod=NULL)

enrich_relaxed_path <- getEnrichedPATH(relaxed_anno, orgAnn="org.Hs.eg.db", pathAnn="KEGG.db", maxP=0.05, minPATHterm=10, multiAdjMethod="BH")
enrich_relaxed_path$gene_name <- mapIds(org.Hs.eg.db, as.vector(enrich_relaxed_path$EntrezID), 'SYMBOL', 'ENTREZID')
write.table(enrich_relaxed_path, 'ALL_RELAXED_regular_KEGG_enrichment.csv', sep='\t',row.names=F, quote=F)



enrich_relaxed_path$gene_name <- mapIds(org.Hs.eg.db, as.vector(enrich_relaxed_path$EntrezID), 'SYMBOL', 'ENTREZID')
enrich_relaxed_path$PATH <- mapIds(reactome.db, as.vector(enrich_relaxed_path$path.id), 'PATHNAME', 'PATHID')
write.table(enrich_relaxed_path, 'ALL_RELAXED_regular_reactome_enrichment.csv', sep='\t',row.names=F, quote=F)



write.table(over_path_tss, 'regular_TSS_KEGG_enrichment.csv', sep='\t',row.names=F, quote=F)


#promoter
annotationData <- promoters(TSS.human.GRCh37, upstream=1000, downstream=500)
annotatedPeak_pr <- annotatePeakInBatch(idr_peaks, 
                                     AnnotationData=annotationData,
                                     output="both")

annotatedPeak_pr <- as.data.frame(annotatedPeak_pr)
annotatedPeak_pr$gene_description <- annotationData$description[match(annotatedPeak_pr$feature,names(annotationData))]

write.table(annotatedPeak_pr, 'regular_Promoter_anno.csv', sep='\t',row.names=F, quote=F)

#enrichment
over_pr <- getEnrichedGO(annotatedPeak_pr, orgAnn="org.Hs.eg.db", 
                    maxP=0.05, minGOterm=10, 
                    multiAdjMethod="BH",
                    condense=FALSE)

write.table(over_pr[["bp"]], 'regular_promoter_BP_enrichment.csv', sep='\t',row.names=F, quote=F)
write.table(over_pr[["cc"]], 'regular_promoter_CC_enrichment.csv', sep='\t',row.names=F, quote=F)
write.table(over_pr[["mf"]], 'regular_promoter_MF_enrichment.csv', sep='\t',row.names=F, quote=F)


library(KEGG.db)
over_path_tss <- getEnrichedPATH(annotatedPeak_pr, orgAnn="org.Hs.eg.db", 
                    pathAnn="KEGG.db", maxP=0.05, minPATHterm=10, 
                    multiAdjMethod=NULL)

write.table(over_path_tss, 'regular_promoter_KEGG_enrichment.csv', sep='\t',row.names=F, quote=F)


####   broad

annoData <- toGRanges(edb,feature="gene")
seqlevelsStyle(annoData) <- "Ensembl"

idr_bro <- read.table('broad_noModel-idr', sep='\t')
colnames(idr_bro) = c('Chr','Start','End','Name','scaledIDR', 'Strand','signalValue', 'pvalue','qvalue','globalIDR', 'localIDR', 'rep1_Start', 'rep1_End', 'rep1_signalValue', 'rep2_Start', 'rep2_End', 'rep2_signalValue')

#myPeakList= GRange object
idr_peaks <- toGRanges(idr_bro, format="BED")

seqlevelsStyle(idr_peaks) <- seqlevelsStyle(annoData)

#macs.anno <- annotatePeakInBatch(macsOutput, AnnotationData=TSS.human.GRCh38)
idr_anno <- annotatePeakInBatch(idr_peaks, AnnotationData=annoData, output="both", FeatureLocForDistance="geneEnd", bindingType=c("startSite", "endSite", "fullRange"), 
               maxgap=10000L, select="all")


head(idr_anno)

idr_anno <- as.data.frame(idr_anno)

idr_anno$gene_name <- annoData$gene_name[match(idr_anno$feature,names(annoData))]

write.table(idr_anno, 'broad_gene_anno.csv', sep='\t',row.names=F, quote=F)
idr_anno <- makeGRangesFromDataFrame(idr_anno, keep.extra.columns=TRUE)


head(idr_anno)

over <- getEnrichedGO(idr_anno, orgAnn="org.Hs.eg.db", 
                    maxP=0.05, minGOterm=10, 
                    multiAdjMethod="BH",
                    condense=FALSE)

write.table(over[["bp"]], 'broad_gene_enrichment.csv', sep='\t',row.names=F, quote=F)
write.table(over[["cc"]], 'broad_gene_CC_enrichment.csv', sep='\t',row.names=F, quote=F)
write.table(over[["mf"]], 'broad_gene_MF_enrichment.csv', sep='\t',row.names=F, quote=F)


library(KEGG.db)
over_path <- getEnrichedPATH(idr_anno, orgAnn="org.Hs.eg.db", 
                    pathAnn="KEGG.db", maxP=0.05, minPATHterm=10, 
                    multiAdjMethod=NULL)

write.table(over_path, 'broad_gene_KEGG_enrichment.csv', sep='\t',row.names=F, quote=F)

#TSS annotations
data(TSS.human.GRCh37)

annotatedPeak <- annotatePeakInBatch(idr_peaks, AnnotationData=TSS.human.GRCh37)

annotatedPeak$gene_name <- annoData$gene_name[match(annotatedPeak$feature,names(annoData))]

write.table(annotatedPeak, 'broad_TSS_anno.csv', sep='\t',row.names=F, quote=F)

#enrichment
over_tss <- getEnrichedGO(annotatedPeak, orgAnn="org.Hs.eg.db", 
                    maxP=0.05, minGOterm=10, 
                    multiAdjMethod="BH",
                    condense=FALSE)

write.table(over_tss[["bp"]], 'broad_TSS_BP_enrichment.csv', sep='\t',row.names=F, quote=F)
write.table(over_tss[["cc"]], 'broad_TSS_CC_enrichment.csv', sep='\t',row.names=F, quote=F)
write.table(over_tss[["mf"]], 'broad_TSS_MF_enrichment.csv', sep='\t',row.names=F, quote=F)


library(KEGG.db)
over_path_tss <- getEnrichedPATH(annotatedPeak, orgAnn="org.Hs.eg.db", 
                    pathAnn="KEGG.db", maxP=0.05, minPATHterm=10, 
                    multiAdjMethod=NULL)

write.table(over_path_tss, 'broad_TSS_KEGG_enrichment.csv', sep='\t',row.names=F, quote=F)



#promoter
annotationData <- promoters(TSS.human.GRCh37, upstream=1000, downstream=500)
annotatedPeak_pr <- annotatePeakInBatch(idr_peaks, 
                                     AnnotationData=annotationData,
                                     output="both")

annotatedPeak_pr$gene_description <- annotationData$description[match(annotatedPeak_pr$feature,names(annotationData))]

write.table(annotatedPeak_pr, 'broad_Promoter_anno.csv', sep='\t',row.names=F, quote=F)

#enrichment
over_pr <- getEnrichedGO(annotatedPeak_pr, orgAnn="org.Hs.eg.db", 
                    maxP=0.05, minGOterm=10, 
                    multiAdjMethod="BH",
                    condense=FALSE)

write.table(over_pr[["bp"]], 'broad_promoter_BP_enrichment.csv', sep='\t',row.names=F, quote=F)
write.table(over_pr[["cc"]], 'broad_promoter_CC_enrichment.csv', sep='\t',row.names=F, quote=F)
write.table(over_pr[["mf"]], 'broad_promoter_MF_enrichment.csv', sep='\t',row.names=F, quote=F)


library(KEGG.db)
over_path_tss <- getEnrichedPATH(annotatedPeak_pr, orgAnn="org.Hs.eg.db", 
                    pathAnn="KEGG.db", maxP=0.05, minPATHterm=10, 
                    multiAdjMethod=NULL)

write.table(over_path_tss, 'broad_promoter_KEGG_enrichment.csv', sep='\t',row.names=F, quote=F)


### PSEUDOREPLICATES
#pr_bro_rnM-idr #2076
#pr_rnM-idr #670 

annoData <- toGRanges(edb,feature="gene")
seqlevelsStyle(annoData) <- "Ensembl"

idr_reg <- read.table('pr_rnM-idr', sep='\t')
colnames(idr_reg) = c('Chr','Start','End','Name','scaledIDR', 'Strand','signalValue', 'pvalue','qvalue','peak','globalIDR', 'loc> IDR', 'rep1_Start', 'rep1_End', 'rep1_signalValue', 'rep1_summit', 'rep2_Start', 'rep2_End', 'rep2_signalValue', 'rep2_summit')

#myPeakList= GRange object
idr_peaks <- toGRanges(idr_reg, format="BED")

seqlevelsStyle(idr_peaks) <- seqlevelsStyle(annoData)

#macs.anno <- annotatePeakInBatch(macsOutput, AnnotationData=TSS.human.GRCh38)
idr_anno <- annotatePeakInBatch(idr_peaks, AnnotationData=annoData, output="both", FeatureLocForDistance="geneEnd", bindingType=c("startSite", "endSite", "fullRange"), 
               maxgap=10000L, select="all")

head(idr_anno)

idr_anno <- as.data.frame(idr_anno)

idr_anno$gene_name <- annoData$gene_name[match(idr_anno$feature,names(annoData))]

write.table(idr_anno, 'pr_regular_gene_anno.csv', sep='\t',row.names=F, quote=F)
idr_anno <- makeGRangesFromDataFrame(idr_anno, keep.extra.columns=TRUE)


head(idr_anno)

over <- getEnrichedGO(idr_anno, orgAnn="org.Hs.eg.db", 
                    maxP=0.05, minGOterm=10, 
                    multiAdjMethod="BH",
                    condense=FALSE)

write.table(over[["bp"]], 'pr_regular_gene_enrichment.csv', sep='\t',row.names=F, quote=F)
write.table(over[["cc"]], 'pr_regular_gene_CC_enrichment.csv', sep='\t',row.names=F, quote=F)
write.table(over[["mf"]], 'pr_regular_gene_MF_enrichment.csv', sep='\t',row.names=F, quote=F)


library(KEGG.db)
over_path <- getEnrichedPATH(idr_anno, orgAnn="org.Hs.eg.db", 
                    pathAnn="KEGG.db", maxP=0.05, minPATHterm=10, 
                    multiAdjMethod=NULL)

write.table(over_path, 'pr_regular_gene_KEGG_enrichment.csv', sep='\t',row.names=F, quote=F)

#TSS annotations
data(TSS.human.GRCh37)

annotatedPeak <- annotatePeakInBatch(idr_peaks, AnnotationData=TSS.human.GRCh37)

annotatedPeak$gene_name <- annoData$gene_name[match(annotatedPeak$feature,names(annoData))]

write.table(annotatedPeak, 'pr_regular_TSS_anno.csv', sep='\t',row.names=F, quote=F)

#enrichment
over_tss <- getEnrichedGO(annotatedPeak, orgAnn="org.Hs.eg.db", 
                    maxP=0.05, minGOterm=10, 
                    multiAdjMethod="BH",
                    condense=FALSE)

write.table(over_tss[["bp"]], 'pr_regular_TSS_BP_enrichment.csv', sep='\t',row.names=F, quote=F)
write.table(over_tss[["cc"]], 'pr_regular_TSS_CC_enrichment.csv', sep='\t',row.names=F, quote=F)
write.table(over_tss[["mf"]], 'pr_regular_TSS_MF_enrichment.csv', sep='\t',row.names=F, quote=F)


library(KEGG.db)
over_path_tss <- getEnrichedPATH(annotatedPeak, orgAnn="org.Hs.eg.db", 
                    pathAnn="KEGG.db", maxP=0.05, minPATHterm=10, 
                    multiAdjMethod=NULL)

write.table(over_path_tss, 'pr_regular_TSS_KEGG_enrichment.csv', sep='\t',row.names=F, quote=F)



#promoter
annotationData <- promoters(TSS.human.GRCh37, upstream=1000, downstream=500)
annotatedPeak_pr <- annotatePeakInBatch(idr_peaks, 
                                     AnnotationData=annotationData,
                                     output="both")

annotatedPeak_pr$gene_description <- annotationData$description[match(annotatedPeak_pr$feature,names(annotationData))]

write.table(annotatedPeak_pr, 'pr_regular_Promoter_anno.csv', sep='\t',row.names=F, quote=F)

#enrichment
over_pr <- getEnrichedGO(annotatedPeak_pr, orgAnn="org.Hs.eg.db", 
                    maxP=0.05, minGOterm=10, 
                    multiAdjMethod="BH",
                    condense=FALSE)

write.table(over_pr[["bp"]], 'pr_regular_promoter_BP_enrichment.csv', sep='\t',row.names=F, quote=F)
write.table(over_pr[["cc"]], 'pr_regular_promoter_CC_enrichment.csv', sep='\t',row.names=F, quote=F)
write.table(over_pr[["mf"]], 'pr_regular_promoter_MF_enrichment.csv', sep='\t',row.names=F, quote=F)


library(KEGG.db)
over_path_tss <- getEnrichedPATH(annotatedPeak_pr, orgAnn="org.Hs.eg.db", 
                    pathAnn="KEGG.db", maxP=0.05, minPATHterm=10, 
                    multiAdjMethod=NULL)

write.table(over_path_tss, 'pr_regular_promoter_KEGG_enrichment.csv', sep='\t',row.names=F, quote=F)


####   broad
annoData <- toGRanges(edb,feature="gene")
seqlevelsStyle(annoData) <- "Ensembl"

idr_bro <- read.table('pr_bro_rnM-idr', sep='\t')
colnames(idr_bro) = c('Chr','Start','End','Name','scaledIDR', 'Strand','signalValue', 'pvalue','qvalue','globalIDR', 'localIDR', 'rep1_Start', 'rep1_End', 'rep1_signalValue', 'rep2_Start', 'rep2_End', 'rep2_signalValue')

#myPeakList= GRange object
idr_peaks <- toGRanges(idr_bro, format="BED")

seqlevelsStyle(idr_peaks) <- seqlevelsStyle(annoData)

#macs.anno <- annotatePeakInBatch(macsOutput, AnnotationData=TSS.human.GRCh38)
idr_anno <- annotatePeakInBatch(idr_peaks, AnnotationData=annoData, output="both", FeatureLocForDistance="geneEnd", bindingType=c("startSite", "endSite", "fullRange"), 
               maxgap=10000L, select="all")


head(idr_anno)

idr_anno <- as.data.frame(idr_anno)

idr_anno$gene_name <- annoData$gene_name[match(idr_anno$feature,names(annoData))]

write.table(idr_anno, 'pr_broad_gene_anno.csv', sep='\t',row.names=F, quote=F)
idr_anno <- makeGRangesFromDataFrame(idr_anno, keep.extra.columns=TRUE)


head(idr_anno)

over <- getEnrichedGO(idr_anno, orgAnn="org.Hs.eg.db", 
                    maxP=0.05, minGOterm=10, 
                    multiAdjMethod="BH",
                    condense=FALSE)

write.table(over[["bp"]], 'pr_broad_gene_enrichment.csv', sep='\t',row.names=F, quote=F)
write.table(over[["cc"]], 'pr_broad_gene_CC_enrichment.csv', sep='\t',row.names=F, quote=F)
write.table(over[["mf"]], 'pr_broad_gene_MF_enrichment.csv', sep='\t',row.names=F, quote=F)


library(KEGG.db)
over_path <- getEnrichedPATH(idr_anno, orgAnn="org.Hs.eg.db", 
                    pathAnn="KEGG.db", maxP=0.05, minPATHterm=10, 
                    multiAdjMethod=NULL)

write.table(over_path, 'pr_broad_gene_KEGG_enrichment.csv', sep='\t',row.names=F, quote=F)

#TSS annotations
data(TSS.human.GRCh37)

annotatedPeak <- annotatePeakInBatch(idr_peaks, AnnotationData=TSS.human.GRCh37)

annotatedPeak$gene_name <- annoData$gene_name[match(annotatedPeak$feature,names(annoData))]

write.table(annotatedPeak, 'pr_broad_TSS_anno.csv', sep='\t',row.names=F, quote=F)

#enrichment
over_tss <- getEnrichedGO(annotatedPeak, orgAnn="org.Hs.eg.db", 
                    maxP=0.05, minGOterm=10, 
                    multiAdjMethod="BH",
                    condense=FALSE)

write.table(over_tss[["bp"]], 'pr_broad_TSS_BP_enrichment.csv', sep='\t',row.names=F, quote=F)
write.table(over_tss[["cc"]], 'pr_broad_TSS_CC_enrichment.csv', sep='\t',row.names=F, quote=F)
write.table(over_tss[["mf"]], 'pr_broad_TSS_MF_enrichment.csv', sep='\t',row.names=F, quote=F)


library(KEGG.db)
over_path_tss <- getEnrichedPATH(annotatedPeak, orgAnn="org.Hs.eg.db", 
                    pathAnn="KEGG.db", maxP=0.05, minPATHterm=10, 
                    multiAdjMethod=NULL)

write.table(over_path_tss, 'pr_broad_TSS_KEGG_enrichment.csv', sep='\t',row.names=F, quote=F)



#promoter
annotationData <- promoters(TSS.human.GRCh37, upstream=1000, downstream=500)
annotatedPeak_pr <- annotatePeakInBatch(idr_peaks, 
                                     AnnotationData=annotationData,
                                     output="both")

annotatedPeak_pr$gene_description <- annotationData$description[match(annotatedPeak_pr$feature,names(annotationData))]

write.table(annotatedPeak_pr, 'pr_broad_Promoter_anno.csv', sep='\t',row.names=F, quote=F)

#enrichment
over_pr <- getEnrichedGO(annotatedPeak_pr, orgAnn="org.Hs.eg.db", 
                    maxP=0.05, minGOterm=10, 
                    multiAdjMethod="BH",
                    condense=FALSE)

write.table(over_pr[["bp"]], 'pr_broad_promoter_BP_enrichment.csv', sep='\t',row.names=F, quote=F)
write.table(over_pr[["cc"]], 'pr_broad_promoter_CC_enrichment.csv', sep='\t',row.names=F, quote=F)
write.table(over_pr[["mf"]], 'pr_broad_promoter_MF_enrichment.csv', sep='\t',row.names=F, quote=F)


library(KEGG.db)
over_path_tss <- getEnrichedPATH(annotatedPeak_pr, orgAnn="org.Hs.eg.db", 
                    pathAnn="KEGG.db", maxP=0.05, minPATHterm=10, 
                    multiAdjMethod=NULL)

write.table(over_path_tss, 'pr_broad_promoter_KEGG_enrichment.csv', sep='\t',row.names=F, quote=F)


### true & pseudo replicates enrichment!

true_reg <- read.table('regular_noModel-idr', sep='\t')
colnames(true_reg) = c('Chr','Start','End','Name','scaledIDR', 'Strand','signalValue', 'pvalue','qvalue','peak','globalIDR', 'localIDR', 'rep1_Start', 'rep1_End', 'rep1_signalValue', 'rep1_summit', 'rep2_Start', 'rep2_End', 'rep2_signalValue', 'rep2_summit')
true_reg <- toGRanges(true_reg, format="BED")
seqlevelsStyle(true_reg) <- seqlevelsStyle(annoData)
true_reg_anno <- annotatePeakInBatch(true_reg, AnnotationData=annoData, output="both", FeatureLocForDistance="geneEnd", bindingType=c("startSite", "endSite", "fullRange"), 
               maxgap=10000L, select="all")

true_reg_anno <- as.data.frame(true_reg_anno)
true_reg_anno$gene_name <- annoData$gene_name[match(true_reg_anno$feature,names(annoData))]


#true_reg_anno <- unname(true_reg_anno)
pseudo_reg <- read.table('pr_rnM-idr', sep='\t')
colnames(pseudo_reg) = c('Chr','Start','End','Name','scaledIDR', 'Strand','signalValue', 'pvalue','qvalue','peak','globalIDR', 'localIDR', 'rep1_Start', 'rep1_End', 'rep1_signalValue', 'rep1_summit', 'rep2_Start', 'rep2_End', 'rep2_signalValue', 'rep2_summit')
pseudo_reg <- toGRanges(pseudo_reg, format="BED")
seqlevelsStyle(pseudo_reg) <- seqlevelsStyle(annoData)
pseudo_reg_anno <- annotatePeakInBatch(pseudo_reg, AnnotationData=annoData, output="both", FeatureLocForDistance="geneEnd", bindingType=c("startSite", "endSite", "fullRange"), 
               maxgap=10000L, select="all")

pseudo_reg_anno <- as.data.frame(pseudo_reg_anno)
pseudo_reg_anno$gene_name <- annoData$gene_name[match(pseudo_reg_anno$feature,names(annoData))]
#pseudo_reg_anno <- unname(pseudo_reg_anno)

all_reg <- rbind(true_reg_anno,pseudo_reg_anno)

all_reg  <- all_reg[!all_reg$seqnames=='chrY',]
all_reg_sign <- subset(all_reg, all_reg$scaledIDR >= 540)

#NO chrMT
all_reg_mt  <- all_reg_sign[!all_reg_sign$seqnames=='chrM',]
all_reg_mt <- makeGRangesFromDataFrame(all_reg_mt, keep.extra.columns=TRUE)

over_all_reg_mt <- getEnrichedGO(all_reg_mt, orgAnn="org.Hs.eg.db", maxP=0.05, minGOterm=10, multiAdjMethod=NULL,condense=FALSE)


over_all_reg_mt[["bp"]]$gene_name <- mapIds(org.Hs.eg.db, as.vector(over_all_reg_mt[["bp"]]$EntrezID), 'SYMBOL', 'ENTREZID')
over_all_reg_mt[["cc"]]$gene_name <- mapIds(org.Hs.eg.db, as.vector(over_all_reg_mt[["cc"]]$EntrezID), 'SYMBOL', 'ENTREZID')
over_all_reg_mt[["mf"]]$gene_name <- mapIds(org.Hs.eg.db, as.vector(over_all_reg_mt[["mf"]]$EntrezID), 'SYMBOL', 'ENTREZID')

write.table(over_all_reg_mt[["bp"]], 'significant-regular_gene_BP_enrichment_noChrM.csv', sep='\t',row.names=F, quote=F)
write.table(over_all_reg_mt[["cc"]], 'significant-regular_gene_CC_enrichment_noChrM.csv', sep='\t',row.names=F, quote=F)
write.table(over_all_reg_mt[["mf"]], 'significant-regular_gene_MF_enrichment_noChrM.csv', sep='\t',row.names=F, quote=F)


#BH not significant!!!
over_reg_path_mt <- getEnrichedPATH(all_reg_mt, orgAnn="org.Hs.eg.db", 
                    pathAnn="reactome.db", maxP=0.05, minPATHterm=10, 
                    multiAdjMethod=NULL)

over_reg_path_mt$gene_name <- mapIds(org.Hs.eg.db, as.vector(over_reg_path_mt$EntrezID), 'SYMBOL', 'ENTREZID')
over_reg_path_mt$PATH <- mapIds(reactome.db, as.vector(over_reg_path_mt$path.id), 'PATHNAME', 'PATHID')

write.table(over_reg_path_mt, 'significant-regular_gene_reactome_noChrM.csv', sep='\t',row.names=F, quote=F)

library(KEGG.db)
over_reg_kegg <- getEnrichedPATH(all_reg_mt, orgAnn="org.Hs.eg.db", 
                    pathAnn="KEGG.db", maxP=0.05, minPATHterm=10, 
                    multiAdjMethod=NULL)

over_reg_kegg$gene_name <- mapIds(org.Hs.eg.db, as.vector(over_reg_kegg$EntrezID), 'SYMBOL', 'ENTREZID')
#over_reg_path_mt$PATH <- mapIds(reactome.db, as.vector(over_reg_path_mt$path.id), 'PATHNAME', 'PATHID')

write.table(over_reg_path, 'significant-regular_gene_KEGG_enrichment.csv', sep='\t',row.names=F, quote=F)



#significant enrichment
#all_reg <- makeGRangesFromDataFrame(all_reg, keep.extra.columns=TRUE)
all_reg_sign <- makeGRangesFromDataFrame(all_reg_sign, keep.extra.columns=TRUE)

over_reg_sign <- getEnrichedGO(all_reg_sign, orgAnn="org.Hs.eg.db", 
                    maxP=0.05, minGOterm=10, 
                    multiAdjMethod="BH",
                    condense=FALSE)
over_reg_sign[["cc"]]$gene_name <- mapIds(org.Hs.eg.db, as.vector(over_reg_sign[["cc"]]$EntrezID), 'SYMBOL', 'ENTREZID')
over_reg_sign[["bp"]]$gene_name <- mapIds(org.Hs.eg.db, as.vector(over_reg_sign[["bp"]]$EntrezID), 'SYMBOL', 'ENTREZID')
over_reg_sign[["mf"]]$gene_name <- mapIds(org.Hs.eg.db, as.vector(over_reg_sign[["mf"]]$EntrezID), 'SYMBOL', 'ENTREZID')

write.table(over_reg_sign[["bp"]], 'significant-regular_gene_enrichment.csv', sep='\t',row.names=F, quote=F)
write.table(over_reg_sign[["cc"]], 'significant-regular_gene_CC_enrichment.csv', sep='\t',row.names=F, quote=F)
write.table(over_reg_sign[["mf"]], 'significant-regular_gene_MF_enrichment.csv', sep='\t',row.names=F, quote=F)


library(KEGG.db)
over_reg_kegg <- getEnrichedPATH(all_reg_sign, orgAnn="org.Hs.eg.db", 
                    pathAnn="reactome.db", maxP=0.05, minPATHterm=10, 
                    multiAdjMethod=NULL)

over_reg_path$gene_name <- mapIds(org.Hs.eg.db, as.vector(over_reg_path$EntrezID), 'SYMBOL', 'ENTREZID')

write.table(over_reg_path, 'significant-regular_gene_KEGG_enrichment.csv', sep='\t',row.names=F, quote=F)


##full enrichment
over_full <- getEnrichedGO(all_reg, orgAnn="org.Hs.eg.db", 
                    maxP=0.05, minGOterm=10, 
                    multiAdjMethod=N,
                    condense=FALSE)

write.table(over_full[["bp"]], 'full-regular_gene_enrichment.csv', sep='\t',row.names=F, quote=F)
write.table(over_full[["cc"]], 'full-regular_gene_CC_enrichment.csv', sep='\t',row.names=F, quote=F)
write.table(over_full[["mf"]], 'full-regular_gene_MF_enrichment.csv', sep='\t',row.names=F, quote=F)


library(KEGG.db)
over_full_path <- getEnrichedPATH(all_reg, orgAnn="org.Hs.eg.db", 
                    pathAnn="KEGG.db", maxP=0.05, minPATHterm=10, 
                    multiAdjMethod=NULL)

write.table(over_full_path, 'full-regular_gene_KEGG_enrichment.csv', sep='\t',row.names=F, quote=F)

###broad

true_bro <- read.table('broad_noModel-idr', sep='\t')
colnames(true_bro) = c('Chr','Start','End','Name','scaledIDR', 'Strand','signalValue', 'pvalue','qvalue','globalIDR', 'localIDR', 'rep1_Start', 'rep1_End', 'rep1_signalValue', 'rep2_Start', 'rep2_End', 'rep2_signalValue')
true_bro <- toGRanges(true_bro, format="BED")
seqlevelsStyle(true_bro) <- seqlevelsStyle(annoData)
true_bro_anno <- annotatePeakInBatch(true_bro, AnnotationData=annoData, output="both", FeatureLocForDistance="geneEnd", bindingType=c("startSite", "endSite", "fullRange"), 
               maxgap=10000L, select="all")
true_bro_anno <- as.data.frame(true_bro_anno)
true_bro_anno$gene_name <- annoData$gene_name[match(true_bro_anno$feature,names(annoData))]


pseudo_bro <- read.table('pr_bro_rnM-idr', sep='\t')
colnames(pseudo_bro) = c('Chr','Start','End','Name','scaledIDR', 'Strand','signalValue', 'pvalue','qvalue','globalIDR', 'localIDR', 'rep1_Start', 'rep1_End', 'rep1_signalValue', 'rep2_Start', 'rep2_End', 'rep2_signalValue')
pseudo_bro <- toGRanges(pseudo_bro, format="BED")
seqlevelsStyle(pseudo_bro) <- seqlevelsStyle(annoData)
pseudo_bro_anno <- annotatePeakInBatch(pseudo_bro, AnnotationData=annoData, output="both", FeatureLocForDistance="geneEnd", bindingType=c("startSite", "endSite", "fullRange"), 
               maxgap=10000L, select="all")
pseudo_bro_anno <- as.data.frame(pseudo_bro_anno)
pseudo_bro_anno$gene_name <- annoData$gene_name[match(pseudo_bro_anno$feature,names(annoData))]

all_bro <- rbind(true_bro_anno,pseudo_bro_anno)
all_bro <- all_bro[!all_bro$seqnames=='chrY',]
all_bro_sign <- subset(all_bro, all_bro$scaledIDR >= 540)
bro_mt <- all_bro_sign[!all_bro_sign$seqnames == 'chrM',]

#NO chrMT!!!

bro_mt <- makeGRangesFromDataFrame(bro_mt, keep.extra.columns=TRUE)
over_bro_mt <- getEnrichedGO(bro_mt, orgAnn="org.Hs.eg.db", 
                    maxP=0.05, minGOterm=10, 
                    multiAdjMethod=NULL,
                    condense=FALSE)

over_bro_mt[["cc"]]$gene_name <- mapIds(org.Hs.eg.db, as.vector(over_bro_mt[["cc"]]$EntrezID), 'SYMBOL', 'ENTREZID')
over_bro_mt[["bp"]]$gene_name <- mapIds(org.Hs.eg.db, as.vector(over_bro_mt[["bp"]]$EntrezID), 'SYMBOL', 'ENTREZID')
over_bro_mt[["mf"]]$gene_name <- mapIds(org.Hs.eg.db, as.vector(over_bro_mt[["mf"]]$EntrezID), 'SYMBOL', 'ENTREZID')


write.table(over_bro_mt[["bp"]], 'significant_broad_gene_BP_enrichment_noChrM.csv', sep='\t',row.names=F, quote=F)
write.table(over_bro_mt[["cc"]], 'significant_broad_gene_CC_enrichment_noChrM.csv', sep='\t',row.names=F, quote=F)
write.table(over_bro_mt[["mf"]], 'significant_broad_gene_MF_enrichment_noChrM.csv', sep='\t',row.names=F, quote=F)


over_bro_kegg <- getEnrichedPATH(bro_mt, orgAnn="org.Hs.eg.db", 
                    pathAnn="KEGG.db", maxP=0.05, minPATHterm=10, 
                    multiAdjMethod=NULL)

write.table(over_bro_sign_path, 'significant_broad_gene_KEGG_noChrM.csv', sep='\t',row.names=F, quote=F)


over_bro_reactome <- getEnrichedPATH(bro_mt, orgAnn="org.Hs.eg.db", 
                    pathAnn="reactome.db", maxP=0.05, minPATHterm=10, 
                    multiAdjMethod=NULL)

over_bro_reactome$gene_name <- mapIds(org.Hs.eg.db, as.vector(over_bro_reactome$EntrezID), 'SYMBOL', 'ENTREZID')
over_bro_reactome$PATH <- mapIds(reactome.db, as.vector(over_bro_reactome$path.id), 'PATHNAME', 'PATHID')

write.table(over_bro_sign_path, 'significant_broad_gene_REACTOME_noChrM.csv', sep='\t',row.names=F, quote=F)

#sgnificant enrichment
all_bro_sign <- makeGRangesFromDataFrame(all_bro_sign, keep.extra.columns=TRUE)
over_bro_sign <- getEnrichedGO(all_bro_sign, orgAnn="org.Hs.eg.db", 
                    maxP=0.05, minGOterm=10, 
                    multiAdjMethod="BH",
                    condense=FALSE)

over_bro_sign[["cc"]]$gene_name <- mapIds(org.Hs.eg.db, as.vector(over_bro_sign[["cc"]]$EntrezID), 'SYMBOL', 'ENTREZID')
over_bro_sign[["bp"]]$gene_name <- mapIds(org.Hs.eg.db, as.vector(over_bro_sign[["bp"]]$EntrezID), 'SYMBOL', 'ENTREZID')
over_bro_sign[["mf"]]$gene_name <- mapIds(org.Hs.eg.db, as.vector(over_bro_sign[["mf"]]$EntrezID), 'SYMBOL', 'ENTREZID')


write.table(over_bro_sign[["bp"]], 'significant_broad_gene_BP_enrichment.csv', sep='\t',row.names=F, quote=F)
write.table(over_bro_sign[["cc"]], 'significant_broad_gene_CC_enrichment.csv', sep='\t',row.names=F, quote=F)
write.table(over_bro_sign[["mf"]], 'significant_broad_gene_MF_enrichment.csv', sep='\t',row.names=F, quote=F)


library(KEGG.db)

over_bro_sign_path <- getEnrichedPATH(all_bro_sign, orgAnn="org.Hs.eg.db", 
                    pathAnn="KEGG.db", maxP=0.05, minPATHterm=10, 
                    multiAdjMethod="BH")


over_bro_sign_path$gene_name <- mapIds(org.Hs.eg.db, as.vector(over_bro_sign_path$EntrezID), 'SYMBOL', 'ENTREZID')
write.table(over_bro_sign_path, 'significant_broad_gene_KEGG_enrichment.csv', sep='\t',row.names=F, quote=F)


over_bro_sign_path <- getEnrichedPATH(all_bro_sign, orgAnn="org.Hs.eg.db", 
                    pathAnn="reactome.db", maxP=0.05, minPATHterm=10, 
                    multiAdjMethod="BH")

over_bro_sign_path$gene_name <- mapIds(org.Hs.eg.db, as.vector(over_bro_sign_path$EntrezID), 'SYMBOL', 'ENTREZID')
write.table(over_bro_sign_path, 'significant_broad_gene_REACTOME_enrichment.csv', sep='\t',row.names=F, quote=F)




over_bro <- getEnrichedGO(all_bro, orgAnn="org.Hs.eg.db", 
                    maxP=0.05, minGOterm=10, 
                    multiAdjMethod="BH",
                    condense=FALSE)

write.table(over_bro[["bp"]], 'full_broad_gene_BP_enrichment.csv', sep='\t',row.names=F, quote=F)

write.table(over_bro[["cc"]], 'full_broad_gene_CC_enrichment.csv', sep='\t',row.names=F, quote=F)

write.table(over_bro[["mf"]], 'full_broad_gene_MF_enrichment.csv', sep='\t',row.names=F, quote=F)


library(KEGG.db)

over_bro_path <- getEnrichedPATH(all_bro, orgAnn="org.Hs.eg.db", 
                    pathAnn="KEGG.db", maxP=0.05, minPATHterm=10, 
                    multiAdjMethod=NULL)


write.table(over_bro_path, 'full_broad_gene_KEGG_enrichment.csv', sep='\t',row.names=F, quote=F)
















###Find possible enhancers depend on DNA interaction data 
 
fantom <- read.table("/marconi_work/uTS18_Schoeftn_0/data_annotation/enhancer_tss_associations.bed", stringsAsFactors = FALSE)

colnames(fantom) <- c("chrom","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","itemRgb","blockCount","blockSizes","chromStarts")

head(fantom,3)

library(splitstackshape)
fantom <- as.data.frame(cSplit(fantom, splitCols = "name", sep = ";", direction = "wide"))

locs <- strsplit(as.character(fantom$name_1),"[:-]")
fantom$chr <- sapply(locs,"[",1)
fantom$start <- as.numeric(sapply(locs,"[",2))
fantom$end <- as.numeric(sapply(locs,"[",3))
fantom$symbol <- fantom$name_3
fantom$corr <- sub("R:","", fantom$name_4)
fantom$fdr <- sub("FDR:","", fantom$name_5)

fantom <- unique(subset(fantom, corr >= 0.25 & fdr <1e-5, select = c("chr","start","end","symbol","corr","name_2")))
write.table(fantom, 'enhancer_promoter_corr025_fdr1e-5.bed', sep='\t',row.names=F, quote=F)


fantom <- makeGRangesFromDataFrame(fantom, keep.extra.columns = TRUE)
#fantom <- unlist(liftOver(fantom, ch)) to hg38!
seqlevelsStyle(fantom) <- "Ensembl"

seqlevelsStyle(fantom) <- seqlevelsStyle(fantom)
#enhancer - promoter pairs with a decent level of correlation and significance and tidy the data at the same time:
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75
organism(edb)
#supportedFilters(edb)

library(ChIPpeakAnno)
library(GenomicFeatures)
library(biomaRt)
library(org.Hs.eg.db)


annoData <- toGRanges(edb,feature="gene")
seqlevelsStyle(annoData) <- "Ensembl"

idr_reg <- read.table('regular_noModel-idr', sep='\t')
colnames(idr_reg) = c('Chr','Start','End','Name','scaledIDR', 'Strand','signalValue', 'pvalue','qvalue','peak','globalIDR', 'localIDR', 'rep1_Start', 'rep1_End', 'rep1_signalValue', 'rep1_summit', 'rep2_Start', 'rep2_End', 'rep2_signalValue', 'rep2_summit')

#myPeakList= GRange object
idr_peaks <- toGRanges(idr_reg, format="BED")
seqlevelsStyle(idr_peaks) <- "Ensembl"


#hits <- findOverlaps(snps_hard, fantom)
hit_fantom <- annotatePeakInBatch(idr_peaks, AnnotationData=fantom, output="both", maxgap=10000L, select="all")


hit_fantom$symbol <- fantom$symbol[match(hit_fantom$feature,names(annoData))]

download.file("http://www.cell.com/cms/attachment/2086554122/2074217047/mmc4.zip", destfile = "mmc4.zip")
# uncomment the following lines to extract zipped files
#unzip("mmc4.zip")
#unzip("DATA_S1.zip")
pchic <- read.delim("ActivePromoterEnhancerLinks.tsv", stringsAsFactors = FALSE)
head(pchic,3)
                   
tsss <- promoters((TSS.human.GRCh37, upstream=1000, downstream=500, columns = "gene_id")
hits <- nearest(promoters, tsss)
pchic$gene_id <- unlist(tsss[hits]$gene_id)
                    


