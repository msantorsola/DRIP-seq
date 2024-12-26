# DRIP summary


#Merge overlapping peaks
#"Chr","start","end","name","score","strand","signalValue","pValue","qValue","width"

sort -k1,1 -k2,2n bro_rep1_sign_q05.bed > bro_rep1_sign_q05_sorted.bed

sort -k1,1 -k2,2n bro_rep2_sign_q05.bed > bro_rep2_sign_q05_sorted.bed


bedtools merge -c 4,7,9,10,10 -o count,sum,collapse,mean,collapse -i bro_rep1_sign_q05_sorted.bed -delim "|"  > bro_rep1_sign_q05_merged.bed

bedtools merge -c 4,7,9,10,10 -o count,sum,collapse,mean,collapse -i bro_rep2_sign_q05_sorted.bed -delim "|"  > bro_rep2_sign_q05_merged.bed



#PEAK size
br2 <- read.table("bro_rep2_sign_q05_sorted.bed", header=F, sep='\t')
rr2 <- read.table("reg_rep2_sign_q05_sorted.bed", header=F, sep='\t')
colnames(br2) <- c("Chr","start","end","name","score","strand","signalValue","pValue","qValue","width")
colnames(rr2) <-  c("Chr","start","end","name","score","strand","signalValue","pValue","qValue","summit","width")


rr2$type <- "Narrow"
rr2$rep <- "siSFPQ_D"
br2$type <- "Broad"
br2$rep <- "siSFPQ_D"


br1 <- read.table("bro_rep1_sign_q05_sorted.bed", header=T, sep='\t')
rr1 <- read.table("reg_rep1_sign_q05_sorted.bed", header=T, sep='\t')
r1 <- rr1[,c(4,11,13)]
b1 <- br1[,c(4,10,12)]
all <- rbind(r1,b1)


png("PeakSizeDistribution_siSFPQM.png", height=20, width=20, res=600, units="cm")

boxplot(all$width~all$type,  col=c( "lightskyblue","lightskyblue2"), ylab="Peak size (bp)", ylim=c(0, 2000), main= "siSFPQM")

dev.off()



### Total read counts

count <- read.table("ReadCountInPeaks.txt", fill=T, header=T)


png("TotalReadCount.png", height=20, width=20, res=600, units="cm")
ggplot(count) +
  geom_bar(aes(x=sample, y=nReads), stat="identity", fill=c("red3", "darkgoldenrod1","mediumblue"), alpha=1, width=0.5, colour="black") +
  ggtitle("") +
  #ggtitle("Total Read Count") +
  ylab("Read count (million)") +
  xlab("") +
  theme_bw() + #blank area of graph
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + #blank grid
  theme(text = element_text(size=40),panel.border = element_blank(), axis.line = element_line()) +
  scale_y_continuous(expand = c(0,0)) + #no space between origin and bars #show only x and y axes and not square around area plot!
  theme(aspect.ratio = 2/1.5, axis.text.x = element_text(angle = 45, hjust = 1))
  
dev.off()
 


###	Gene-level pie chart

library(dplyr)
library(scales)
 
count <- read.table("Gene-level_regular.txt", fill=T)
count$V1 <- factor(count$V1, levels = c("promoter", "terminal", "promoter-terminal", "gene-body", "intergenic"))
  
count %>%
arrange(desc(count$V2)) %>%
mutate(prop = percent(count$V2 / sum(count$V2))) -> count 


png("Gene-level_regular_pie.png", height=40, width=40, res=600, units="cm")

ggplot(count, aes(x = "", y = count$V2, fill = count$V1)) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  coord_polar("y", start = 0)+
  #geom_text_repel(aes(label = prop)) +
  #geom_text(aes(y = prop, label = count$V2)  +
  geom_text(aes(x=1.6, label = count$V2), position = position_stack(vjust = 0.5), color = "black", size=15)+
  theme_void() +
  theme(legend.title = element_blank(), legend.text=element_text(size=9))

dev.off()


#      Relationship to genomic context



#PEAK NUMBER
count <- read.table("Rep1_Rep2_peak_number.txt", fill=T, header=T, sep='\t')
colnames(count) <- c("Sample","Peak_type","Peak_number")


png("Rep1_Rep2_peak_number.png", height=20, width=20, res=600, units="cm")


ggplot(count, aes(x=Sample, y=Peak_number, fill=Peak_type)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  ggtitle(" ") +
  ylab("Peak #") +
  xlab("") +
  theme_bw() + #blank area of graph
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Blues") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + #blank grid
  theme(text = element_text(size=20),panel.border = element_blank(), axis.line = element_line()) +
  scale_y_continuous(expand = c(0,0)) + #no space between origin and bars #show only x and y axes and not square around area plot!
  theme(aspect.ratio = 2/1, axis.text.x = element_text(angle = 90, hjust = 1))
  
dev.off()
  
#Intersection
count <- read.table("Rep1_Rep2_Intersection_peak_number.txt", fill=T, header=T, sep='\t')
colnames(count) <- c("Sample","Peak_type","Peak_number")


png("Rep1_Rep2_Intersection_peak_number.png", height=20, width=20, res=600, units="cm")


ggplot(count, aes(x=Sample, y=Peak_number, fill=Peak_type)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  ggtitle(" ") +
  ylab("Peak #") +
  xlab("") +
  theme_bw() + #blank area of graph
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Greens") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + #blank grid
  theme(text = element_text(size=20),panel.border = element_blank(), axis.line = element_line()) +
  scale_y_continuous(expand = c(0,0)) + #no space between origin and bars #show only x and y axes and not square around area plot!
  theme(aspect.ratio = 2/1, axis.text.x = element_text(angle = 45, hjust = 1))
  
dev.off()
  
#IDR
count <- read.table("Rep1_Rep2_IDR_peak_number.txt", fill=T, header=T, sep='\t')
colnames(count) <- c("Sample","Peak_type","Peak_number")


png("Rep1_Rep2_IDR_peak_number.png", height=20, width=20, res=600, units="cm")


ggplot(count, aes(x=Sample, y=Peak_number, fill=Peak_type)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  ggtitle(" ") +
  ylab("Peak #") +
  xlab("") +
  theme_bw() + #blank area of graph
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Reds") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + #blank grid
  theme(text = element_text(size=20),panel.border = element_blank(), axis.line = element_line()) +
  scale_y_continuous(expand = c(0,0)) + #no space between origin and bars #show only x and y axes and not square around area plot!
  theme(aspect.ratio = 2/1, axis.text.x = element_text(angle = 45, hjust = 1))
  
dev.off()


# regular
chrom <- read.table("Chromatin_regular_intersection_GAT.txt", header=T, sep='\t')

chrom$Track <- factor(chrom$Track, levels = rev(chrom$Track)) 


png("regular_intersection_ChromHMM_GAT-analysis.png", height=30, width=40, res=600, units='cm')
ggplot(chrom, aes(Track,l2fold)) +
    geom_bar(aes(fill = qvalue), stat="identity",position="identity", color = "black") +
    coord_flip() +
    scale_fill_distiller(palette = "RdBu") +
    theme(text = element_text(size=20)) +
    theme_bw() + #blank area of graph
    theme(#panel.grid.major = element_blank(), 
          #panel.grid.minor = element_blank(), 
          #axis.line.x=element_blank(),
          axis.text.x=element_text(size = 12),
          #axis.ticks.x=element_blank(),
          axis.title.x=element_text(size = 20),
          #panel.grid.minor.x=element_blank(), #blank grid
          #panel.grid.major.x=element_blank(),
          axis.text.y = element_text(size = 25)) + 
    ylab("log2(fc)") +
    xlab("") +
    ggtitle("")


dev.off()




#Simple Repeats


chrom <- read.table("Repeats_regular_intersection_GAT.txt", header=T, sep='\t')
head(chrom)
chrom <- chrom[order(chrom$l2fold,decreasing = TRUE),]
chrom$Class <- factor(chrom$Class, levels = rev(chrom$Class))  #Then turn it back into a factor with the levels in the correct order (the same order in df, FDR >!)


png("regular_intersection_SimpleRepeats_GAT_Analysis.png", height=30, width=40, res=600, units='cm')
ggplot(chrom, aes(Class,l2fold)) +
    geom_bar(aes(fill = qvalue), stat="identity",position="identity", color = "black") +
    coord_flip() +
    scale_fill_distiller(palette = "RdBu") +
    theme(text = element_text(size=15)) +
    theme_bw() + #blank area of graph
    theme(#panel.grid.major = element_blank(), 
          #panel.grid.minor = element_blank(), 
          #axis.line.x=element_blank(),
          axis.text.x=element_text(size = 12),
          #axis.ticks.x=element_blank(),
          axis.title.x=element_text(size = 20),
          #panel.grid.minor.x=element_blank(), #blank grid
          #panel.grid.major.x=element_blank(),
          axis.text.y = element_text(size = 25)) + 
    ylab("log2(fc)") +
    xlab("") +
    ggtitle("")


dev.off()





# Peak annotation and Pathway Enrichment
#annotated files
        #broad_intersection_q05_merged_Final.bed
        #regular_intersection_q05_merged_Final.csv 


library(EnsDb.Hsapiens.v75) #Citation (from within R, enter citation("EnsDb.Hsapiens.v75")): Rainer J (2017). EnsDb.Hsapiens.v75: Ensembl based annotation package. R package version 2.99.0.
edb <- EnsDb.Hsapiens.v75
organism(edb)
library(ChIPpeakAnno)
library(GenomicFeatures)
library(biomaRt)
library(org.Hs.eg.db)
library(rtracklayer)


annoData <- toGRanges(edb,feature="gene")
seqlevelsStyle(annoData) <- "Ensembl"




# 1 broad
broad <- read.table("broad_intersection_q05_merged_Final.bed", sep='\t')
#Chr start end merged_peaks signalValue qValue mean_width single_width
colnames(broad) = c('Chr','Start','End', 'merged_peaks','peak_name','signalValue_sum', 'qValue','mean_width', 'width')
summary(broad)


broad <- toGRanges(broad, format="BED")
seqlevelsStyle(broad) <- seqlevelsStyle(annoData)
anno <- annotatePeakInBatch(broad, AnnotationData=annoData, FeatureLocForDistance="TSS", PeakLocForDistance="start", output="both", maxgap=10000L)
#“TSS” (default) means using the start of the feature when the feature is on plus strand and using the end of feature when the feature is on minus strand
anno$gene_name <- annoData$gene_name[match(anno$feature, names(annoData))]
head(anno)
write.table(anno, 'broad_intersection_annotations.csv', sep='\t',row.names=F, quote=F)


#Functional annotation and enrichment plots
#KEGG
#### Broad


BP <- read.table('broad_intersection_KEGG_Enrichment', sep='\t', header=T, quote=NULL)
BP_up <- subset(BP, BP$FDR.100 <0.05) #FDR by david have to be /100
#BP_up$X.log10.FDR. <- with(BP_up, reorder(X.log10.FDR., X.log10.FDR., function(x) -length(x)))
#BP_up$Term <- as.character(BP_up$Term)
#BP_up$X.log10.FDR. <- as.numeric(BP_up$X.log10.FDR.)
library(tidyr)
BP_up <- separate(data = BP_up, col = Term, into = c("ID", "Description"), sep = "\\:")
BP_up$Description<- factor(BP_up$Description, levels = rev(BP_up$Description))  #Then turn it back into a factor with the levels in the correct order (the same order in df, FDR >!)

png("broad_KEGG_Enrichment_15_dic.png", height=30, width=50, res=600, units='cm')
ggplot(BP_up, aes(Description,X.log10.FDR.)) +
    geom_bar(fill = "lightsteelblue", stat="identity", position = "dodge", color = "black", width = 0.6) +
    #geom_bar(width = 1, stat = "identity", color = "black")
    coord_flip() +
    #scale_fill_distiller(palette = "RdBu") +
    theme(text = element_text(size=35), axis.text.y = element_text(angle=0, hjust=1, size=35), plot.title = element_text(hjust = 0.5, size=15)) +
    theme_bw() + #blank area of graph
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))  +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line.x=element_blank(),
          #axis.text.x=element_blank(),
          #axis.ticks.x=element_blank(),
          axis.title.x=element_text(size = 20),
          panel.grid.minor.x=element_blank(), #blank grid
          panel.grid.major.x=element_blank(),
          #aspect.ratio = 2/1.5,
          axis.text.y = element_text(size = 40)) + 
    theme(legend.position = "none") +
    #scale_y_continuous(breaks = round(seq(min(as.numeric(BP_up$X.log10.FDR.)), max(as.numeric(BP_up$X.log10.FDR.)), by = 0.5),1)) +
    #scale_x_continuous(breaks = seq(0, 2.5, len = 0.5)) +
    ylab("-log10(FDR)") +
    xlab("") +
    ggtitle("")
 
dev.off()       


#Category        Term        Count        %        PValue        Genes        List Total        Pop Hits        Pop Total        Fold Enrichment        Bonferroni        Benjamini        FDR        FDR/100        -log10(FDR)
BP <- read.table('/Users/mariangela/Documents/DRIP-Seq/annotation_july19/broad_intersection_GO-BP_enrichment.csv', sep='\t', header=T, quote=NULL)
BP_up <- subset(BP, BP$FDR.100 <0.05) #FDR by david have to be /100
#BP_up <- BP_up[1:20,]
#BP_up$X.log10.FDR. <- as.numeric(BP_up$X.log10.FDR.)
#BP_up$X.log10.FDR. <- with(BP_up, reorder(X.log10.FDR., X.log10.FDR., function(x) -length(x)))
#BP_up$Term <- as.character(BP_up$Term)
library(tidyr)
BP_up <- separate(data = BP_up, col = Term, into = c("ID", "Description"), sep = "\\~")
BP_up$Description<- factor(BP_up$Description, levels = rev(BP_up$Description))


png("broad_GO-BP_Enrichment_15_dic.png", height=30, width=50, res=600, units='cm')
ggplot(BP_up, aes(Description,X.log10.FDR.)) +
    geom_bar(fill = "darkolivegreen3", stat="identity", position = "dodge", color = "black") +
    #geom_bar(width = 1, stat = "identity", color = "black")
    coord_flip() +
    #scale_fill_distiller(palette = "RdBu") +
    theme(text = element_text(size=35), axis.text.y = element_text(angle=0, hjust=1, size=35), plot.title = element_text(hjust = 0.5, size=15)) +
    theme_bw() + #blank area of graph
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))  +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line.x=element_blank(),
          #axis.text.x=element_blank(),
          #axis.ticks.x=element_blank(),
          axis.title.x=element_text(size = 20),
          panel.grid.minor.x=element_blank(), #blank grid
          panel.grid.major.x=element_blank(),
          #aspect.ratio = 2/1.5,
          axis.text.y = element_text(size = 40)) + 
    theme(legend.position = "none") +
    #scale_y_continuous(breaks = round(seq(min(as.numeric(BP_up$X.log10.FDR.)), max(as.numeric(BP_up$X.log10.FDR.)), by = 0.5),1)) +
    #scale_x_continuous(breaks = seq(0, 2.5, len = 0.5)) +
    ylab("-log10(FDR)") +
    xlab("") +
    ggtitle("")
 
dev.off() 




BP <- read.table('broad_intersection_GO-MF_enrichment.csv', sep='\t', header=T, quote=NULL)
BP_up <- subset(BP, BP$FDR.100 <0.05) #FDR by david have to be /100
#BP_up <- BP_up[1:20,]
#BP_up$X.log10.FDR. <- as.numeric(BP_up$X.log10.FDR.)
#BP_up$X.log10.FDR. <- with(BP_up, reorder(X.log10.FDR., X.log10.FDR., function(x) -length(x)))
#BP_up$Term <- as.character(BP_up$Term)
BP_up <- separate(data = BP_up, col = Term, into = c("ID", "Description"), sep = "\\~")
BP_up$Description<- factor(BP_up$Description, levels = rev(BP_up$Description))  #Then turn it back into a factor with the levels in the correct order (the same order in df, FDR >!)


#png("broad_GO-MF_Enrichment.png", height=15, width=30, res=600, units='cm')
go_mf <- ggplot(BP_up, aes(Description,X.log10.FDR.)) +
    geom_bar(fill = "lightsteelblue2", stat="identity", position = "dodge", color = "black") +
    #geom_bar(width = 1, stat = "identity", color = "black")
    coord_flip() +
    #scale_fill_distiller(palette = "RdBu") +
    theme(text = element_text(size=6), axis.text.y = element_text(angle=0, hjust=1, size=20), plot.title = element_text(hjust = 0.5, size=15)) +
    theme_bw() + #blank area of graph
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))  +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line.x=element_blank(),
          #axis.text.x=element_blank(),
          #axis.ticks.x=element_blank(),
          #axis.title.x=element_blank(),
          panel.grid.minor.x=element_blank(), #blank grid
          panel.grid.major.x=element_blank(),
          axis.text.y = element_text(size = 12)) + 
    theme(legend.position = "none") +
    #scale_y_continuous(breaks = round(seq(min(as.numeric(BP_up$X.log10.FDR.)), max(as.numeric(BP_up$X.log10.FDR.)), by = 0.5),1)) +
    #scale_x_continuous(breaks = seq(0, 2.5, len = 0.5)) +
    ylab("-log10(FDR)") +
    xlab("") +
    ggtitle("GO Molecular Function")
 
dev.off() 







