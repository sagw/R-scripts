rm(list=ls())
#open DESeq2
require(DESeq2)

setwd("~/Documents/seastar_transcriptome")

#Import your counts and goslim info
counts_goslim=read.csv("counts_then_goslim.csv", header=T)

#Import your metadata (site, treatment, etc)
data=read.csv("ss_data.csv", header=T, row.names=1)

#Import counts
counts=read.table("Phel_clc_RNAseq_count.txt", header=T, row.names=1)

#cut out the goslim info
geneinfo=counts_goslim[7:NCOL(counts_goslim)]

#Create Countdataset
dds <-DESeqDataSetFromMatrix(countData = counts, colData=data, 
                             design= ~ innoculated)
#set variables as factors
colData(dds)$innoculated <- factor(colData(dds)$innoculated, levels=c("C","T"))

colData(dds)$innoculated <- factor(colData(dds)$Location_collected, levels=c("DB","PH","FHL"))

colData(dds)$innoculated <- factor(colData(dds)$days_to_sick, levels=c("NA","1","20"))

#estimate size factors (normalize)
dds <- estimateSizeFactors(dds)

#estimate dispersions
dds <- estimateDispersions(dds)

#check dispersions
plotDispEsts(dds)

#run stats
dds <- nbinomWaldTest(dds, pAdjustMethod="fdr")

#plot sig DE
plotMA(dds, ylim=c(-10,10))

#make results file
res<-results(dds)
#order by padj value
res <-res[order(res$padj),]
#make results a data frame
res<-as.data.frame(res)

#Subset significant values only
ressig <- subset(res, padj<0.05)
#make data frame
ressig<-as.data.frame(ressig)
#merge list of significant genes with GO information
ressig_go<-merge(ressig, geneinfo, by.x="row.names", by.y="Column1", all.x=TRUE)
#write csv
write.csv(as.data.frame(ressig_go), file="treatment_sig_goslim.csv")