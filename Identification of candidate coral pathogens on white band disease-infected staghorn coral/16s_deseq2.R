rm(list=ls())
require(DESeq2)

setwd("~/Desktop/cbbs4/cbbs4_97open_otus_coldose2/")

#this is the otu table from qiime converted to .csv
c=read.csv("complete_otu_table_frombiom.csv", h=T, row.names=1)

#information about samples
s=read.csv("dosecol_data.csv", h=T, as.is=T)
b=t(c)
d=cbind.data.frame(s,b)

#check order of samples
check<-cbind(row.names(d), d[1:1])

# Remove first column
d=d[2:(NCOL(d))]

#subset dosage and collection data
dosage<-subset(d, Experiment=="Dosage")
collection<-subset(d, Experiment=="Collection" & Sample_name!="bad")

#analyze collection  data
d=collection
c=d[6:NCOL(d)]
c=t(c)
s=d[0:5]

#filter out OTUS present only in other dataset
rs = rowSums (c)
use = (rs > 0)
c = c [ use, ]

#create count data set 
dds <-DESeqDataSetFromMatrix(countData = c, colData=s, 
                             design= ~ Disease_state+Year)

#set as factors
colData(dds)$Disease_state <- factor(colData(dds)$Disease_state, levels=c("Healthy","Disease"))

colData(dds)$Year <- factor(colData(dds)$Year, levels=c("nine","ten"))



#estimate size factors (normalize)
dds <- estimateSizeFactors(dds)

#filter otus <10 times
rs = rowSums ( counts ( dds))
use = (rs > 10)
dds = dds [ use, ]

#estimate dispersions
dds <- estimateDispersions(dds)

#look at dispersons 
plotDispEsts(dds)

#perform wald test for significance
dds <- nbinomWaldTest(dds, pAdjustMethod="fdr")

#plot significant OTUS 
plotMA(dds, ylim=c(-10,10))
rm(list=ls())

resultsNames(dds)

#significant OTUs dvh
rescol<-results(dds, "Disease_state_Disease_vs_Healthy")
ressigcol <- subset(rescol, padj<0.05)
ressigcol<-as.data.frame(ressigcol)

#rename columns
names(ressiginter)[names(ressiginter)=="log2FoldChange"] <- "log2FC_collection"
names(ressiginter)[names(ressiginter)=="baseMean"] <- "baseMean_colleciton"
names(ressiginter)[names(ressiginter)=="padj"] <- "padj_collection"
ressiginter$pvalue<-NULL
ressiginter$lfcSE<-NULL

#significant OTUS 9v10
resyear<-results(dds, "Year_ten_vs_nine")
ressigyear <- subset(resyear, padj<0.05)
ressigyear<-as.data.frame(ressigyear)

#rename columns
names(ressigyearint)[names(ressigyearint)=="log2FoldChange"] <- "log2FC_year"
names(ressigyearint)[names(ressigyearint)=="baseMean"] <- "baseMean_year"
names(ressigyearint)[names(ressigyearint)=="padj"] <- "padj_year"
ressigyearint$pvalue<-NULL
ressigyearint$lfcSE<-NULL

#significant OTUS interaction
resinter<-results(dds, "Disease_stateDisease.Yearten")
ressiginter <- subset(resinter, padj<0.05)
ressiginter<-as.data.frame(ressiginter)

#rename columns
names(ressiginter)[names(ressiginter)=="log2FoldChange"] <- "log2FC_interaction"
names(ressiginter)[names(ressiginter)=="baseMean"] <- "baseMean_interaction"
names(ressiginter)[names(ressiginter)=="padj"] <- "padj_interaction"
ressiginter$pvalue<-NULL
ressiginter$lfcSE<-NULL

# Create and write normalized counts
comm.normalized=t(counts(dds, normalized=TRUE))
# Save the normalized data
write.csv(file="counts_norm_collection.csv", comm.normalized, row.names=T)

#analyze dose data
d=dosage
c=d[6:NCOL(d)]
c=t(c)
s=d[0:5]

#create count data set and run DESEQ2 on dosage
dds <-DESeqDataSetFromMatrix(countData = c, colData=s, 
                             design= ~ Disease_state)

#set as factors
colData(dds)$Disease_state <- factor(colData(dds)$Disease_state, levels=c("Healthy","Disease"))

#estimate size factors (normalize)
dds <- estimateSizeFactors(dds)

#filter otus <10 times
rs = rowSums ( counts ( dds))
use = (rs > 10)
dds = dds [ use, ]

#estimate dispersions
dds <- estimateDispersions(dds)

#look at dispersons 
plotDispEsts(dds)

#perform wald test for significance
dds <- nbinomWaldTest(dds, pAdjustMethod="fdr")

#plot significant OTUS 
plotMA(dds, ylim=c(-10,10))
rm(list=ls())

resdose<-results(dds)
ressigdose<-subset(resdose, padj<0.05)
ressigdose<-as.data.frame(ressigdose)

#rename columns
names(ressigdose)[names(ressigdose)=="log2FoldChange"] <- "log2FC_tank"
names(ressigdose)[names(ressigdose)=="baseMean"] <- "baseMean_tank"
names(ressigdose)[names(ressigdose)=="padj"] <- "padj_tank"
ressigdose$pvalue<-NULL
ressigdose$lfcSE<-NULL

# Create and write normalized counts
comm.normalized=t(counts(dds, normalized=TRUE))

# Save the normalized data
write.csv(file="counts_norm_dose.csv", comm.normalized, row.names=T)

#find overlap of col and dose
coldose<-merge(ressigcol, ressigdose, by.x="row.names", by.y="row.names", all.x=TRUE, all.y=TRUE)
yearint<-merge(ressigyearint, ressiginter, by.x="row.names", by.y="row.names", all.x=TRUE, all.y=TRUE)
coldoseyearint<-merge(coldose, yearint, by.x="Row.names", by.y="Row.names", all.x=TRUE, all.y=TRUE)

#import and assigntaxonomy
t=read.csv("new_repset_tax_assignments.csv", h=T)
coldoseyearinttax<-merge(coldoseyearint, t, by.x = "Row.names", by.y="OTUID", all.x= TRUE)

write.csv(as.data.frame(coldoseyearinttax), file="supplementary_table_1.csv")







