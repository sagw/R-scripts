rm(list=ls())
require(vegan)
require(DESeq2)
setwd("~/Documents/SD1/SD1_analyses/")
c=read.csv("input_files//SD1_OTU_table_97.csv", h=T, row.names=1)

#information about samples
s=read.csv("input_files//SD1_map.csv", h=T, as.is=T)
b=t(c6)
d=cbind.data.frame(s,b)


#check order of samples
check<-cbind(row.names(d), d[1:1])

t=read.csv("input_files//rep_set_tax_assignments.csv", header=T)

# Remove first column

d=d[2:(NCOL(d))]

#filter out chloroplast DNA
dtax<-merge(t, c6, by.x="OTUID", by.y="row.names", all.y=T)

dtax1<-dtax[- grep("Chloroplast", dtax$Tax),]
row.names(dtax1)<-dtax1$OTUID
c<-dtax1[5:NCOL(dtax1)]

write.csv(c, file="dtax_test.csv")

b=t(c)
d=cbind.data.frame(s,b)

timetwo<-subset(d, Timepoint=="two")   

d<-timetwo
c<-d[10:NCOL(d)]
rs = rowSums (c)
use = (rs > 1)
c = c [ use, ]
c<-t(c)
s<-d[0:9]

dds <-DESeqDataSetFromMatrix(countData = c, colData=s, 
                             design= ~ Site*Dose_disease_state*Dose.Site)
colData(dds)$Dose_disease_state <- factor(colData(dds)$Dose_disease_state, levels=c("Healthy","Diseased"))
colData(dds)$Site <- factor(colData(dds)$Site, levels=c("CK4","CK14"))
colData(dds)$Dose.Site <- factor(colData(dds)$Dose.Site, levels=c("CK4","CK14"))

dds<-estimateSizeFactors(dds)
rs = rowSums (counts(dds))
use = (rs > 1)
dds = dds [ use, ]

dds<-DESeq(dds)
dds<-estimateDispersions(dds)
plotDispEsts(dds)
dds <- nbinomWaldTest(dds, maxit=100)

resultsNames(dds)

plotMA(dds, ylim=c(-10,10), "Dose_disease_stateDiseased.Dose.SiteCK14", pvalCutoff=0.05)

plotMA(dds, ylim=c(-10,10), "Dose_disease_state_Diseased_vs_Healthy", pvalCutoff=0.05)

comm.normalized=(counts(dds, normalized=TRUE))

res<-results(dds, name="Dose_disease_state_Diseased_vs_Healthy")
ressig <- subset(res, padj<0.05)
ressig<-as.data.frame(ressig)
ressigtax<-merge(ressig, t, by.x="row.names", by.y="OTUID", all.x=T)
write.csv(as.data.frame(ressigtax), file="Output_files/T10_Dose_disease_state.csv")

ressigotu<-merge(ressigtax, comm.normalized, by.x="Row.names", by.y="row.names")

write.csv(as.data.frame(ressigotu), file="Output_files/T10_Dose_disease_state_otutable.csv")

resdosesite<-results(dds, name="Dose_disease_stateDiseased.Dose.SiteCK14")
resdosesitesig <- subset(resdosesite, padj<0.05)
resdosesitesig<-as.data.frame(resdosesitesig)
resdosesitetaxsig<-merge(resdosesitesig, t, by.x="row.names", by.y="OTUID", all.x=T)

all<-merge(resdosesitesig, ressig, by="row.names")
alltax<-merge(all, t, by.x="Row.names", by.y="OTUID")