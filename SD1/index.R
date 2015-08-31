rm(list=ls())
require(DESeq2)
require(genefilter)
require(ape)
setwd("~/Documents/SD1/SD1_analyses/")
c=read.csv("input_files//filtered_otutable.csv", h=T, row.names=1)
c=c[5:NCOL(c)]

#information about samples
s=read.csv("input_files//SD1_map.csv", h=T, as.is=T)
b=t(c)
d=cbind.data.frame(s,b)
d=d[!d$Sample_name=="SD1.CK14D2.CK14.O",]
d=d[d$Timepoint %in% c("one","two","three"),]

c=d[15:NCOL(d)]
c=t(c)
rs = rowSums (c)
use = (rs > 0)
c = c [ use, ]

s=d[0:14]

#normalize using deseq
dds <-DESeqDataSetFromMatrix(countData = c, colData=s, 
                             design= ~ T3_disease_state)

colData(dds)$T3_disease_state <- factor(colData(dds)$T3_disease_state, levels=c("Healthy","Diseased"))

#calculate different geometric means because of zeros.
gm_mean=function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0])))
geoMeans = apply(c, 1, gm_mean) 

#use locfunc=shorth instead of median because of low counts.
dds = estimateSizeFactors(dds, geoMeans=geoMeans, locfunc=shorth)

comm.normalized=t(counts(dds, normalized=TRUE))

row.names(comm.normalized)=colnames(c)
comm.normalized=t(comm.normalized)

#import the filtered to 20,000 otus
filtered=read.csv("input_files/otus_filtered20000.csv", h=T, row.names=1)

commfilt=merge(filtered, comm.normalized, by.x="OTUID", by.y="row.names")
row.names(commfilt)=commfilt$OTUID
commfilt=commfilt[6:NCOL(commfilt)]
commfilt=t(commfilt)

#create normalized counts with info
dnorm=merge(s, commfilt, by.x="Sample_name", by.y="row.names")

#subset time one
one=dnorm[dnorm$Timepoint=="one",]
c1=one[14:NCOL(one)]
row.names(c1)=one$Sample_name
s1=one[0:13]

#subset time two d
twod=dnorm[dnorm$Timepoint=="two"& dnorm$T3_disease_state=="Diseased",]
cd2=twod[14:NCOL(twod)]
row.names(cd2)=twod$Sample_name
sd2=twod[0:13]

#filter out otus not in those two times
combined=rbind(c1, cd2)
combined=t(combined)
rs = rowSums (combined)
use = (rs > 0)
combined = combined [ use, ]

#separate one and two again
c1=combined[,0:82]
c1=t(c1)
cd2=combined[,83:NCOL(combined)]
cd2=t(cd2)

#calculate means and standard dev for each timepoint
means1=apply(c1, 2, mean)
means1=as.data.frame(means1)

dev1=apply(c1, 2, sd)
dev1=as.data.frame(dev1)

meansd2=apply(cd2, 2, mean)
meansd2=as.data.frame(meansd2)

devd2=apply(cd2, 2, sd)
devd2=as.data.frame(devd2)

#combine and calculate difference between means
all=cbind(means1, meansd2, dev1, devd2)
all$difference=(all$meansd2-all$means1)
all=all[order(all$devd2, -all$difference ),]

#only more abundant in diseased
allpos=all[all$difference>=0,]

#create index based on sd and diff
allpos$index=(allpos$dev1+allpos$devd2)-log(allpos$difference)
allpos=allpos[with(allpos, order(index)),]

plot(all$difference, all$devd2, ylab="Time2 diseased sd", xlab="difference", main="Difference of T2 d and T1 vs sd of T2 d")

#read tax info
t=read.csv("input_files//rep_set_tax_assignments.csv", header=T)

#do the saame thing with healthy
c1=one[14:NCOL(one)]
row.names(c1)=one$Sample_name

twoh=dnorm[dnorm$Timepoint=="two"& dnorm$Dose_disease_state=="Healthy",]
ch2=twoh[14:NCOL(twoh)]
row.names(ch2)=twoh$Sample_name
sh2=twoh[0:13]

combinedh=rbind(c1, ch2)

combinedh2=t(combinedh)
colnames(combinedh2)=rownames(combinedh)
rs = rowSums (combinedh2)
use = (rs > 0)
combinedh2 = combinedh2 [ use, ]


ch2=combinedh2[,83:NCOL(combinedh2)]
ch2=t(ch2)

#redo t1 with all otus
c1=t(combinedh2[,1:82])


means1=apply(c1, 2, mean)
means1=as.data.frame(means1)

dev1=apply(c1, 2, sd)
dev1=as.data.frame(dev1)

#calculate mean and sd for healthy
meansh2=apply(ch2, 2, mean)
meansh2=as.data.frame(meansh2)

devh2=apply(ch2, 2, sd)
devh2=as.data.frame(devh2)

all=cbind(means1, meansh2, dev1, devh2)
all$difference=(all$meansh2-all$means1)
all=all[order(all$devh2, all$difference ),]

#plot difference and sd
plot(all$difference, all$devh2, ylab="Time2 diseased sd", xlab="difference", main="Difference of T2 h and T1 vs sd of T2 h")

#determine those otus that are more abundant in healthy T2 than T1
filth=all[all$difference>0& all$devh2<200, ]

#remove those otus mor abun in healthy from the more abun in diseased
allposnoh=subset(allpos, !(row.names(allpos) %in% row.names(filth)))

allposnoht=merge(allposnoh, t, by.x="row.names", by.y="OTUID")

allposnoht=allposnoht[order(allposnoht$index ),]

allposnoht$rank=1:NROW(allposnoht)

#read in the phylogenetic tree of all 20000 otus
tree=read.tree(file="input_files/filtered_20000_tree.tre")

dropotus=subset(filtered$OTUID, !(filtered$OTUID %in% row.names(allposnoh)))

dropotus=as.vector(dropotus)

#prune tree so only otus in allposnoh are in tree
apnhtree=drop.tip(tree, dropotus)

#calculate pairwise phylogenetic distance
dm=cophenetic.phylo(apnhtree)

#calculate minimum phylogenetic distance for each otu
minphy=apply(dm, 1, FUN=function(x) {min(x[x>0])})
minphy=as.data.frame(minphy)

row.names(minphy)=row.names(dm)

allposnohtphy=merge(allposnoht, minphy, by.x="Row.names", by.y="row.names")

#plot minphy by index
plot(allposnohtphy$minphy, allposnohtphy$rank)

allposnohtphy200=allposnohtphy[allposnohtphy$minphy<0.1&allposnohtphy$rank<200,]

dropotus=subset(filtered$OTUID, !(filtered$OTUID %in% row.names(allposnoh)))

filtphyotu200=subset(combined, (row.names(combined) %in% allposnohtphy200$Row.names))

write.csv(filtphyotu200, file="Output_files/filtphyotutable200.csv")

#create full otutable

threed=dnorm[dnorm$Timepoint=="three"& dnorm$T3_disease_state=="Diseased",]
cd3=threed[14:NCOL(threed)]
row.names(cd3)=threed$Sample_name

cd3=t(cd3)

filtphyotuthree200=subset(cd3, (row.names(cd3) %in% allposnohtphy200$Row.names))

filtphyotuall200=cbind(filtphyotuthree200, filtphyotu200)
write.csv(filtphyotuall200, file="Output_files/filtphyotutableall200.csv")

"New.ReferenceOTU21033", "New.ReferenceOTU25985"39471, "2373263"39460, "New.ReferenceOTU28484"19986, "New.ReferenceOTU26314"19986
"New.ReferenceOTU14529"32868, "New.ReferenceOTU16534"39482
getMRCA(tree, c( "New.ReferenceOTU21033","New.ReferenceOTU25985", "2373263", "New.ReferenceOTU16534"))

getMRCA(tree, c( "New.ReferenceOTU25985", "New.ReferenceOTU16534"))
24627

pasttree=extract.clade(tree, 39471)
plot.phylo(pasttree, show.tip.label=TRUE, cex=0.5)
tiplabels(text=c("21033", ""), tip=c(20, 138))

pastclade=pasttree$tip.label
pastclade=as.data.frame(pastclade)
pastcladet=merge(pastclade, t, by.x="pastclade", by.y="OTUID")

write.tree(pasttree, file="Output_files/pasttree.tre")

pastcladeotu3=subset(cd3, (row.names(cd3) %in% pastclade$pastclade))
pastcladeotu2=subset(combined, (row.names(combined) %in% pastclade$pastclade))
pastcladeallotu=cbind(pastcladeotu2, pastcladeotu3)

dnorm=dnorm[order(dnorm$T3_disease_state, dnorm$Timepoint ),]
otu=dnorm[,14:NCOL(dnorm)]
row.names(otu)=dnorm$Sample_name
otu=t(otu)
pastcladeotu=subset(combined, (row.names(combined) %in% pastclade$pastclade))
write.csv(pastcladeotu, file="Output_files/pastcladeotu12_s.csv")


##redo with silva tree 
trees=read.tree(file="input_files/rep_set_filtered20000_silva_phylogeny.tre")
  19986, "New.ReferenceOTU26314"19986
32868
getMRCA(trees, c( "New.ReferenceOTU21033","New.ReferenceOTU25985", "2373263", 
                  "New.ReferenceOTU16534", "New.ReferenceOTU28484", "New.ReferenceOTU26314",
                  "New.ReferenceOTU14529"))

getMRCA(trees, c(  "2373263", "New.ReferenceOTU14529", "New.ReferenceOTU16534" ))

pasttree=extract.clade(trees, 31178)
plot.phylo(pasttree, show.tip.label=TRUE, cex=0.5)
tiplabels(text=c("o", "o", "o", "o", "o", "o", "o"), tip=c(3, 14, 47, 51, 119, 181, 186), cex=0.5)

pastclade=pasttree$tip.label
pastclade=as.data.frame(pastclade)

ts=read.csv("input_files/rep_set_tax_assignments_silva2.csv", header=F)
colnames(ts)=c("OTUID", "kingdom", "phylum", "class", "order", "family", "genus", "species", "evalue", "ref")

pastcladet=merge(pastclade, ts, by.x="pastclade", by.y="OTUID")

write.tree(pasttree, file="Output_files/pasttree_s.tre")
