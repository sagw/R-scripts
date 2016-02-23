rm(list=ls())
require(vegan)
require(DESeq2)
require(genefilter)
require(GUniFrac)
setwd("~/Documents/SD1/SD1_analyses/")
c=read.csv("input_files/filtered_otutable.csv", h=T, row.names=1)
c=c[7:NCOL(c)]

#information about samples
s=read.csv("input_files//SD1_map.csv", h=T, as.is=T)
s=s[!s$Timepoint=="control",]
b=t(c)
d=cbind.data.frame(s,b)
row.names(d)=row.names(b)

d=d[!d$Sample_name=="SD1.CK14D2.CK14.O",]
d=d[!d$Sample_name=="SD1.CK14.D2.CK4.W.T57.D",]

d$Dose_disease_state=ifelse(d$Timepoint=="one", "Healthy", d$Dose_disease_state)

#check order of samples
check=cbind(row.names(d), d[1:1])

#remove past data
d[13]=NULL

experimental=d[d$Timepoint %in% c("one", "two", "three"),]

d1=experimental
c1=d1[13:NCOL(d1)]
s1=d1[0:12]

c1=t(c1)
rs = rowSums (c1)
use = (rs > 0)
c1 = c1 [ use, ]

dds <-DESeqDataSetFromMatrix(countData = c1, colData=s1, design= ~ Genotype)

#calculate arithemtic means because of zeros.
arith_mean=function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0])))
arithmeans = apply(c1, 1, arith_mean) 

#estimate size factors
dds = estimateSizeFactors(dds, geoMeans=arithmeans,locfunc=shorth)

comm.normalized=t(counts(dds, normalized=TRUE))

d2=cbind.data.frame(s1,comm.normalized)

#check order of samples
check=cbind(row.names(d2), d2[1:1])

one=d2[d2$Timepoint %in% c("one"),]

d3=one
c2=d3[13:NCOL(d3)]
s2=d3[0:12]

c2=t(c2)
rs = rowSums (c2)
use = (rs > 0)
c2norm = c2 [ use, ]
c2norm=t(c2norm)
         
perm=adonis(c2norm ~ Genotype, data=s2, permutations=999, method="bray")

twoThree=d2[d2$Timepoint %in% c("two", "three"),]

d4=twoThree
c3=d4[13:NCOL(d4)]
s3=d4[0:12]

c3=t(c3)
rs = rowSums (c3)
use = (rs > 0)
c3norm = c3 [ use, ]
c3norm=t(c2norm)
s3=s1

perm23=adonis(c3norm ~ T3_disease_state + Site*Dose_disease_state*Timepoint*Dose_Site, data=s3, permutations=999, method="bray")

