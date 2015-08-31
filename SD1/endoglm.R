rm(list=ls())
require(vegan)
require(DESeq2)
require(MASS)
require(ggplot2)
require(plyr)
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
#check order of samples
check=cbind(row.names(d), d[1:1])

t=read.csv("input_files//rep_set_tax_assignments_split2.csv", header=T)

t=read.csv("input_files//rep_set_tax_assignments_silva2.csv")
colnames(t)=c("OTUID", "kingdom", "phylum", "class", "order", "family", "genus", "species", "evalue", "ref")


d[14]=NULL
d[13]=NULL

experimental=d[d$Timepoint %in% c("one", "two", "three"),]

d1=experimental
c1=d1[13:NCOL(d1)]
s1=d1[0:12]

c1=t(c1)
rs = rowSums (c1)
use = (rs > 0)
c1 = c1 [ use, ]

dds <-DESeqDataSetFromMatrix(countData = c1, colData=s1, 
                             design= ~ T3_disease_state)

#calculate different geometric means because of zeros.
gm_mean=function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0])))
geoMeans = apply(c1, 1, gm_mean) 

#use locfunc=shorth instead of median because of low counts.
dds = estimateSizeFactors(dds, geoMeans=geoMeans, locfunc=shorth)

comm.normalized=(counts(dds, normalized=TRUE))
colnames(comm.normalized)=colnames(c1)

c2=comm.normalized

rs = rowSums (c2)
use = (rs > 0)
c2 = c2 [ use, ]

ct=merge(t, c2,  by.x="OTUID", by.y="row.names")
row.names(ct)=ct$Row.names
ctendo=ct[ct$genus=="Endozoicomonas",]
ctendoc=ctendo[,11:NCOL(ctendo)]
row.names(ctendoc)=ctendo$OTUID

ctendoc$genus="E"
ctendoagg=aggregate(. ~ genus, data=ctendoc, FUN=sum)

max(ctendoagg[,2:263])

endofull=t(ctendoagg)

endofull=endofull[-(1),]

endofull=as.data.frame(endofull)
endofull=apply(endofull, 1, as.numeric)
endofull=as.matrix(endofull)
colnames(endofull)="Abundance"
endofulls=cbind(s1, endofull)

hist(endofulls$Abundance)

is.numeric(endofulls$Abundance)

endofulls$Genotype=as.factor(endofulls$Genotype)

#glm quasipoisson
endoglm=glm(Abundance ~ T3_disease_state*Timepoint*Genotype, data=endofulls, 
            family="quasipoisson")
summary(endoglm)

genos=unique(endofulls$Genotype)
genos=as.data.frame(genos)

results = data.frame(Estimate= numeric(0), SE= numeric(0), z = numeric(0), pvalue=numeric(0))

for (i in 1:NROW(genos)) {
  level=print(as.character(genos[i,]))
  endofulls$Genotype=relevel(endofulls$Genotype, ref=level)
  endoglmnb=glm.nb(Abundance ~ (Timepoint/T3_disease_state)*Genotype, data=endofulls)
  endoglmsum=summary(endoglmnb)
  pvalues=endoglmsum$coefficients
  colnames(pvalues)[4]= "p"
  pvalues1=apply(pvalues, 2, as.numeric)
  pvalues1=apply(pvalues1, 2, signif, digits=2)
  row.names(pvalues1)=row.names(pvalues)
  pvalues1=as.data.frame(pvalues1)
  sigpvalues=subset(pvalues1, pvalues1$p<.05)
  sigpvalues=rbind(level, sigpvalues)
  results=rbind(results,sigpvalues, row.names=TRUE)
}


write.csv(results, file="endoglm_genotype_sig_nb.csv")


#glm nb
endoglmnb=glm.nb(Abundance ~ Timepoint/T3_disease_state, data=endofulls, contrasts=FALSE)
summary(endoglmnb)


endoglmsum=summary(endoglmnb)
pvalues=endoglmsum$coefficients
pvalues=signif(pvalues, digits=3)
write.csv(pvalues, file="endoglmnb.csv")

#barplots with mean and se
meanse=ddply(endofulls, c("T3_disease_state", "Timepoint", "Dose_disease_state"), 
             summarise, mean=mean(Abundance), sem=sd(Abundance)/sqrt(length(Abundance)))
meanse= transform(meanse, lower=mean-sem, upper=mean+sem)

eplot=qplot(x=Timepoint , y=mean, fill=T3_disease_state,
            data=meanse, geom="bar", stat="identity",
            position="dodge", facets=.~Dose_disease_state)


eplot+ geom_errorbar(aes(ymax=upper,
                         ymin=lower),
                     position=position_dodge(0.9),
                     data=meanse) +theme_bw()

#+ scale_fill_grey(start = 0, end = .9)



