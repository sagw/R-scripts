rm(list=ls())
require(vegan)
require(DESeq2)
require(MASS)
require(ggplot2)
require(plyr)
setwd("~/Documents/cbbs4_te/cbbs4_te/")

c=read.csv("cbbs4_te_allotus.csv", h=T, row.names=1)

s=read.csv("cbbs4_te_data.csv", h=T, as.is=T)

b=t(c)

d=cbind.data.frame(s,b)

check=cbind(row.names(d), d[1:1])

d=d[d$Disease_state!="?",]

d=d[d$base_mid!="M",]

c=d[9:NCOL(d)]

s=d[0:8]

# Remove first column
d=d[2:(NCOL(d))]

c=t(c)

dds <-DESeqDataSetFromMatrix(countData = c, colData=s, 
                             design= ~ resistant_susceptible)

dds=estimateSizeFactors(dds)

comm.normalized=(counts(dds, normalized=TRUE))

t=read.csv("rep_set_tax_assignments3.csv", h=T)

c1=c

rs = rowSums (c1)
use = (rs > 0)
c1 = c1 [ use, ]

ct=merge(c1, t, by.x="row.names", by.y="OTUID")
row.names(ct)=ct$Row.names
ctendo=ct[ct$family=="Endozoicimonaceae",]
ctendoc=ctendo[,2:55]
ctendoc$family="E"
ctendoagg=aggregate(. ~ family, data=ctendoc, FUN=sum)

endofull=t(ctendoagg)

endofull=endofull[-1,]
endofull=as.data.frame(endofull)
endofull=apply(endofull, 1, as.numeric)
endofull=as.matrix(endofull)
colnames(endofull)="Abundance"
endofulls=cbind(s, endofull)

hist(endofulls$Abundance)

#glm nb
endoglmnb=glm.nb(Abundance ~ resistant_susceptible*(Exposure/Disease_state), data=endofulls)
endoglmsum=summary(endoglmnb)
pvalues=endoglmsum$coefficients
pvalues=signif(pvalues, digits=3)
write.csv(pvalues, file="endoglmnb.csv")

#glm quasipoisson
endoglm=glm(Abundance ~ Disease_state*Day*(resistant_susceptible), data=endofulls, 
            family="quasipoisson")
endoglmsum=summary(endoglm)


  
#barplots with mean and se
meanse=ddply(endofulls, c("resistant_susceptible", "Disease_state", "Exposure"), 
             summarise, mean=mean(Abundance), sem=sd(Abundance)/sqrt(length(Abundance)))
meanse= transform(meanse, lower=mean-sem, upper=mean+sem)

eplot=qplot(x=resistant_susceptible , y=mean, fill=resistant_susceptible,
data=meanse, geom="bar", stat="identity",
position="dodge", facets=Exposure~Disease_state)


eplot+ geom_errorbar(aes(ymax=upper,
                  ymin=lower),
              position=position_dodge(0.9),
              data=meanse) +theme_bw()+ scale_fill_grey(start = 0, end = .9)


