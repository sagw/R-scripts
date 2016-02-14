rm(list=ls())
require(vegan)
require(DESeq2)
require(genefilter)
require(car)
require(MASS)
require(lme4)


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

experimental=d[d$Timepoint %in% c("dose", "one", "two", "three"),]

d1=experimental
c1=d1[13:NCOL(d1)]
s1=d1[0:12]

c1=t(c1)
rs = rowSums (c1)
use = (rs > 0)
c1 = c1 [ use, ]

#create dds
dds =DESeqDataSetFromMatrix(countData= c1, colData=s1, design= ~Genotype)

#calculate arithemtic means because of zeros.
arith_mean=function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0])))
arithmeans = apply(c1, 1, arith_mean) 

#estimate size factors
dds = estimateSizeFactors(dds, geoMeans=arithmeans,locfunc=shorth)

#get normalized counts from DESEQ
c2normfull=counts(dds, normalized=TRUE)
write.csv(c2normfull, file="input_files/c2norm.csv")
d2=cbind.data.frame(s1,c2normfull)

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

#pseudorsquare function
pseudo.rsq <- function(m) {
  r2 <- cor.test(fitted(m), m$data$vec)
  return(as.numeric(r2$estimate))
}

#make empty dataframe for pvalues
allpvalues_g =matrix(NA,nrow=1,ncol=nrow(c2norm))
colnames(allpvalues_g)=row.names(c2norm)

allrsquare_g =matrix(NA,nrow=1,ncol=nrow(c2norm))
colnames(allrsquare_g)=row.names(c2norm)

##run glmm and Anova
for (i in 1:nrow(c2norm)) {
  tryCatch({
    vec = as.numeric(c2norm[i,])
    combined=cbind.data.frame(s2, vec)
    glmtest2=glmmPQL(vec ~ Genotype,
                     random = ~1|Dose, data=combined, family = "quasipoisson", verbose=TRUE)
    anovaglm=as.data.frame(Anova(glmtest2, test="Chisq", error.estimate="dispersion"))
    colnames(anovaglm)=c("chisq", "df", "pvalue")
    rsquare=pseudo.rsq(glmtest2)
    allrsquare_g[,i]=rsquare
    anovaglm$effect=row.names(anovaglm)
    pvalues=anovaglm$pvalue
    pvalues=as.data.frame(pvalues)
    pvalues=apply(pvalues, 2, as.numeric)
    pvalues=as.data.frame(pvalues)
    row.names(pvalues)=row.names(anovaglm)
    pvalues=t(pvalues)
    Genotype_p=pvalues[,"Genotype"]
    allpvalues_g[,i] = c(Genotype_p)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

rownames(allpvalues_g) = c("Genotype_p")

allrsquare_g=as.data.frame(allrsquare_g)

#look at rsquare values
allrsquare_g=t(allrsquare_g)
allrsquare_g_high=subset(allrsquare_g, allrsquare_g>0.7)

#separate p values
allpvalues_g=t(allpvalues_g)
allpvalues_g=apply(allpvalues_g, 2, as.numeric)
rownames(allpvalues_g)=rownames(allrsquare_g)

tgg=read.csv("input_files/rep_set_tax_assignments_split2.csv", header=T)
ts=read.csv("input_files/rep_set_tax_assignments_silva_clean2.csv", header=F)
colnames(ts)=c("OTUID", "kingdom_s", "phylum_s", "class_s", "order_s", "family_s", "genus_s", "species_s", "evalue", "ref")

#pull out t3 disease sig otus
Genotype= as.matrix(allpvalues_g[,1])
Genotype_padj=as.matrix(p.adjust(Genotype, method="fdr"))
colnames(Genotype_padj)=c("padj")
row.names(Genotype_padj)=row.names(allpvalues_g)
Genotype_sigOTUs=subset(Genotype_padj, Genotype_padj<0.05)
Genotype_sigOTUs=as.data.frame(Genotype_sigOTUs)
Genotype_sigOTUs$OTUID=row.names(Genotype_sigOTUs)
Genotype_sigOTUs_tax=merge(Genotype_sigOTUs, ts, by="OTUID")
Genotype_sigOTUs_alltax=merge(Genotype_sigOTUs_tax, tgg, by="OTUID")
Genotype_sigOTUs_alltaxr=merge(Genotype_sigOTUs_alltax, allrsquare_g, by.x="OTUID", by.y="row.names")
colnames(Genotype_sigOTUs_alltaxr)[19]=c("pseudor2")
write.csv(Genotype_sigOTUs_alltaxr, file="Output_files/GLMM_results/Genotype_sigOTUs_alltax.csv")
