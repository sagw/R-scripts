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

experimental=d[d$Timepoint %in% c("two", "three"),]

d1=experimental
c1=d1[13:NCOL(d1)]
s1=d1[0:12]

c1=t(c1)
rs = rowSums (c1)
use = (rs > 0)
c1 = c1 [ use, ]

#create dds
dds =DESeqDataSetFromMatrix(countData= c1, colData=s1, design= ~T3_disease_state)
colData(dds)$T3_disease_state= factor(colData(dds)$T3_disease_state, levels=c("Healthy","Diseased"))

#calculate arithemtic means because of zeros.
arith_mean=function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0])))
arithmeans = apply(c1, 1, arith_mean) 

#estimate size factors
dds = estimateSizeFactors(dds, geoMeans=arithmeans, locfunc=shorth)

#get normalized counts from DESEQ
c2norm=counts(dds, normalized=TRUE)

#remove otus found fewer than 10 times
c2normfiltered=c2norm[apply(c2norm, MARGIN=1, function(x) (sum(x>0))>=3),]

c2normfilteredsub=head(c2normfiltered, n=100)

  
pseudo.rsq <- function(m) {
  r2 <- cor.test(fitted(m), m$data$vec)
  return(as.numeric(r2$estimate))
}

#make empty dataframe for pvalues
allpvalues_qp2 =matrix(NA,nrow=17,ncol=nrow(c2norm))
colnames(allpvalues_qp2)=row.names(c2norm)

allrsquare_qp2 =matrix(NA,nrow=1,ncol=nrow(c2norm))
colnames(allrsquare_qp2)=row.names(c2norm)
)
##test with three otus
#full model
for (i in 69964:nrow(c2norm)) {
  tryCatch({
  vec = as.numeric(c2norm[i,])
  combined=cbind.data.frame(s1, vec)
  glmtest2=glmmPQL(vec ~ T3_disease_state*Dose_disease_state + Site*Dose_disease_state*Timepoint*Dose_Site,
                   random = ~1|Dose/Timepoint, data=combined, family = "quasipoisson")
  anovaglm=as.data.frame(Anova(glmtest2, test="Chisq", error.estimate="dispersion"))
  colnames(anovaglm)=c("chisq", "df", "pvalue")
  rsquare=pseudo.rsq(glmtest2)
  allrsquare_qp[,i]=rsquare
  anovaglm$effect=row.names(anovaglm)
  pvalues=anovaglm$pvalue
  pvalues=as.data.frame(pvalues)
  pvalues=apply(pvalues, 2, as.numeric)
  pvalues=as.data.frame(pvalues)
  row.names(pvalues)=row.names(anovaglm)
  pvalues=t(pvalues)
  T3_disease_state_p=pvalues[,"T3_disease_state"]
  Dose_disease_state_p=pvalues[,"Dose_disease_state"]
  Site_p=pvalues[,"Site"]
  Timepoint_p=pvalues[,"Timepoint"]
  Dose_Site_p=pvalues[,"Dose_Site"]
  T3_disease_state_Dose_disease_state_p=pvalues[,"T3_disease_state:Dose_disease_state"]
  Dose_disease_state_Site_p=pvalues[,"Dose_disease_state:Site"]
  Site_Timepoint_p=pvalues[,"Site:Timepoint"]
  Dose_disease_state_Timepoint_p=pvalues[,"Dose_disease_state:Timepoint"]
  Site_Dose_Site_p=pvalues[,"Site:Dose_Site"]
  Dose_disease_state_Dose_Site_p=pvalues[,"Dose_disease_state:Dose_Site"]
  Timepoint_Dose_Site_p=pvalues[,"Timepoint:Dose_Site"]
  Dose_disease_state_Site_Timepoint_p=pvalues[,"Dose_disease_state:Site:Timepoint"]
  Dose_disease_state_Site_Dose_Site_p=pvalues[,"Dose_disease_state:Site:Dose_Site"]
  Site_Timepoint_Dose_Site_p=pvalues[,"Site:Timepoint:Dose_Site"]
  Dose_disease_state_Timepoint_Dose_Site_p=pvalues[,"Dose_disease_state:Timepoint:Dose_Site"]
  Dose_disease_state_Site_Timepoint_Dose_Site_p=pvalues[,"Dose_disease_state:Site:Timepoint:Dose_Site"]
  
  allpvalues_qp2[,i] = c(T3_disease_state_p, Dose_disease_state_p, Site_p,Timepoint_p,
                        Dose_Site_p, T3_disease_state_Dose_disease_state_p,  Dose_disease_state_Site_p,
                        Site_Timepoint_p, Dose_disease_state_Timepoint_p, Site_Dose_Site_p,
                        Dose_disease_state_Dose_Site_p, Timepoint_Dose_Site_p, Dose_disease_state_Site_Timepoint_p,
                        Dose_disease_state_Site_Dose_Site_p, Site_Timepoint_Dose_Site_p, 
                        Dose_disease_state_Timepoint_Dose_Site_p, Dose_disease_state_Site_Timepoint_Dose_Site_p)
  
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

rownames(allpvalues_qp2) = c("T3_disease_state_p", "Dose_disease_state_p", "Site_p","Timepoint_p",
                            "Dose_Site_p", "T3_disease_state_Dose_disease_state_p",  "Dose_disease_state_Site_p",
                            "Site_Timepoint_p", "Dose_disease_state_Timepoint_p", "Site_Dose_Site_p",
                            "Dose_disease_state_Dose_Site_p", "Timepoint_Dose_Site_p", "Dose_disease_state_Site_Timepoint_p",
                            "Dose_disease_state_Site_Dose_Site_p", "Site_Timepoint_Dose_Site_p", 
                            "Dose_disease_state_Timepoint_Dose_Site_p", "Dose_disease_state_Site_Timepoint_Dose_Site_p")

allpvalues_complete=cbind(allpvalues_qp[,1:69964], allpvalues_qp2[,69965:NCOL(allpvalues_qp2)])
allrsquare_qp=as.data.frame(allrsquare_qp)
allrsquare_qp2=as.data.frame(allrsquare_qp2)
allrsquare_complete=cbind(allrsquare_qp[,1:69964], allrsquare_qp2[,69965:NCOL(allrsquare_qp2)])

#look at rsquare values
allrsquare_complete=t(allrsquare_complete)
allrsquare_complete_high=subset(allrsquare_complete, allrsquare_complete>0.1)

#separate p values
allpvalues_complete=t(allpvalues_complete)
allpvalues_complete=apply(allpvalues_complete, 2, as.numeric)

tgg=read.csv("input_files/rep_set_tax_assignments_split2.csv", header=T)
ts=read.csv("input_files/rep_set_tax_assignments_silva_clean2.csv", header=F)
colnames(ts)=c("OTUID", "kingdom_s", "phylum_s", "class_s", "order_s", "family_s", "genus_s", "species_s", "evalue", "ref")

#pull out t3 disease sig otus
T3_disease_state= as.matrix(allpvalues_complete[,1])
T3_disease_state_padj=as.matrix(p.adjust(T3_disease_state, method="fdr"))
colnames(T3_disease_state_padj)=c("padj")
row.names(T3_disease_state_padj)=row.names(allpvalues_complete)
T3_disease_state_sigOTUs=subset(T3_disease_state_padj, T3_disease_state_padj<0.05)
T3_disease_state_sigOTUs=as.data.frame(T3_disease_state_sigOTUs)
T3_disease_state_sigOTUs$OTUID=row.names(T3_disease_state_sigOTUs)
T3_disease_state_sigOTUs_tax=merge(T3_disease_state_sigOTUs, ts, by="OTUID")
T3_disease_state_sigOTUs_alltax=merge(T3_disease_state_sigOTUs_tax, tgg, by="OTUID")
T3_disease_state_sigOTUs_alltaxr=merge(T3_disease_state_sigOTUs_alltax, allrsquare_complete, by.x="OTUID", by.y="row.names")
colnames(T3_disease_state_sigOTUs_alltaxr)[19]=c("pseudor2")
write.csv(T3_disease_state_sigOTUs_alltaxr, file="Output_files/T3_disease_state_sigOTUs_alltax.csv")

#Dose disease state
Dose_disease_state= as.matrix(allpvalues_complete[,2])
Dose_disease_state_padj=as.matrix(p.adjust(Dose_disease_state, method="fdr"))
colnames(Dose_disease_state_padj)=c("padj")
row.names(Dose_disease_state_padj)=row.names(allpvalues_complete)
Dose_disease_state_sigOTUs=subset(Dose_disease_state_padj, Dose_disease_state_padj<0.05)
Dose_disease_state_sigOTUs=as.data.frame(Dose_disease_state_sigOTUs)
Dose_disease_state_sigOTUs$OTUID=row.names(Dose_disease_state_sigOTUs)
Dose_disease_state_sigOTUs_tax=merge(Dose_disease_state_sigOTUs, ts, by="OTUID")
Dose_disease_state_sigOTUs_alltax=merge(Dose_disease_state_sigOTUs_tax, tgg, by="OTUID")
Dose_disease_state_sigOTUs_alltaxr=merge(Dose_disease_state_sigOTUs_alltax, allrsquare_complete, by.x="OTUID", by.y="row.names")
colnames(Dose_disease_state_sigOTUs_alltaxr)[19]=c("pseudor2")
write.csv(Dose_disease_state_sigOTUs_alltax, file="Output_files/Dose_disease_state_sigOTUs_alltax.csv")

#Site
Site= as.matrix(allpvalues_complete[,3])
Site_padj=as.matrix(p.adjust(Site, method="fdr"))
colnames(Site_padj)=c("padj")
row.names(Site_padj)=row.names(allpvalues_complete)
Site_sigOTUs=subset(Site_padj, Site_padj<0.05)
Site_sigOTUs=as.data.frame(Site_sigOTUs)
Site_sigOTUs$OTUID=row.names(Site_sigOTUs)
Site_sigOTUs_tax=merge(Site_sigOTUs, ts, by="OTUID")
Site_sigOTUs_alltax=merge(Site_sigOTUs_tax, tgg, by="OTUID")
Site_sigOTUs_alltaxr=merge(Site_sigOTUs_alltax, allrsquare_complete, by.x="OTUID", by.y="row.names")
colnames(Site_sigOTUs_alltaxr)[19]=c("pseudor2")
write.csv(Site_sigOTUs_alltax, file="Output_files/Site_sigOTUs_alltax.csv")

#Timepoint
Timepoint= as.matrix(allpvalues_complete[,4])
Timepoint_padj=as.matrix(p.adjust(Timepoint, method="fdr"))
colnames(Timepoint_padj)=c("padj")
row.names(Timepoint_padj)=row.names(allpvalues_complete)
Timepoint_sigOTUs=subset(Timepoint_padj, Timepoint_padj<0.05)
Timepoint_sigOTUs=as.data.frame(Timepoint_sigOTUs)
Timepoint_sigOTUs$OTUID=row.names(Timepoint_sigOTUs)
Timepoint_sigOTUs_tax=merge(Timepoint_sigOTUs, ts, by="OTUID")
Timepoint_sigOTUs_alltax=merge(Timepoint_sigOTUs_tax, tgg, by="OTUID")
T3_disease_state_sigOTUs_alltaxr=merge(T3_disease_state_sigOTUs_alltax, allrsquare_complete, by.x="OTUID", by.y="row.names")
colnames(T3_disease_state_sigOTUs_alltaxr)[19]=c("pseudor2")
write.csv(Timepoint_sigOTUs_alltax, file="Output_files/Timepoint_sigOTUs_alltax.csv")

#Dose_Site
Dose_Site= as.matrix(allpvalues_complete[,5])
Dose_Site_padj=as.matrix(p.adjust(Dose_Site, method="fdr"))
colnames(Dose_Site_padj)=c("padj")
row.names(Dose_Site_padj)=row.names(allpvalues_complete)
Dose_Site_sigOTUs=subset(Dose_Site_padj, Dose_Site_padj<0.05)
Dose_Site_sigOTUs=as.data.frame(Dose_Site_sigOTUs)
Dose_Site_sigOTUs$OTUID=row.names(Dose_Site_sigOTUs)
Dose_Site_sigOTUs_tax=merge(Dose_Site_sigOTUs, ts, by="OTUID")
Dose_Site_sigOTUs_alltax=merge(Dose_Site_sigOTUs_tax, tgg, by="OTUID")
Dose_Site_sigOTUs_alltaxr=merge(Dose_Site_sigOTUs_alltax, allrsquare_complete, by.x="OTUID", by.y="row.names")
colnames(Dose_Site_sigOTUs_alltaxr)[19]=c("pseudor2")
write.csv(Dose_Site_sigOTUs_alltax, file="Output_files/Dose_Site_sigOTUs_alltax.csv")

#T3_disease_state_Dose_disease_state
T3_disease_state_Dose_disease_state= as.matrix(allpvalues_complete[,6])
T3_disease_state_Dose_disease_state_padj=as.matrix(p.adjust(T3_disease_state_Dose_disease_state, method="fdr"))
colnames(T3_disease_state_Dose_disease_state_padj)=c("padj")
row.names(T3_disease_state_Dose_disease_state_padj)=row.names(allpvalues_complete)
T3_disease_state_Dose_disease_state_sigOTUs=subset(T3_disease_state_Dose_disease_state_padj, T3_disease_state_Dose_disease_state_padj<0.05)
T3_disease_state_Dose_disease_state_sigOTUs=as.data.frame(T3_disease_state_Dose_disease_state_sigOTUs)
T3_disease_state_Dose_disease_state_sigOTUs$OTUID=row.names(T3_disease_state_Dose_disease_state_sigOTUs)
T3_disease_state_Dose_disease_state_sigOTUs_tax=merge(T3_disease_state_Dose_disease_state_sigOTUs, ts, by="OTUID")
T3_disease_state_Dose_disease_state_sigOTUs_alltax=merge(T3_disease_state_Dose_disease_state_sigOTUs_tax, tgg, by="OTUID")
T3_disease_state_Dose_disease_state_sigOTUs_alltaxr=merge(T3_disease_state_Dose_disease_state_sigOTUs_alltax, allrsquare_complete, by.x="OTUID", by.y="row.names")
colnames(T3_disease_state_Dose_disease_state_sigOTUs_alltaxr)[19]=c("pseudor2")
write.csv(T3_disease_state_Dose_disease_state_sigOTUs_alltax, file="Output_files/T3_disease_state_Dose_disease_state_sigOTUs_alltax.csv")

#Dose_disease_state_Site
Dose_disease_state_Site= as.matrix(allpvalues_complete[,7])
Dose_disease_state_Site_padj=as.matrix(p.adjust(Dose_disease_state_Site, method="fdr"))
colnames(Dose_disease_state_Site_padj)=c("padj")
row.names(Dose_disease_state_Site_padj)=row.names(allpvalues_complete)
Dose_disease_state_Site_sigOTUs=subset(Dose_disease_state_Site_padj, Dose_disease_state_Site_padj<0.05)
Dose_disease_state_Site_sigOTUs=as.data.frame(Dose_disease_state_Site_sigOTUs)
Dose_disease_state_Site_sigOTUs$OTUID=row.names(Dose_disease_state_Site_sigOTUs)
Dose_disease_state_Site_sigOTUs_tax=merge(Dose_disease_state_Site_sigOTUs, ts, by="OTUID")
Dose_disease_state_Site_sigOTUs_alltax=merge(Dose_disease_state_Site_sigOTUs_tax, tgg, by="OTUID")
Dose_disease_state_Site_sigOTUs_alltaxr=merge(Dose_disease_state_Site_sigOTUs_alltax, allrsquare_complete, by.x="OTUID", by.y="row.names")
colnames(Dose_disease_state_Site_sigOTUs_alltaxr)[19]=c("pseudor2")
write.csv(Dose_disease_state_Site_sigOTUs_alltax, file="Output_files/Dose_disease_state_Site_sigOTUs_alltax.csv")

#Dose_disease_state_Site
Dose_disease_state_Site= as.matrix(allpvalues_complete[,8])
Dose_disease_state_Site_padj=as.matrix(p.adjust(Dose_disease_state_Site, method="fdr"))
colnames(Dose_disease_state_Site_padj)=c("padj")
row.names(Dose_disease_state_Site_padj)=row.names(allpvalues_complete)
Dose_disease_state_Site_sigOTUs=subset(Dose_disease_state_Site_padj, Dose_disease_state_Site_padj<0.05)
Dose_disease_state_Site_sigOTUs=as.data.frame(Dose_disease_state_Site_sigOTUs)
Dose_disease_state_Site_sigOTUs$OTUID=row.names(Dose_disease_state_Site_sigOTUs)
Dose_disease_state_Site_sigOTUs_tax=merge(Dose_disease_state_Site_sigOTUs, ts, by="OTUID")
Dose_disease_state_Site_sigOTUs_alltax=merge(Dose_disease_state_Site_sigOTUs_tax, tgg, by="OTUID")
Dose_disease_state_Site_sigOTUs_alltaxr=merge(Dose_disease_state_Site_sigOTUs_alltax, allrsquare_complete, by.x="OTUID", by.y="row.names")
colnames(Dose_disease_state_Site_sigOTUs_alltaxr)[19]=c("pseudor2")
write.csv(Dose_disease_state_Site_sigOTUs_alltax, file="Output_files/Dose_disease_state_Site_sigOTUs_alltax.csv")


#Dose_disease_state_Timepoint
Dose_disease_state_Timepoint= as.matrix(allpvalues_complete[,9])
Dose_disease_state_Timepoint_padj=as.matrix(p.adjust(Dose_disease_state_Timepoint, method="fdr"))
colnames(Dose_disease_state_Timepoint_padj)=c("padj")
row.names(Dose_disease_state_Timepoint_padj)=row.names(allpvalues_complete)
Dose_disease_state_Timepoint_sigOTUs=subset(Dose_disease_state_Timepoint_padj, Dose_disease_state_Timepoint_padj<0.05)
Dose_disease_state_Timepoint_sigOTUs=as.data.frame(Dose_disease_state_Timepoint_sigOTUs)
Dose_disease_state_Timepoint_sigOTUs$OTUID=row.names(Dose_disease_state_Timepoint_sigOTUs)
Dose_disease_state_Timepoint_sigOTUs_tax=merge(Dose_disease_state_Timepoint_sigOTUs, ts, by="OTUID")
Dose_disease_state_Timepoint_sigOTUs_alltax=merge(Dose_disease_state_Timepoint_sigOTUs_tax, tgg, by="OTUID")
Dose_disease_state_Timepoint_sigOTUs_alltaxr=merge(Dose_disease_state_Timepoint_sigOTUs_alltax, allrsquare_complete, by.x="OTUID", by.y="row.names")
colnames(Dose_disease_state_Timepoint_sigOTUs_alltaxr)[19]=c("pseudor2")
write.csv(Dose_disease_state_Timepoint_sigOTUs_alltax, file="Output_files/Dose_disease_state_Timepoint_sigOTUs_alltax.csv")

#Site_Dose_Site
Site_Dose_Site= as.matrix(allpvalues_complete[,10])
Site_Dose_Site_padj=as.matrix(p.adjust(Site_Dose_Site, method="fdr"))
colnames(Site_Dose_Site_padj)=c("padj")
row.names(Site_Dose_Site_padj)=row.names(allpvalues_complete)
Site_Dose_Site_sigOTUs=subset(Site_Dose_Site_padj, Site_Dose_Site_padj<0.05)
Site_Dose_Site_sigOTUs=as.data.frame(Site_Dose_Site_sigOTUs)
Site_Dose_Site_sigOTUs$OTUID=row.names(Site_Dose_Site_sigOTUs)
Site_Dose_Site_sigOTUs_tax=merge(Site_Dose_Site_sigOTUs, ts, by="OTUID")
Site_Dose_Site_sigOTUs_alltax=merge(Site_Dose_Site_sigOTUs_tax, tgg, by="OTUID")
Site_Dose_Site_sigOTUs_alltaxr=merge(Site_Dose_Site_sigOTUs_alltax, allrsquare_complete, by.x="OTUID", by.y="row.names")
colnames(Site_Dose_Site_sigOTUs_alltaxr)[19]=c("pseudor2")
write.csv(Site_Dose_Site_sigOTUs_alltax, file="Output_files/Site_Dose_Site_sigOTUs_alltax.csv")

#Dose_disease_state_Dose_Site
Dose_disease_state_Dose_Site= as.matrix(allpvalues_complete[,11])
Dose_disease_state_Dose_Site_padj=as.matrix(p.adjust(Dose_disease_state_Dose_Site, method="fdr"))
colnames(Dose_disease_state_Dose_Site_padj)=c("padj")
row.names(Dose_disease_state_Dose_Site_padj)=row.names(allpvalues_complete)
Dose_disease_state_Dose_Site_sigOTUs=subset(Dose_disease_state_Dose_Site_padj, Dose_disease_state_Dose_Site_padj<0.05)
Dose_disease_state_Dose_Site_sigOTUs=as.data.frame(Dose_disease_state_Dose_Site_sigOTUs)
Dose_disease_state_Dose_Site_sigOTUs$OTUID=row.names(Dose_disease_state_Dose_Site_sigOTUs)
Dose_disease_state_Dose_Site_sigOTUs_tax=merge(Dose_disease_state_Dose_Site_sigOTUs, ts, by="OTUID")
Dose_disease_state_Dose_Site_sigOTUs_alltax=merge(Dose_disease_state_Dose_Site_sigOTUs_tax, tgg, by="OTUID")
Dose_disease_state_Dose_Site_sigOTUs_alltaxr=merge(Dose_disease_state_Dose_Site_sigOTUs_alltax, allrsquare_complete, by.x="OTUID", by.y="row.names")
colnames(Dose_disease_state_Dose_Site_sigOTUs_alltaxr)[19]=c("pseudor2")
write.csv(Dose_disease_state_Dose_Site_sigOTUs_alltax, file="Output_files/Dose_disease_state_Dose_Site_sigOTUs_alltax.csv")

#Timepoint_Dose_Site
Timepoint_Dose_Site= as.matrix(allpvalues_complete[,12])
Timepoint_Dose_Site_padj=as.matrix(p.adjust(Timepoint_Dose_Site, method="fdr"))
colnames(Timepoint_Dose_Site_padj)=c("padj")
row.names(Timepoint_Dose_Site_padj)=row.names(allpvalues_complete)
Timepoint_Dose_Site_sigOTUs=subset(Timepoint_Dose_Site_padj, Timepoint_Dose_Site_padj<0.05)
Timepoint_Dose_Site_sigOTUs=as.data.frame(Timepoint_Dose_Site_sigOTUs)
Timepoint_Dose_Site_sigOTUs$OTUID=row.names(Timepoint_Dose_Site_sigOTUs)
Timepoint_Dose_Site_sigOTUs_tax=merge(Timepoint_Dose_Site_sigOTUs, ts, by="OTUID")
Timepoint_Dose_Site_sigOTUs_alltax=merge(Timepoint_Dose_Site_sigOTUs_tax, tgg, by="OTUID")
Timepoint_Dose_Site_sigOTUs_alltaxr=merge(Timepoint_Dose_Site_sigOTUs_alltax, allrsquare_complete, by.x="OTUID", by.y="row.names")
colnames(Timepoint_Dose_Site_sigOTUs_alltaxr)[19]=c("pseudor2")
write.csv(Timepoint_Dose_Site_sigOTUs_alltax, file="Output_files/Timepoint_Dose_Site_sigOTUs_alltax.csv")



#Dose_disease_state_Site_Timepoint
Dose_disease_state_Site_Timepoint= as.matrix(allpvalues_complete[,13])
Dose_disease_state_Site_Timepoint_padj=as.matrix(p.adjust(Dose_disease_state_Site_Timepoint, method="fdr"))
colnames(Dose_disease_state_Site_Timepoint_padj)=c("padj")
row.names(Dose_disease_state_Site_Timepoint_padj)=row.names(allpvalues_complete)
Dose_disease_state_Site_Timepoint_sigOTUs=subset(Dose_disease_state_Site_Timepoint_padj, Dose_disease_state_Site_Timepoint_padj<0.05)
Dose_disease_state_Site_Timepoint_sigOTUs=as.data.frame(Dose_disease_state_Site_Timepoint_sigOTUs)
Dose_disease_state_Site_Timepoint_sigOTUs$OTUID=row.names(Dose_disease_state_Site_Timepoint_sigOTUs)
Dose_disease_state_Site_Timepoint_sigOTUs_tax=merge(Dose_disease_state_Site_Timepoint_sigOTUs, ts, by="OTUID")
Dose_disease_state_Site_Timepoint_sigOTUs_alltax=merge(Dose_disease_state_Site_Timepoint_sigOTUs_tax, tgg, by="OTUID")
Dose_disease_state_Site_Timepoint_sigOTUs_alltaxr=merge(Dose_disease_state_Site_Timepoint_sigOTUs_alltax, allrsquare_complete, by.x="OTUID", by.y="row.names")
colnames(Dose_disease_state_Site_Timepoint_sigOTUs_alltaxr)[19]=c("pseudor2")
write.csv(Dose_disease_state_Site_Timepoint_sigOTUs_alltax, file="Output_files/Dose_disease_state_Site_Timepoint_sigOTUs_alltax.csv")

#Dose_disease_state_Site_Dose_Site
Dose_disease_state_Site_Dose_Site= as.matrix(allpvalues_complete[,14])
Dose_disease_state_Site_Dose_Site_padj=as.matrix(p.adjust(Dose_disease_state_Site_Dose_Site, method="fdr"))
colnames(Dose_disease_state_Site_Dose_Site_padj)=c("padj")
row.names(Dose_disease_state_Site_Dose_Site_padj)=row.names(allpvalues_complete)
Dose_disease_state_Site_Dose_Site_sigOTUs=subset(Dose_disease_state_Site_Dose_Site_padj, Dose_disease_state_Site_Dose_Site_padj<0.05)
Dose_disease_state_Site_Dose_Site_sigOTUs=as.data.frame(Dose_disease_state_Site_Dose_Site_sigOTUs)
Dose_disease_state_Site_Dose_Site_sigOTUs$OTUID=row.names(Dose_disease_state_Site_Dose_Site_sigOTUs)
Dose_disease_state_Site_Dose_Site_sigOTUs_tax=merge(Dose_disease_state_Site_Dose_Site_sigOTUs, ts, by="OTUID")
Dose_disease_state_Site_Dose_Site_sigOTUs_alltax=merge(Dose_disease_state_Site_Dose_Site_sigOTUs_tax, tgg, by="OTUID")
Dose_disease_state_Site_Dose_Site_sigOTUs_alltaxr=merge(Dose_disease_state_Site_Dose_Site_sigOTUs_alltax, allrsquare_complete, by.x="OTUID", by.y="row.names")
colnames(Dose_disease_state_Site_Dose_Site_sigOTUs_alltaxr)[19]=c("pseudor2")
write.csv(Dose_disease_state_Site_Dose_Site_sigOTUs_alltax, file="Output_files/Dose_disease_state_Site_Dose_Site_sigOTUs_alltax.csv")

#Site_Timepoint_Dose_Site
Site_Timepoint_Dose_Site= as.matrix(allpvalues_complete[,15])
Site_Timepoint_Dose_Site_padj=as.matrix(p.adjust(Site_Timepoint_Dose_Site, method="fdr"))
colnames(Site_Timepoint_Dose_Site_padj)=c("padj")
row.names(Site_Timepoint_Dose_Site_padj)=row.names(allpvalues_complete)
Site_Timepoint_Dose_Site_sigOTUs=subset(Site_Timepoint_Dose_Site_padj, Site_Timepoint_Dose_Site_padj<0.05)
Site_Timepoint_Dose_Site_sigOTUs=as.data.frame(Site_Timepoint_Dose_Site_sigOTUs)
Site_Timepoint_Dose_Site_sigOTUs$OTUID=row.names(Site_Timepoint_Dose_Site_sigOTUs)
Site_Timepoint_Dose_Site_sigOTUs_tax=merge(Site_Timepoint_Dose_Site_sigOTUs, ts, by="OTUID")
Site_Timepoint_Dose_Site_sigOTUs_alltax=merge(Site_Timepoint_Dose_Site_sigOTUs_tax, tgg, by="OTUID")
Site_Timepoint_Dose_Site_sigOTUs_alltaxr=merge(Site_Timepoint_Dose_Site_sigOTUs_alltax, allrsquare_complete, by.x="OTUID", by.y="row.names")
colnames(Site_Timepoint_Dose_Site_sigOTUs_alltaxr)[19]=c("pseudor2")
write.csv(Site_Timepoint_Dose_Site_sigOTUs_alltax, file="Output_files/Site_Timepoint_Dose_Site_sigOTUs_alltax.csv")

#Dose_disease_state_Timepoint_Dose_Site
Dose_disease_state_Timepoint_Dose_Site= as.matrix(allpvalues_complete[,16])
Dose_disease_state_Timepoint_Dose_Site_padj=as.matrix(p.adjust(Dose_disease_state_Timepoint_Dose_Site, method="fdr"))
colnames(Dose_disease_state_Timepoint_Dose_Site_padj)=c("padj")
row.names(Dose_disease_state_Timepoint_Dose_Site_padj)=row.names(allpvalues_complete)
Dose_disease_state_Timepoint_Dose_Site_sigOTUs=subset(Dose_disease_state_Timepoint_Dose_Site_padj, Dose_disease_state_Timepoint_Dose_Site_padj<0.05)
Dose_disease_state_Timepoint_Dose_Site_sigOTUs=as.data.frame(Dose_disease_state_Timepoint_Dose_Site_sigOTUs)
Dose_disease_state_Timepoint_Dose_Site_sigOTUs$OTUID=row.names(Dose_disease_state_Timepoint_Dose_Site_sigOTUs)
Dose_disease_state_Timepoint_Dose_Site_sigOTUs_tax=merge(Dose_disease_state_Timepoint_Dose_Site_sigOTUs, ts, by="OTUID")
Dose_disease_state_Timepoint_Dose_Site_sigOTUs_alltax=merge(Dose_disease_state_Timepoint_Dose_Site_sigOTUs_tax, tgg, by="OTUID")
Dose_disease_state_Timepoint_Dose_Site_sigOTUs_alltaxr=merge(Dose_disease_state_Timepoint_Dose_Site_sigOTUs_alltax, allrsquare_complete, by.x="OTUID", by.y="row.names")
colnames(Dose_disease_state_Timepoint_Dose_Site_sigOTUs_alltaxr)[19]=c("pseudor2")
write.csv(Dose_disease_state_Timepoint_Dose_Site_sigOTUs_alltax, file="Output_files/Dose_disease_state_Timepoint_Dose_Site_sigOTUs_alltax.csv")

#Dose_disease_state_Site_Timepoint_Dose_Site
Dose_disease_state_Site_Timepoint_Dose_Site= as.matrix(allpvalues_complete[,17])
Dose_disease_state_Site_Timepoint_Dose_Site_padj=as.matrix(p.adjust(Dose_disease_state_Site_Timepoint_Dose_Site, method="fdr"))
colnames(Dose_disease_state_Site_Timepoint_Dose_Site_padj)=c("padj")
row.names(Dose_disease_state_Site_Timepoint_Dose_Site_padj)=row.names(allpvalues_complete)
Dose_disease_state_Site_Timepoint_Dose_Site_sigOTUs=subset(Dose_disease_state_Site_Timepoint_Dose_Site_padj, Dose_disease_state_Site_Timepoint_Dose_Site_padj<0.05)
Dose_disease_state_Site_Timepoint_Dose_Site_sigOTUs=as.data.frame(Dose_disease_state_Site_Timepoint_Dose_Site_sigOTUs)
Dose_disease_state_Site_Timepoint_Dose_Site_sigOTUs$OTUID=row.names(Dose_disease_state_Site_Timepoint_Dose_Site_sigOTUs)
Dose_disease_state_Site_Timepoint_Dose_Site_sigOTUs_tax=merge(Dose_disease_state_Site_Timepoint_Dose_Site_sigOTUs, ts, by="OTUID")
Dose_disease_state_Site_Timepoint_Dose_Site_sigOTUs_alltax=merge(Dose_disease_state_Site_Timepoint_Dose_Site_sigOTUs_tax, tgg, by="OTUID")
Dose_disease_state_Site_Timepoint_Dose_Site_sigOTUs_alltaxr=merge(Dose_disease_state_Site_Timepoint_Dose_Site_sigOTUs_alltax, allrsquare_complete, by.x="OTUID", by.y="row.names")
colnames(Dose_disease_state_Site_Timepoint_Dose_Site_sigOTUs_alltaxr)[19]=c("pseudor2")
write.csv(Dose_disease_state_Site_Timepoint_Dose_Site_sigOTUs_alltax, file="Output_files/Dose_disease_state_Site_Timepoint_Dose_Site_sigOTUs_alltax.csv")


write.csv(t(allpvalues_complete), file="Output_files/allpvalues_complete.csv")
write.csv(t(allrsquare_complete), file="Output_files/allrsquare_complete.csv")





