rm(list=ls())
require(vegan)
require(DESeq2)
require(genefilter)
require(car)
require(MASS)

setwd("~/Documents/SD1/SD1_analyses/")

c=read.csv("input_files//filtered_otutable.csv", h=T, row.names=1)
c=c[5:NCOL(c)]

#information about samples
s=read.csv("input_files//SD1_map.csv", h=T, as.is=T)

b=t(c)
d=cbind.data.frame(s,b)

#check order of samples
check=cbind(row.names(d), d[1:1])

#get rid of dose disease stayed healthy
d2d=d[d$Timepoint== "two"& d$Dose_disease_state=="Diseased" & d$T3_disease_state=="Diseased",]
d2h=d[d$Timepoint== "two"& d$Dose_disease_state=="Healthy" & d$T3_disease_state=="Healthy",]
d2=rbind(d2d, d2h)

c2=d2[14:NCOL(d2)]
c2=t(c2)
rs = rowSums (c2)
use = (rs > 0)
c2 = c2 [ use, ]

s2=d2[,1:13]

#create dds
dds =DESeqDataSetFromMatrix(countData= c2, colData=s2, design= ~T3_disease_state)
colData(dds)$T3_disease_state= factor(colData(dds)$T3_disease_state, levels=c("Healthy","Diseased"))

#calculate arithemtic means because of zeros.
arith_mean=function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0])))
arithmeans = apply(c2, 1, arith_mean) 

#estimate size factors
dds = estimateSizeFactors(dds, geoMeans=arithmeans, locfunc=shorth)

#estimate dispersions
dds=estimateDispersions(dds)
plotDispEsts(dds)

#calculate p values using lrt through deseq
dds=nbinomLRT(dds, reduced=~1)

#get normalized counts from DESEQ
c2norm=counts(dds, normalized=TRUE)

#set filter parameters
filter=rowMeans(c2norm)

#extract deseq results
res=results(dds, alpha=0.05, cooksCutoff=FALSE, pAdjustMethod="fdr", filter=filter)

padjdeseq=res$padj
padjdeseq=as.data.frame(padjdeseq)
rownames(padjdeseq)=rownames(res)

#non-adjusted p values
pvalsdeseq=as.data.frame(res$pvalue)

##################deseq code for filtering############
lowerQuantile <- mean(filter == 0)
if (lowerQuantile < 0.95) upperQuantile <- 0.95 else upperQuantile <- 1
theta <- seq(0, .95, length=20)

filtPadj <- filtered_p(filter=filter, test=res$pvalue,
                       theta=theta, method="fdr") 
numRej <- colSums(filtPadj < .05, na.rm = TRUE)

lo.fit <- lowess(numRej ~ theta, f=1/5)
# two ways of calculating j, our deseq uses second one
#one
#if (max(numRej) <= 10) {
#  j <- 1
#} else { 
#  residual <- if (all(numRej==0)) {
#    0
#  } else {
#    numRej[numRej > 0] - lo.fit$y[numRej > 0]
#  }
#  thresh <- max(lo.fit$y) - sqrt(mean(residual^2))
#  j <- if (any(numRej > thresh)) {
#    which(numRej > thresh)[1]
#  } else {
#    1  
#  }
#}
#two (use this for now)
j <- which.max(numRej) 

#get padj values and other parameters for comparison
padj <- filtPadj[, j, drop=TRUE]
cutoffs <- quantile(filter, theta)
filterThreshold <- cutoffs[j]
filterNumRej <- data.frame(theta=theta, numRej=numRej)
filterTheta <- theta[j]

#make sure same threshold caluclated using straight up deseq vs code above
attr(res,"filterThreshold")
filterThreshold

####Quasipoisson####################################################################
#make empty dataframe for pvalues
allpvalues_qp =matrix(NA,nrow=1,ncol=nrow(c2norm))
colnames(allpvalues_qp)=row.names(c2norm)

#run quasiposson glm and anova on each otu
for (i in 1:ncol(allpvalues_qp)) {
  vec = as.numeric(c2norm[i,])
  combined=cbind.data.frame(s2, vec)
  glmtest=glm(vec~T3_disease_state, data=combined, family="quasipoisson")
  anovaglm=as.data.frame(Anova(glmtest, type=c("II"), test.statistic="LR"))
  colnames(anovaglm)=c("chisq", "df", "pvalue")
  anovaglm$effect=row.names(anovaglm)
  pvalues=anovaglm$pvalue
  pvalues=as.data.frame(pvalues)
  pvalues=apply(pvalues, 2, as.numeric)
  pvalues=as.data.frame(pvalues)
  row.names(pvalues)=row.names(anovaglm)
  pvalues=t(pvalues)
  T3_disease_state_p=pvalues[,"T3_disease_state"]
  allpvalues_qp[,i] = T3_disease_state_p
}

#extract p values and combine with deseq
row.names(allpvalues_qp)="qp_p"
qp_pvalues=t(allpvalues_qp)

deseq_qp_p=merge(pvalsdeseq, qp_pvalues, by.x="row.names", by.y="row.names")
colnames(deseq_qp_p)=c("OTU", "deseq_p", "qp_p")

##filter according to deseq and adjust using fdr
##################deseq code for determining theta############
lowerQuantile = mean(filter == 0)
if (lowerQuantile < 0.95) upperQuantile = 0.95 else upperQuantile = 1
theta = seq(0, .95, length=20)

filtPadj = filtered_p(filter=filter, test=allpvalues_qp,
                       theta=theta, method="fdr") 
numRej = colSums(filtPadj < .05, na.rm = TRUE)

lo.fit = lowess(numRej ~ theta, f=1/5)
# two ways of calculating j, our deseq uses second one
#one
#if (max(numRej) <= 10) {
#  j <- 1
#} else { 
#  residual <- if (all(numRej==0)) {
#    0
#  } else {
#    numRej[numRej > 0] - lo.fit$y[numRej > 0]
#  }
#  thresh <- max(lo.fit$y) - sqrt(mean(residual^2))
#  j <- if (any(numRej > thresh)) {
#    which(numRej > thresh)[1]
#  } else {
#    1  
#  }
#}
#two (use this for now)
j <- which.max(numRej) 

#get filtered p values
qp_padjfilt = filtPadj[, j, drop=TRUE]
filterThreshold = cutoffs[j]

#qp_padjfilt=filtered_p(filter=filter, test=allpvalues, theta=.85, method="fdr")
qp_padjfilt=as.data.frame(qp_padjfilt)
colnames(qp_padjfilt)="qp"
row.names(qp_padjfilt)=colnames(allpvalues_qp)

#make dataframe with deseq res
deseq_qp=merge(padjdeseq, qp_padjfilt, by.x="row.names", by.y="row.names")


############run nb glm################################################
#make empty dataframe for pvalues
allpvalues =matrix(NA,nrow=1,ncol=nrow(c2norm))
colnames(allpvalues)=row.names(c2norm)

#run negativebinomial glm and anova on each otu
for (i in 1:ncol(allpvalues)) {
  vec = as.numeric(c2norm[i,])
  combined=cbind.data.frame(s2, vec)
  glmtest=glm.nb(vec~T3_disease_state, data=combined)
  anovaglm=as.data.frame(Anova(glmtest, type=c("II"), test.statistic="LR"))
  colnames(anovaglm)=c("chisq", "df", "pvalue")
  anovaglm$effect=row.names(anovaglm)
  pvalues=anovaglm$pvalue
  pvalues=as.data.frame(pvalues)
  pvalues=apply(pvalues, 2, as.numeric)
  pvalues=as.data.frame(pvalues)
  row.names(pvalues)=row.names(anovaglm)
  pvalues=t(pvalues)
  T3_disease_state_p=pvalues[,"T3_disease_state"]
  allpvalues[,i] = T3_disease_state_p
}

#extract p values  
row.names(allpvalues)="nb"
nb_pvalues=t(allpvalues)
deseq_nb_p=merge(pvalsdeseq, nb_pvalues, by.x="row.names", by.y="row.names")
colnames(deseq_nb_p)=c("OTUID", "deseq_p", "nb_p")

##filter according to deseq and adjust using fdr
##################deseq code for determining theta############
lowerQuantile = mean(filter == 0)
if (lowerQuantile < 0.95) upperQuantile = 0.95 else upperQuantile = 1
theta = seq(0, .95, length=20)

filtPadj = filtered_p(filter=filter, test=allpvalues,
                      theta=theta, method="fdr") 
numRej = colSums(filtPadj < .05, na.rm = TRUE)

lo.fit = lowess(numRej ~ theta, f=1/5)
# two ways of calculating j our deseq uses the second one
##adjust using fdr 
#if (max(numRej) <= 10) {
#  j = 1
#} else { 
#  residual = if (all(numRej==0)) {
#    0
#  } else {
#    numRej[numRej > 0] - lo.fit$y[numRej > 0]
#  }
#  thresh = max(lo.fit$y) - sqrt(mean(residual^2))
#  j = if (any(numRej > thresh)) {
#    which(numRej > thresh)[1]
#  } else {
#    1  
#  }
}
#use this j for now
j = which.max(numRej) 

nb_padjfilt = filtPadj[, j, drop=TRUE]
filterThreshold = cutoffs[j]
nb_padjfilt=as.data.frame(nb_padjfilt)
colnames(nb_padjfilt)="nb"
row.names(nb_padjfilt)=colnames(allpvalues)

deseq_nb=merge(padjdeseq, nb_padjfilt, by.x="row.names", by.y="row.names")

##run comparisons of all three methods
##deseq vs qp p val
plot(deseq_qp_p$qp_p, deseq_qp_p$deseq_p)
deseq_qp_lm=lm(deseq_p~qp_p, data=deseq_qp_p)
summary(deseq_qp_lm)

##deseq vs nb p val
plot(deseq_nb_p$nb_p, deseq_nb_p$deseq_p)
deseq_nb_lm=lm(deseq_p~nb_p, data=deseq_nb_p)
summary(deseq_nb_lm)

##deseq vs qppadj
plot(deseq_qp$qp, deseq_qp$padjdeseq)
deseq_qp_lm=lm(padjdeseq~qp, data=deseq_qp)
summary(deseq_qp_lm)
deseq_qp$qp=as.numeric(deseq_qp$qp)

##deseq vs nbpadj
plot(deseq_nb$nb, deseq_nb$padjdeseq)
deseq_nb_lm=lm(padjdeseq~nb, data=deseq_nb)
summary(deseq_nb_lm)
deseq_nb$nb=as.numeric(deseq_nb$nb)

##Identify problem otus
problems=deseq_nb[is.na(deseq_nb$padjdeseq)==TRUE & is.na(deseq_nb$nb)==FALSE & deseq_nb$nb<0.05,]
sig_qp=deseq_qp[deseq_qp$qp<0.05 & !is.na(deseq_qp$qp),]
sig_deseq=deseq_qp[deseq_qp$padjdeseq<0.05 & !is.na(deseq_qp$padjdeseq),]
sig_nb=deseq_nb[deseq_nb$nb<0.05 & !is.na(deseq_nb$nb),]
problems=deseq_qp[is.na(deseq_qp$padjdeseq) & !is.na(deseq_qp$qp) & deseq_qp$qp<.05,]
sum(is.na(deseq_qp$qp))
sum(is.na(deseq_qp$padjdeseq))