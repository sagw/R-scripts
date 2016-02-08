rm(list=ls())
require(MASS)
require(lme4)

setwd("~/Documents/SD1/SD1_analyses/")

s1=read.csv("Input_files/sample_data.csv", header=TRUE, row.names=1)
c2norm=read.csv("Input_files/c2norm.csv", header=TRUE, row.names=1)
c2norm=as.matrix(sapply(c2norm, as.numeric)) 

#remove otus found fewer than 10 times
c2normfiltered=c2norm[apply(c2norm, MARGIN=1, function(x) (sum(x>0))>=10),]

#subset 100 otus (or however many)
c2normfilteredsub=head(c2normfiltered, n=100)

##test with three otus
#full model
for (i in c(1,3,4)) { #:nrow(c2normfilteredsub)) {
  tryCatch({
    vec = as.numeric(c2normfilteredsub[i,])
    combined=cbind.data.frame(s1, vec)
    glmtest=glmer(vec~T3_disease_state*Dose_disease_state+ Site*Timepoint*Dose_Site*Dose_disease_state+(1|Dose/Timepoint), 
                     data=combined, glmerControl(optCtrl = list(maxfun = 1000000)), family=quasipoisson)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


#remove dose disease state
for (i in c(1,3,4)) { #:nrow(c2normfilteredsub)) {
  tryCatch({
    vec = as.numeric(c2normfilteredsub[i,])
    combined=cbind.data.frame(s1, vec)
    glmtest=glmer.nb(vec~T3_disease_state*Site*Timepoint*Dose_Site+(1|Timepoint)+(1|Dose/Genotype), 
                     data=combined, glmerControl(optCtrl = list(maxfun = 1000000)))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

#remove nesting
for (i in c(1,3,4)) { #:nrow(c2normfilteredsub)) {
  tryCatch({
    vec = as.numeric(c2normfilteredsub[i,])
    combined=cbind.data.frame(s1, vec)
    glmtest=glmer.nb(vec~T3_disease_state*Site*Timepoint*Dose_Site+(1|Timepoint)+(1|Dose)+(1|Genotype), 
                     data=combined, glmerControl(optCtrl = list(maxfun = 1000000)))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

#remove Dose
for (i in c(1,3,4)) { #:nrow(c2normfilteredsub)) {
  tryCatch({
    vec = as.numeric(c2normfilteredsub[i,])
    combined=cbind.data.frame(s1, vec)
    glmtest=glmer.nb(vec~T3_disease_state*Site*Timepoint*Dose_Site+(1|Timepoint)+(1|Genotype), 
                     data=combined, glmerControl(optCtrl = list(maxfun = 1000000)))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

#remove genotype
for (i in c(1,3,4)) { #:nrow(c2normfilteredsub)) {
  tryCatch({
    vec = as.numeric(c2normfilteredsub[i,])
    combined=cbind.data.frame(s1, vec)
    glmtest=glmer.nb(vec~T3_disease_state*Site*Timepoint*Dose_Site+(1|Timepoint)+(1|Dose)+(1|Genotype), 
                     data=combined, glmerControl(optCtrl = list(maxfun = 1000000)))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

#remove timepoint
for (i in c(1,3,4)) { #:nrow(c2normfilteredsub)) {
  tryCatch({
    vec = as.numeric(c2normfilteredsub[i,])
    combined=cbind.data.frame(s1, vec)
    glmtest=glmer.nb(vec~T3_disease_state*Site*Timepoint*Dose_Site+(1|Dose/Genotype), 
                     data=combined, glmerControl(optCtrl = list(maxfun = 1000000)))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}



#remove both genotype and dose
for (i in 1:nrow(c2normfilteredsub)) {
  tryCatch({
    vec = as.numeric(c2normfilteredsub[i,])
    combined=cbind.data.frame(s1, vec)
    glmtest=glmer.nb(vec~T3_disease_state*Site*Timepoint*Dose_Site+(1|Timepoint), 
                     data=combined, glmerControl(optCtrl = list(maxfun = 1000000)))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}



#remove both genotype and dose try poisson dist?
for (i in 1:nrow(c2normfilteredsub)) {
  tryCatch({
    vec = as.numeric(c2normfilteredsub[i,])
    combined=cbind.data.frame(s1, vec)
    glmtest=glmer(vec~T3_disease_state*Site*Timepoint*Dose_Site+(1|Timepoint), 
                  data=combined, family="poisson", glmerControl(optCtrl = list(maxfun = 1000000)))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

#full model, poisson distribution
for (i in c(1,3,4)) { #:nrow(c2normfilteredsub)) {
  tryCatch({
    vec = as.numeric(c2normfilteredsub[i,])
    combined=cbind.data.frame(s1, vec)
    glmtest=glmer(vec~T3_disease_state*Site*Timepoint*Dose_Site+(1|Timepoint)+(1|Dose/Genotype), 
                  data=combined, family="poisson", glmerControl(optCtrl = list(maxfun = 1000000)))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}





for (i in c(1,3,4)) { #:nrow(c2normfilteredsub)) {
  tryCatch({
    vec = as.numeric(c2normfilteredsub[i,])
    combined=cbind.data.frame(s1, vec)
    glmtest=glmer.nb(vec~T3_disease_state*Site*Timepoint*Dose_Site+(1|Dose)+(1|Genotype)+(1|Timepoint), 
                     data=combined, glmerControl(optCtrl = list(maxfun = 1000000)))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

#full model raise the tol threshold
for (i in c(1,3,4)) { #:nrow(c2normfilteredsub)) {
  tryCatch({
    vec = as.numeric(c2normfilteredsub[i,])
    combined=cbind.data.frame(s1, vec)
    glmtest=glmer.nb(vec~T3_disease_state*Site*Dose_disease_state*Timepoint*Dose_Site+(1|Timepoint)+(1|Dose/Genotype), 
                     data=combined, glmerControl(tolPwrss=1e-1))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}



