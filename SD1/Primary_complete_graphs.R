rm(list=ls())
require("ggplot2")
require("reshape")
library(wesanderson)
require("plyr")
require("grid")

setwd("~/Documents/SD1/SD1_analyses/")

ts=read.csv("input_files/rep_set_tax_assignments_silva_clean2.csv")
colnames(ts)=c("OTUID", "kingdom_s", "phylum_s", "class_s", "order_s", "family_s", "genus_s", "species_s", "evalue", "ref")

tgg=read.csv("input_files/rep_set_tax_assignments_split2.csv")

s1=read.csv("input_files/SD1_map_clean.csv", header=TRUE, row.names=1)
s1=as.data.frame(s1)
c2norm=read.csv("input_files/c2norm.csv", header=TRUE, row.names=1)
#c2norm=as.matrix(sapply(c2norm, as.numeric)) 
check=cbind(rownames(s1), colnames(c2norm))

Primary=read.csv(file="Output_files/GLMM_results/Primary_pathogens_2.csv", header=T, row.names=1)

##Francisellaceae
P_F=Primary[Primary$family_s=="Francisellaceae",]
c2norm_P_F=c2norm[row.names(c2norm) %in% P_F$OTUID,]
b_P_F=t(c2norm_P_F)
b_P_F=as.data.frame(b_P_F)
s1_s=cbind.data.frame(s1$T3_disease_state, s1$Dose_disease_state, s1$Timepoint)
d_P_F=cbind.data.frame(s1_s, b_P_F)
colnames(d_P_F)[1]="T3_disease_state"
colnames(d_P_F)[2]="Dose_disease_state"
colnames(d_P_F)[3]="Timepoint"

#check order of samples
check=cbind(row.names(d_T3ex3), s1[1:1])

d_P_F_split=split(d_P_F, with(d_P_F, interaction(Dose_disease_state,T3_disease_state,Timepoint)), drop = TRUE)
d_P_F_split=lapply(d_P_F_split, function(x) {x[1:3]=list(NULL);x})
list2env(d_P_F_split, environment())
full_mean=data.frame(NA)
full_mean$D_D_Dose=mean(as.matrix(Diseased.Diseased.dose))
full_mean$D_D_T3=mean(as.matrix(Diseased.Diseased.three))
full_mean$D_D_T2=mean(as.matrix(Diseased.Diseased.two))
full_mean$D_H_T3=mean(as.matrix(Diseased.Healthy.three))
full_mean$D_H_T2=mean(as.matrix(Diseased.Healthy.two))
full_mean$H_H_dose=mean(as.matrix(Healthy.Healthy.dose))
full_mean$H_H_T1=mean(as.matrix(Healthy.Healthy.one))
full_mean$H_H_T3=mean(as.matrix(Healthy.Healthy.three))
full_mean$H_H_T2=mean(as.matrix(Healthy.Healthy.two))
full_mean=full_mean[-1]

se=function(x) sd(x)/sqrt(length(x))
full_se=data.frame(NA)
full_se$D_D_Dose=se(as.matrix(Diseased.Diseased.dose))
full_se$D_D_T3=se(as.matrix(Diseased.Diseased.three))
full_se$D_D_T2=se(as.matrix(Diseased.Diseased.two))
full_se$D_H_T3=se(as.matrix(Diseased.Healthy.three))
full_se$D_H_T2=se(as.matrix(Diseased.Healthy.two))
full_se$H_H_dose=se(as.matrix(Healthy.Healthy.dose))
full_se$H_H_T1=se(as.matrix(Healthy.Healthy.one))
full_se$H_H_T3=se(as.matrix(Healthy.Healthy.three))
full_se$H_H_T2=se(as.matrix(Healthy.Healthy.two))
full_se=full_se[-1]

full_f=rbind(full_se, full_mean)
row.names(full_f)=c("SE", "mean")
full_f=t(full_f)
full_f=as.data.frame(full_f)
full_f$Result=c("Dose", "Diseased", "Diseased", "Healthy", "Healthy", 
                "Control", "Healthy", "Healthy", "Healthy")
full_f$Treatments=c("Dosed", "Dosed", "Dosed", "Dosed", "Dosed", 
                    "Control", "Control", "Control", "Control")
full_f$Time=c("Inoculant", "T3", "T2", "T3", "T2", "Inoculant", "T1", "T3", "T2")
full_f$Family="Francisellaceae"
full_f$Number_OTUs="38"

##Pasteurellaceae####################################
P_p=Primary[Primary$family_s=="Pasteurellaceae",]
c2norm_P_p=c2norm[row.names(c2norm) %in% P_p$OTUID,]
b_P_p=t(c2norm_P_p)
b_P_p=as.data.frame(b_P_p)
s1_s=cbind.data.frame(s1$T3_disease_state, s1$Dose_disease_state, s1$Timepoint)
d_P_p=cbind.data.frame(s1_s, b_P_p)
colnames(d_P_p)[1]="T3_disease_state"
colnames(d_P_p)[2]="Dose_disease_state"
colnames(d_P_p)[3]="Timepoint"

d_P_p_split=split(d_P_p, with(d_P_p, interaction(Dose_disease_state,T3_disease_state,Timepoint)), drop = TRUE)
d_P_p_split=lapply(d_P_p_split, function(x) {x[1:3]=list(NULL);x})
list2env(d_P_p_split, environment())
full_mean=data.frame(NA)
full_mean$D_D_Dose=mean(as.matrix(Diseased.Diseased.dose))
full_mean$D_D_T3=mean(as.matrix(Diseased.Diseased.three))
full_mean$D_D_T2=mean(as.matrix(Diseased.Diseased.two))
full_mean$D_H_T3=mean(as.matrix(Diseased.Healthy.three))
full_mean$D_H_T2=mean(as.matrix(Diseased.Healthy.two))
full_mean$H_H_dose=mean(as.matrix(Healthy.Healthy.dose))
full_mean$H_H_T1=mean(as.matrix(Healthy.Healthy.one))
full_mean$H_H_T3=mean(as.matrix(Healthy.Healthy.three))
full_mean$H_H_T2=mean(as.matrix(Healthy.Healthy.two))
full_mean=full_mean[-1]

se=function(x) sd(x)/sqrt(length(x))
full_se=data.frame(NA)
full_se$D_D_Dose=se(as.matrix(Diseased.Diseased.dose))
full_se$D_D_T3=se(as.matrix(Diseased.Diseased.three))
full_se$D_D_T2=se(as.matrix(Diseased.Diseased.two))
full_se$D_H_T3=se(as.matrix(Diseased.Healthy.three))
full_se$D_H_T2=se(as.matrix(Diseased.Healthy.two))
full_se$H_H_dose=se(as.matrix(Healthy.Healthy.dose))
full_se$H_H_T1=se(as.matrix(Healthy.Healthy.one))
full_se$H_H_T3=se(as.matrix(Healthy.Healthy.three))
full_se$H_H_T2=se(as.matrix(Healthy.Healthy.two))
full_se=full_se[-1]

full_p=rbind(full_se, full_mean)
row.names(full_p)=c("SE", "mean")
full_p=t(full_p)
full_p=as.data.frame(full_p)
full_p$Result=c("Dose", "Diseased", "Diseased", "Healthy", "Healthy", 
                "Control", "Healthy", "Healthy", "Healthy")
full_p$Treatments=c("Dosed", "Dosed", "Dosed", "Dosed", "Dosed", 
                    "Control", "Control", "Control", "Control")
full_p$Time=c("Inoculant", "T3", "T2", "T3", "T2", "Inoculant", "T1", "T3", "T2")
full_p$Family="Pasteurellaceae"
full_p$Number_OTUs="27"


####Campylobacteraceae#######
P_c=Primary[Primary$family_s=="Campylobacteraceae",]
c2norm_P_c=c2norm[row.names(c2norm) %in% P_c$OTUID,]
b_P_c=t(c2norm_P_c)
b_P_c=as.data.frame(b_P_c)
s1_s=cbind.data.frame(s1$T3_disease_state, s1$Dose_disease_state, s1$Timepoint)
d_P_c=cbind.data.frame(s1_s, b_P_c)
colnames(d_P_c)[1]="T3_disease_state"
colnames(d_P_c)[2]="Dose_disease_state"
colnames(d_P_c)[3]="Timepoint"

d_P_c_split=split(d_P_c, with(d_P_c, interaction(Dose_disease_state,T3_disease_state,Timepoint)), drop = TRUE)
d_P_c_split=lapply(d_P_c_split, function(x) {x[1:3]=list(NULL);x})
list2env(d_P_c_split, environment())
full_mean=data.frame(NA)
full_mean$D_D_Dose=mean(as.matrix(Diseased.Diseased.dose))
full_mean$D_D_T3=mean(as.matrix(Diseased.Diseased.three))
full_mean$D_D_T2=mean(as.matrix(Diseased.Diseased.two))
full_mean$D_H_T3=mean(as.matrix(Diseased.Healthy.three))
full_mean$D_H_T2=mean(as.matrix(Diseased.Healthy.two))
full_mean$H_H_dose=mean(as.matrix(Healthy.Healthy.dose))
full_mean$H_H_T1=mean(as.matrix(Healthy.Healthy.one))
full_mean$H_H_T3=mean(as.matrix(Healthy.Healthy.three))
full_mean$H_H_T2=mean(as.matrix(Healthy.Healthy.two))
full_mean=full_mean[-1]

se=function(x) sd(x)/sqrt(length(x))
full_se=data.frame(NA)
full_se$D_D_Dose=se(as.matrix(Diseased.Diseased.dose))
full_se$D_D_T3=se(as.matrix(Diseased.Diseased.three))
full_se$D_D_T2=se(as.matrix(Diseased.Diseased.two))
full_se$D_H_T3=se(as.matrix(Diseased.Healthy.three))
full_se$D_H_T2=se(as.matrix(Diseased.Healthy.two))
full_se$H_H_dose=se(as.matrix(Healthy.Healthy.dose))
full_se$H_H_T1=se(as.matrix(Healthy.Healthy.one))
full_se$H_H_T3=se(as.matrix(Healthy.Healthy.three))
full_se$H_H_T2=se(as.matrix(Healthy.Healthy.two))
full_se=full_se[-1]

full_c=rbind(full_se, full_mean)
row.names(full_c)=c("SE", "mean")
full_c=t(full_c)
full_c=as.data.frame(full_c)
full_c$Result=c("Dose", "Diseased", "Diseased", "Healthy", "Healthy", 
                "Control", "Healthy", "Healthy", "Healthy")
full_c$Treatments=c("Dosed", "Dosed", "Dosed", "Dosed", "Dosed", 
                    "Control", "Control", "Control", "Control")
full_c$Time=c("Inoculant", "T3", "T2", "T3", "T2", "Inoculant", "T1", "T3", "T2")
full_c$Family="Campylobacteraceae"
full_c$Number_OTUs="25"

#####Alteromonadaceae
P_a=Primary[Primary$family_s=="Alteromonadaceae",]
c2norm_P_a=c2norm[row.names(c2norm) %in% P_a$OTUID,]
b_P_a=t(c2norm_P_a)
b_P_a=as.data.frame(b_P_a)
s1_s=cbind.data.frame(s1$T3_disease_state, s1$Dose_disease_state, s1$Timepoint)
d_P_a=cbind.data.frame(s1_s, b_P_a)
colnames(d_P_a)[1]="T3_disease_state"
colnames(d_P_a)[2]="Dose_disease_state"
colnames(d_P_a)[3]="Timepoint"

#check order of samples
check=cbind(row.names(d_T3ex3), s1[1:1])

d_P_a_split=split(d_P_a, with(d_P_a, interaction(Dose_disease_state,T3_disease_state,Timepoint)), drop = TRUE)
d_P_a_split=lapply(d_P_a_split, function(x) {x[1:3]=list(NULL);x})
list2env(d_P_a_split, environment())
full_mean=data.frame(NA)
full_mean$D_D_Dose=mean(as.matrix(Diseased.Diseased.dose))
full_mean$D_D_T3=mean(as.matrix(Diseased.Diseased.three))
full_mean$D_D_T2=mean(as.matrix(Diseased.Diseased.two))
full_mean$D_H_T3=mean(as.matrix(Diseased.Healthy.three))
full_mean$D_H_T2=mean(as.matrix(Diseased.Healthy.two))
full_mean$H_H_dose=mean(as.matrix(Healthy.Healthy.dose))
full_mean$H_H_T1=mean(as.matrix(Healthy.Healthy.one))
full_mean$H_H_T3=mean(as.matrix(Healthy.Healthy.three))
full_mean$H_H_T2=mean(as.matrix(Healthy.Healthy.two))
full_mean=full_mean[-1]

se=function(x) sd(x)/sqrt(length(x))
full_se=data.frame(NA)
full_se$D_D_Dose=se(as.matrix(Diseased.Diseased.dose))
full_se$D_D_T3=se(as.matrix(Diseased.Diseased.three))
full_se$D_D_T2=se(as.matrix(Diseased.Diseased.two))
full_se$D_H_T3=se(as.matrix(Diseased.Healthy.three))
full_se$D_H_T2=se(as.matrix(Diseased.Healthy.two))
full_se$H_H_dose=se(as.matrix(Healthy.Healthy.dose))
full_se$H_H_T1=se(as.matrix(Healthy.Healthy.one))
full_se$H_H_T3=se(as.matrix(Healthy.Healthy.three))
full_se$H_H_T2=se(as.matrix(Healthy.Healthy.two))
full_se=full_se[-1]

full_a=rbind(full_se, full_mean)
row.names(full_a)=c("SE", "mean")
full_a=t(full_a)
full_a=as.data.frame(full_a)
full_a$Result=c("Dose", "Diseased", "Diseased", "Healthy", "Healthy", 
                "Control", "Healthy", "Healthy", "Healthy")
full_a$Treatments=c("Dosed", "Dosed", "Dosed", "Dosed", "Dosed", 
                    "Control", "Control", "Control", "Control")
full_a$Time=c("Inoculant", "T3", "T2", "T3", "T2", "Inoculant", "T1", "T3", "T2")
full_a$Family="Alteromonadaceae"
full_a$Number_OTUs="22"


####Flavobacteriaceae
P_fl=Primary[Primary$family_s=="Flavobacteriaceae",]
c2norm_P_fl=c2norm[row.names(c2norm) %in% P_fl$OTUID,]
b_P_fl=t(c2norm_P_fl)
b_P_fl=as.data.frame(b_P_fl)
s1_s=cbind.data.frame(s1$T3_disease_state, s1$Dose_disease_state, s1$Timepoint)
d_P_fl=cbind.data.frame(s1_s, b_P_fl)
colnames(d_P_fl)[1]="T3_disease_state"
colnames(d_P_fl)[2]="Dose_disease_state"
colnames(d_P_fl)[3]="Timepoint"

d_P_fl_split=split(d_P_fl, with(d_P_fl, interaction(Dose_disease_state,T3_disease_state,Timepoint)), drop = TRUE)
d_P_fl_split=lapply(d_P_fl_split, function(x) {x[1:3]=list(NULL);x})
list2env(d_P_fl_split, environment())
full_mean=data.frame(NA)
full_mean$D_D_Dose=mean(as.matrix(Diseased.Diseased.dose))
full_mean$D_D_T3=mean(as.matrix(Diseased.Diseased.three))
full_mean$D_D_T2=mean(as.matrix(Diseased.Diseased.two))
full_mean$D_H_T3=mean(as.matrix(Diseased.Healthy.three))
full_mean$D_H_T2=mean(as.matrix(Diseased.Healthy.two))
full_mean$H_H_dose=mean(as.matrix(Healthy.Healthy.dose))
full_mean$H_H_T1=mean(as.matrix(Healthy.Healthy.one))
full_mean$H_H_T3=mean(as.matrix(Healthy.Healthy.three))
full_mean$H_H_T2=mean(as.matrix(Healthy.Healthy.two))
full_mean=full_mean[-1]

se=function(x) sd(x)/sqrt(length(x))
full_se=data.frame(NA)
full_se$D_D_Dose=se(as.matrix(Diseased.Diseased.dose))
full_se$D_D_T3=se(as.matrix(Diseased.Diseased.three))
full_se$D_D_T2=se(as.matrix(Diseased.Diseased.two))
full_se$D_H_T3=se(as.matrix(Diseased.Healthy.three))
full_se$D_H_T2=se(as.matrix(Diseased.Healthy.two))
full_se$H_H_dose=se(as.matrix(Healthy.Healthy.dose))
full_se$H_H_T1=se(as.matrix(Healthy.Healthy.one))
full_se$H_H_T3=se(as.matrix(Healthy.Healthy.three))
full_se$H_H_T2=se(as.matrix(Healthy.Healthy.two))
full_se=full_se[-1]

full_fl=rbind(full_se, full_mean)
row.names(full_fl)=c("SE", "mean")
full_fl=t(full_fl)
full_fl=as.data.frame(full_fl)
full_fl$Result=c("Dose", "Diseased", "Diseased", "Healthy", "Healthy", 
                "Control", "Healthy", "Healthy", "Healthy")
full_fl$Treatments=c("Dosed", "Dosed", "Dosed", "Dosed", "Dosed", 
                    "Control", "Control", "Control", "Control")
full_fl$Time=c("Inoculant", "T3", "T2", "T3", "T2", "Inoculant", "T1", "T3", "T2")
full_fl$Family="Flavobacteriaceae"
full_fl$Number_OTUs="22"


####Enterobacteriaceae
P_e=Primary[Primary$family_s=="Enterobacteriaceae",]
c2norm_P_e=c2norm[row.names(c2norm) %in% P_e$OTUID,]
b_P_e=t(c2norm_P_e)
b_P_e=as.data.frame(b_P_e)
s1_s=cbind.data.frame(s1$T3_disease_state, s1$Dose_disease_state, s1$Timepoint)
d_P_e=cbind.data.frame(s1_s, b_P_e)
colnames(d_P_e)[1]="T3_disease_state"
colnames(d_P_e)[2]="Dose_disease_state"
colnames(d_P_e)[3]="Timepoint"


d_P_e_split=split(d_P_e, with(d_P_e, interaction(Dose_disease_state,T3_disease_state,Timepoint)), drop = TRUE)
d_P_e_split=lapply(d_P_e_split, function(x) {x[1:3]=list(NULL);x})
list2env(d_P_e_split, environment())
full_mean=data.frame(NA)
full_mean$D_D_Dose=mean(as.matrix(Diseased.Diseased.dose))
full_mean$D_D_T3=mean(as.matrix(Diseased.Diseased.three))
full_mean$D_D_T2=mean(as.matrix(Diseased.Diseased.two))
full_mean$D_H_T3=mean(as.matrix(Diseased.Healthy.three))
full_mean$D_H_T2=mean(as.matrix(Diseased.Healthy.two))
full_mean$H_H_dose=mean(as.matrix(Healthy.Healthy.dose))
full_mean$H_H_T1=mean(as.matrix(Healthy.Healthy.one))
full_mean$H_H_T3=mean(as.matrix(Healthy.Healthy.three))
full_mean$H_H_T2=mean(as.matrix(Healthy.Healthy.two))
full_mean=full_mean[-1]

se=function(x) sd(x)/sqrt(length(x))
full_se=data.frame(NA)
full_se$D_D_Dose=se(as.matrix(Diseased.Diseased.dose))
full_se$D_D_T3=se(as.matrix(Diseased.Diseased.three))
full_se$D_D_T2=se(as.matrix(Diseased.Diseased.two))
full_se$D_H_T3=se(as.matrix(Diseased.Healthy.three))
full_se$D_H_T2=se(as.matrix(Diseased.Healthy.two))
full_se$H_H_dose=se(as.matrix(Healthy.Healthy.dose))
full_se$H_H_T1=se(as.matrix(Healthy.Healthy.one))
full_se$H_H_T3=se(as.matrix(Healthy.Healthy.three))
full_se$H_H_T2=se(as.matrix(Healthy.Healthy.two))
full_se=full_se[-1]

full_e=rbind(full_se, full_mean)
row.names(full_e)=c("SE", "mean")
full_e=t(full_e)
full_e=as.data.frame(full_e)
full_e$Result=c("Dose", "Diseased", "Diseased", "Healthy", "Healthy", 
                "Control", "Healthy", "Healthy", "Healthy")
full_e$Treatments=c("Dosed", "Dosed", "Dosed", "Dosed", "Dosed", 
                    "Control", "Control", "Control", "Control")
full_e$Time=c("Inoculant", "T3", "T2", "T3", "T2", "Inoculant", "T1", "T3", "T2")
full_e$Family="Enterobacteriaceae"
full_e$Number_OTUs="9"

####Rhodobacteraceae
P_r=Primary[Primary$family_s=="Rhodobacteraceae",]
c2norm_P_r=c2norm[row.names(c2norm) %in% P_r$OTUID,]
b_P_r=t(c2norm_P_r)
b_P_r=as.data.frame(b_P_r)
s1_s=cbind.data.frame(s1$T3_disease_state, s1$Dose_disease_state, s1$Timepoint)
d_P_r=cbind.data.frame(s1_s, b_P_r)
colnames(d_P_r)[1]="T3_disease_state"
colnames(d_P_r)[2]="Dose_disease_state"
colnames(d_P_r)[3]="Timepoint"


d_P_r_split=split(d_P_r, with(d_P_r, interaction(Dose_disease_state,T3_disease_state,Timepoint)), drop = TRUE)
d_P_r_split=lapply(d_P_r_split, function(x) {x[1:3]=list(NULL);x})
list2env(d_P_r_split, environment())
full_mean=data.frame(NA)
full_mean$D_D_Dose=mean(as.matrix(Diseased.Diseased.dose))
full_mean$D_D_T3=mean(as.matrix(Diseased.Diseased.three))
full_mean$D_D_T2=mean(as.matrix(Diseased.Diseased.two))
full_mean$D_H_T3=mean(as.matrix(Diseased.Healthy.three))
full_mean$D_H_T2=mean(as.matrix(Diseased.Healthy.two))
full_mean$H_H_dose=mean(as.matrix(Healthy.Healthy.dose))
full_mean$H_H_T1=mean(as.matrix(Healthy.Healthy.one))
full_mean$H_H_T3=mean(as.matrix(Healthy.Healthy.three))
full_mean$H_H_T2=mean(as.matrix(Healthy.Healthy.two))
full_mean=full_mean[-1]

se=function(x) sd(x)/sqrt(length(x))
full_se=data.frame(NA)
full_se$D_D_Dose=se(as.matrix(Diseased.Diseased.dose))
full_se$D_D_T3=se(as.matrix(Diseased.Diseased.three))
full_se$D_D_T2=se(as.matrix(Diseased.Diseased.two))
full_se$D_H_T3=se(as.matrix(Diseased.Healthy.three))
full_se$D_H_T2=se(as.matrix(Diseased.Healthy.two))
full_se$H_H_dose=se(as.matrix(Healthy.Healthy.dose))
full_se$H_H_T1=se(as.matrix(Healthy.Healthy.one))
full_se$H_H_T3=se(as.matrix(Healthy.Healthy.three))
full_se$H_H_T2=se(as.matrix(Healthy.Healthy.two))
full_se=full_se[-1]

full_r=rbind(full_se, full_mean)
row.names(full_r)=c("SE", "mean")
full_r=t(full_r)
full_r=as.data.frame(full_r)
full_r$Result=c("Dose", "Diseased", "Diseased", "Healthy", "Healthy", 
                "Control", "Healthy", "Healthy", "Healthy")
full_r$Treatments=c("Dosed", "Dosed", "Dosed", "Dosed", "Dosed", 
                    "Control", "Control", "Control", "Control")
full_r$Time=c("Inoculant", "T3", "T2", "T3", "T2", "Inoculant", "T1", "T3", "T2")
full_r$Family="Rhodobacteraceae"
full_r$Number_OTUs="7"

##Rickettsiaceae
####Rickettsiaceae
P_ri=Primary[Primary$family_s=="mitochondria",]
c2norm_P_ri=c2norm[row.names(c2norm) %in% P_ri$OTUID,]
b_P_ri=t(c2norm_P_ri)
b_P_ri=as.data.frame(b_P_ri)
s1_s=cbind.data.frame(s1$T3_disease_state, s1$Dose_disease_state, s1$Timepoint)
d_P_ri=cbind.data.frame(s1_s, b_P_ri)
colnames(d_P_ri)[1]="T3_disease_state"
colnames(d_P_ri)[2]="Dose_disease_state"
colnames(d_P_ri)[3]="Timepoint"


d_P_ri_split=split(d_P_ri, with(d_P_ri, interaction(Dose_disease_state,T3_disease_state,Timepoint)), drop = TRUE)
d_P_ri_split=lapply(d_P_ri_split, function(x) {x[1:3]=list(NULL);x})
list2env(d_P_ri_split, environment())
full_mean=data.frame(NA)
full_mean$D_D_Dose=mean(as.matrix(Diseased.Diseased.dose))
full_mean$D_D_T3=mean(as.matrix(Diseased.Diseased.three))
full_mean$D_D_T2=mean(as.matrix(Diseased.Diseased.two))
full_mean$D_H_T3=mean(as.matrix(Diseased.Healthy.three))
full_mean$D_H_T2=mean(as.matrix(Diseased.Healthy.two))
full_mean$H_H_dose=mean(as.matrix(Healthy.Healthy.dose))
full_mean$H_H_T1=mean(as.matrix(Healthy.Healthy.one))
full_mean$H_H_T3=mean(as.matrix(Healthy.Healthy.three))
full_mean$H_H_T2=mean(as.matrix(Healthy.Healthy.two))
full_mean=full_mean[-1]

se=function(x) sd(x)/sqrt(length(x))
full_se=data.frame(NA)
full_se$D_D_Dose=se(as.matrix(Diseased.Diseased.dose))
full_se$D_D_T3=se(as.matrix(Diseased.Diseased.three))
full_se$D_D_T2=se(as.matrix(Diseased.Diseased.two))
full_se$D_H_T3=se(as.matrix(Diseased.Healthy.three))
full_se$D_H_T2=se(as.matrix(Diseased.Healthy.two))
full_se$H_H_dose=se(as.matrix(Healthy.Healthy.dose))
full_se$H_H_T1=se(as.matrix(Healthy.Healthy.one))
full_se$H_H_T3=se(as.matrix(Healthy.Healthy.three))
full_se$H_H_T2=se(as.matrix(Healthy.Healthy.two))
full_se=full_se[-1]

full_ri=rbind(full_se, full_mean)
row.names(full_ri)=c("SE", "mean")
full_ri=t(full_ri)
full_ri=as.data.frame(full_ri)
full_ri$Result=c("Dose", "Diseased", "Diseased", "Healthy", "Healthy", 
                "Control", "Healthy", "Healthy", "Healthy")
full_ri$Treatments=c("Dosed", "Dosed", "Dosed", "Dosed", "Dosed", 
                    "Control", "Control", "Control", "Control")
full_ri$Time=c("Inoculant", "T3", "T2", "T3", "T2", "Inoculant", "T1", "T3", "T2")

full_ri$Family="Rickettsiaceae"
full_ri$Number_OTUs="5"

###Vibrio
####Vibrionaceae
P_v=Primary[Primary$family_s=="Vibrionaceae",]
c2norm_P_v=c2norm[row.names(c2norm) %in% P_v$OTUID,]
b_P_v=t(c2norm_P_v)
b_P_v=as.data.frame(b_P_v)
s1_s=cbind.data.frame(s1$T3_disease_state, s1$Dose_disease_state, s1$Timepoint)
d_P_v=cbind.data.frame(s1_s, b_P_v)
colnames(d_P_v)[1]="T3_disease_state"
colnames(d_P_v)[2]="Dose_disease_state"
colnames(d_P_v)[3]="Timepoint"


d_P_v_split=split(d_P_v, with(d_P_v, interaction(Dose_disease_state,T3_disease_state,Timepoint)), drop = TRUE)
d_P_v_split=lapply(d_P_v_split, function(x) {x[1:3]=list(NULL);x})
list2env(d_P_v_split, environment())
full_mean=data.frame(NA)
full_mean$D_D_Dose=mean(as.matrix(Diseased.Diseased.dose))
full_mean$D_D_T3=mean(as.matrix(Diseased.Diseased.three))
full_mean$D_D_T2=mean(as.matrix(Diseased.Diseased.two))
full_mean$D_H_T3=mean(as.matrix(Diseased.Healthy.three))
full_mean$D_H_T2=mean(as.matrix(Diseased.Healthy.two))
full_mean$H_H_dose=mean(as.matrix(Healthy.Healthy.dose))
full_mean$H_H_T1=mean(as.matrix(Healthy.Healthy.one))
full_mean$H_H_T3=mean(as.matrix(Healthy.Healthy.three))
full_mean$H_H_T2=mean(as.matrix(Healthy.Healthy.two))
full_mean=full_mean[-1]

se=function(x) sd(x)/sqrt(length(x))
full_se=data.frame(NA)
full_se$D_D_Dose=se(as.matrix(Diseased.Diseased.dose))
full_se$D_D_T3=se(as.matrix(Diseased.Diseased.three))
full_se$D_D_T2=se(as.matrix(Diseased.Diseased.two))
full_se$D_H_T3=se(as.matrix(Diseased.Healthy.three))
full_se$D_H_T2=se(as.matrix(Diseased.Healthy.two))
full_se$H_H_dose=se(as.matrix(Healthy.Healthy.dose))
full_se$H_H_T1=se(as.matrix(Healthy.Healthy.one))
full_se$H_H_T3=se(as.matrix(Healthy.Healthy.three))
full_se$H_H_T2=se(as.matrix(Healthy.Healthy.two))
full_se=full_se[-1]

full_v=rbind(full_se, full_mean)
row.names(full_v)=c("SE", "mean")
full_v=t(full_v)
full_v=as.data.frame(full_v)
full_v$Result=c("Dose", "Diseased", "Diseased", "Healthy", "Healthy", 
                    "Control", "Healthy", "Healthy", "Healthy")
full_v$Treatments=c("Dosed", "Dosed", "Dosed", "Dosed", "Dosed", 
                                        "Control", "Control", "Control", "Control")
full_v$Time=c("Inoculant", "T3", "T2", "T3", "T2", "Inoculant", "T1", "T3", "T2")
full_v$Family="Vibrionaceae"
full_v$Number_OTUs="1"

##concatenate
full_full=rbind(full_a, full_c, full_e, full_f, full_fl, full_p, full_r, full_ri, full_v)
full_full$Number_OTUs=as.numeric(full_full$Number_OTUs)

full_full_dose=full_full[full_full$Time=="Inoculant",]
full_full_one=full_full[full_full$Time=="T1",]
full_full_one$Treatments="Dosed"
full_full_T123=full_full[!full_full$Time=="Inoculant",]
full_full_1234=rbind(full_full_T123, full_full_one)
##plot
dodge=position_dodge(width=0.9)
pdf(file="Output_files/GLMM_results/plots/Primary_colonizers_complete.pdf", width=6, height=20)
plot=ggplot(data=full_full_1234, aes(x=Time, y=mean, group=Result, shape=Result, colour=Family)) + 
  geom_point(aes(size=Number_OTUs))+
  geom_errorbar(aes(ymin=mean-SE, ymax=mean+SE), width=0.05, colour="black")+
  theme_bw()+ facet_grid(Family+Treatments~., scale="free_y")+ ylim(-1,22)+
  geom_line()+ scale_shape_manual(values=c( 19, 1)) + 
  scale_color_manual(values=c("#E6A0C4", "#F1BB7B", "#0B775E", "#5B1A18", "#3B9AB2", "#798E87", 
                              "#E58601", "#7294D4", "#FD6467")) +
  theme(plot.margin=unit(c(5,5,5,5), "cm")) +
  scale_x_discrete(expand=c(0.5, 0))
plot
ggsave(plot=plot,height=30,width=8,dpi=1000, filename="Output_files/GLMM_results/plots/Primary_colonizers_Time.pdf", useDingbats=FALSE)

##plot doses
dodge=position_dodge(width=0.9)
plot=ggplot(data=full_full_dose, aes(x=Time, y=mean, shape=Result, colour=Family)) + 
  geom_point(aes(size=Number_OTUs))+
  geom_errorbar(aes(ymin=mean-SE, ymax=mean+SE), width=0.05, colour="black")+
  theme_bw()+ facet_grid(Family+Treatments~., scale="free_y")+ ylim(-1,15)+
  scale_shape_manual(values=c( 1, 19)) + 
  scale_color_manual(values=c("#E6A0C4", "#F1BB7B", "#0B775E", "#5B1A18", "#3B9AB2", "#798E87", 
                              "#E58601", "#7294D4", "#FD6467")) +
  scale_x_discrete(expand=c(0.5, 0))
plot
ggsave(plot=plot,height=15,width=4,dpi=1000, filename="Output_files/GLMM_results/plots/Primary_colonizers_Doses.pdf", useDingbats=FALSE)
