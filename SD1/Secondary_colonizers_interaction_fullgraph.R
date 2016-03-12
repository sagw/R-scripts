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

Secondary=read.csv(file="Output_files/GLMM_results/Secondary_colonizers_interaction.csv", header=T, row.names=1)


##Francisellaceae
P_F=Secondary[Secondary$family_s=="Francisellaceae",]
c2norm_P_F=c2norm[row.names(c2norm) %in% P_F$OTUID,]
b_P_F=t(c2norm_P_F)
b_P_F=as.data.frame(b_P_F)
s1_s=cbind.data.frame(s1$T3_disease_state, s1$Dose_disease_state, s1$Timepoint, s1$Dose_Site, s1$Site)
d_P_F=cbind.data.frame(s1_s, b_P_F)
colnames(d_P_F)[1]="T3_disease_state"
colnames(d_P_F)[2]="Dose_disease_state"
colnames(d_P_F)[3]="Timepoint"
colnames(d_P_F)[4]="Dose_Site"
colnames(d_P_F)[5]="Site"


#split by combos
d_P_F_split=split(d_P_F, with(d_P_F, interaction(Dose_disease_state,T3_disease_state,Timepoint, Dose_Site, Site)), drop = TRUE)
d_P_F_split=lapply(d_P_F_split, function(x) {x[1:5]=list(NULL);x})
list2env(d_P_F_split, environment())
full_mean=data.frame(NA)
full_mean$D_D_Dose_CK4=mean(as.matrix(Diseased.Diseased.dose.CK4.CK4))
full_mean$D_D_Dose_CK14=mean(as.matrix(Diseased.Diseased.dose.CK14.CK14))
full_mean$H_H_T1_CK4=mean(as.matrix(Healthy.Healthy.one.CK4.CK14))
full_mean$H_H_T1_CK14=mean(as.matrix(Healthy.Healthy.one.CK14.CK4))
full_mean$D_D_T3_CK4_CK4=mean(as.matrix(Diseased.Diseased.three.CK4.CK4))
full_mean$D_D_T3_CK14_CK4=mean(as.matrix(Diseased.Diseased.three.CK14.CK4))
full_mean$D_D_T3_CK14_CK14=mean(as.matrix(Diseased.Diseased.three.CK14.CK14))
full_mean$D_D_T3_CK4_CK14=mean(as.matrix(Diseased.Diseased.three.CK4.CK14))
full_mean=full_mean[-1]

se=function(x) sd(x)/sqrt(length(x))
full_se=data.frame(NA)
full_se$D_D_Dose_CK4=se(as.matrix(Diseased.Diseased.dose.CK4.CK4))
full_se$D_D_Dose_CK14=se(as.matrix(Diseased.Diseased.dose.CK14.CK14))
full_se$H_H_T1_CK4=se(as.matrix(Healthy.Healthy.one.CK4.CK14))
full_se$H_H_T1_CK14=se(as.matrix(Healthy.Healthy.one.CK14.CK4))
full_se$D_D_T3_CK4_CK4=se(as.matrix(Diseased.Diseased.three.CK4.CK4))
full_se$D_D_T3_CK14_CK4=se(as.matrix(Diseased.Diseased.three.CK14.CK4))
full_se$D_D_T3_CK14_CK14=se(as.matrix(Diseased.Diseased.three.CK14.CK14))
full_se$D_D_T3_CK4_CK14=se(as.matrix(Diseased.Diseased.three.CK4.CK14))
full_se=full_se[-1]

full_f=rbind(full_se, full_mean)
row.names(full_f)=c("SE", "mean")
full_f=t(full_f)
full_f=as.data.frame(full_f)
full_f$Dose_Site=c("CK4", "CK14", "NA", "NA", "CK4", 
                "CK14", "CK14", "CK4")
full_f$Site=c("NA", "NA", "CK4", "CK14", "CK4", 
                    "CK4", "CK14", "CK14")
full_f$Time=c("Inoculant", "Inoculant", "T1", "T1", "T3", "T3", "T3", "T3")
full_f$Family="Francisellaceae"
full_f$Number_OTUs_CK14_CK4="55"
full_f$Number_OTUs_CK4_CK14="2"

##Campylobacteraceae
P_c=Secondary[Secondary$family_s=="Campylobacteraceae",]
c2norm_P_c=c2norm[row.names(c2norm) %in% P_c$OTUID,]
b_P_c=t(c2norm_P_c)
b_P_c=as.data.frame(b_P_c)
s1_s=cbind.data.frame(s1$T3_disease_state, s1$Dose_disease_state, s1$Timepoint, s1$Dose_Site, s1$Site)
d_P_c=cbind.data.frame(s1_s, b_P_c)
colnames(d_P_c)[1]="T3_disease_state"
colnames(d_P_c)[2]="Dose_disease_state"
colnames(d_P_c)[3]="Timepoint"
colnames(d_P_c)[4]="Dose_Site"
colnames(d_P_c)[5]="Site"

#split by combos
d_P_c_split=split(d_P_c, with(d_P_c, interaction(Dose_disease_state,T3_disease_state,Timepoint, Dose_Site, Site)), drop = TRUE)
d_P_c_split=lapply(d_P_c_split, function(x) {x[1:5]=list(NULL);x})
list2env(d_P_c_split, environment())
#calculate means and se
full_mean=data.frame(NA)
full_mean$D_D_Dose_CK4=mean(as.matrix(Diseased.Diseased.dose.CK4.CK4))
full_mean$D_D_Dose_CK14=mean(as.matrix(Diseased.Diseased.dose.CK14.CK14))
full_mean$H_H_T1_CK4=mean(as.matrix(Healthy.Healthy.one.CK4.CK14))
full_mean$H_H_T1_CK14=mean(as.matrix(Healthy.Healthy.one.CK14.CK4))
full_mean$D_D_T3_CK4_CK4=mean(as.matrix(Diseased.Diseased.three.CK4.CK4))
full_mean$D_D_T3_CK14_CK4=mean(as.matrix(Diseased.Diseased.three.CK14.CK4))
full_mean$D_D_T3_CK14_CK14=mean(as.matrix(Diseased.Diseased.three.CK14.CK14))
full_mean$D_D_T3_CK4_CK14=mean(as.matrix(Diseased.Diseased.three.CK4.CK14))
full_mean=full_mean[-1]

se=function(x) sd(x)/sqrt(length(x))
full_se=data.frame(NA)
full_se$D_D_Dose_CK4=se(as.matrix(Diseased.Diseased.dose.CK4.CK4))
full_se$D_D_Dose_CK14=se(as.matrix(Diseased.Diseased.dose.CK14.CK14))
full_se$H_H_T1_CK4=se(as.matrix(Healthy.Healthy.one.CK4.CK14))
full_se$H_H_T1_CK14=se(as.matrix(Healthy.Healthy.one.CK14.CK4))
full_se$D_D_T3_CK4_CK4=se(as.matrix(Diseased.Diseased.three.CK4.CK4))
full_se$D_D_T3_CK14_CK4=se(as.matrix(Diseased.Diseased.three.CK14.CK4))
full_se$D_D_T3_CK14_CK14=se(as.matrix(Diseased.Diseased.three.CK14.CK14))
full_se$D_D_T3_CK4_CK14=se(as.matrix(Diseased.Diseased.three.CK4.CK14))
full_se=full_se[-1]

full_c=rbind(full_se, full_mean)
row.names(full_c)=c("SE", "mean")
full_c=t(full_c)
full_c=as.data.frame(full_c)
full_c$Dose_Site=c("CK4", "CK14", "NA", "NA", "CK4", 
                   "CK14", "CK14", "CK4")
full_c$Site=c("NA", "NA", "CK4", "CK14", "CK4", 
              "CK4", "CK14", "CK14")
full_c$Time=c("Inoculant", "Inoculant", "T1", "T1", "T3", "T3", "T3", "T3")
full_c$Family="Campylobacteraceae"
full_c$Number_OTUs_CK14_CK4="3"
full_c$Number_OTUs_CK4_CK14="358"

##Flavobacteriaceae
P_fl=Secondary[Secondary$family_s=="Flavobacteriaceae",]
c2norm_P_fl=c2norm[row.names(c2norm) %in% P_fl$OTUID,]
b_P_fl=t(c2norm_P_fl)
b_P_fl=as.data.frame(b_P_fl)
s1_s=cbind.data.frame(s1$T3_disease_state, s1$Dose_disease_state, s1$Timepoint, s1$Dose_Site, s1$Site)
d_P_fl=cbind.data.frame(s1_s, b_P_fl)
colnames(d_P_fl)[1]="T3_disease_state"
colnames(d_P_fl)[2]="Dose_disease_state"
colnames(d_P_fl)[3]="Timepoint"
colnames(d_P_fl)[4]="Dose_Site"
colnames(d_P_fl)[5]="Site"

#split by combos
d_P_fl_split=split(d_P_fl, with(d_P_fl, interaction(Dose_disease_state,T3_disease_state,Timepoint, Dose_Site, Site)), drop = TRUE)
d_P_fl_split=lapply(d_P_fl_split, function(x) {x[1:5]=list(NULL);x})
list2env(d_P_fl_split, environment())
#calculate means and se
full_mean=data.frame(NA)
full_mean$D_D_Dose_CK4=mean(as.matrix(Diseased.Diseased.dose.CK4.CK4))
full_mean$D_D_Dose_CK14=mean(as.matrix(Diseased.Diseased.dose.CK14.CK14))
full_mean$H_H_T1_CK4=mean(as.matrix(Healthy.Healthy.one.CK4.CK14))
full_mean$H_H_T1_CK14=mean(as.matrix(Healthy.Healthy.one.CK14.CK4))
full_mean$D_D_T3_CK4_CK4=mean(as.matrix(Diseased.Diseased.three.CK4.CK4))
full_mean$D_D_T3_CK14_CK4=mean(as.matrix(Diseased.Diseased.three.CK14.CK4))
full_mean$D_D_T3_CK14_CK14=mean(as.matrix(Diseased.Diseased.three.CK14.CK14))
full_mean$D_D_T3_CK4_CK14=mean(as.matrix(Diseased.Diseased.three.CK4.CK14))
full_mean=full_mean[-1]

se=function(x) sd(x)/sqrt(length(x))
full_se=data.frame(NA)
full_se$D_D_Dose_CK4=se(as.matrix(Diseased.Diseased.dose.CK4.CK4))
full_se$D_D_Dose_CK14=se(as.matrix(Diseased.Diseased.dose.CK14.CK14))
full_se$H_H_T1_CK4=se(as.matrix(Healthy.Healthy.one.CK4.CK14))
full_se$H_H_T1_CK14=se(as.matrix(Healthy.Healthy.one.CK14.CK4))
full_se$D_D_T3_CK4_CK4=se(as.matrix(Diseased.Diseased.three.CK4.CK4))
full_se$D_D_T3_CK14_CK4=se(as.matrix(Diseased.Diseased.three.CK14.CK4))
full_se$D_D_T3_CK14_CK14=se(as.matrix(Diseased.Diseased.three.CK14.CK14))
full_se$D_D_T3_CK4_CK14=se(as.matrix(Diseased.Diseased.three.CK4.CK14))
full_se=full_se[-1]

full_fl=rbind(full_se, full_mean)
row.names(full_fl)=c("SE", "mean")
full_fl=t(full_fl)
full_fl=as.data.frame(full_fl)
full_fl$Dose_Site=c("CK4", "CK14", "NA", "NA", "CK4", 
                   "CK14", "CK14", "CK4")
full_fl$Site=c("NA", "NA", "CK4", "CK14", "CK4", 
              "CK4", "CK14", "CK14")
full_fl$Time=c("Inoculant", "Inoculant", "T1", "T1", "T3", "T3", "T3", "T3")
full_fl$Family="Flavobacteriaceae"
full_fl$Number_OTUs_CK14_CK4="3"
full_fl$Number_OTUs_CK4_CK14="9"


##Alteromonadaceae
P_a=Secondary[Secondary$family_s=="Alteromonadaceae",]
c2norm_P_a=c2norm[row.names(c2norm) %in% P_a$OTUID,]
b_P_a=t(c2norm_P_a)
b_P_a=as.data.frame(b_P_a)
s1_s=cbind.data.frame(s1$T3_disease_state, s1$Dose_disease_state, s1$Timepoint, s1$Dose_Site, s1$Site)
d_P_a=cbind.data.frame(s1_s, b_P_a)
colnames(d_P_a)[1]="T3_disease_state"
colnames(d_P_a)[2]="Dose_disease_state"
colnames(d_P_a)[3]="Timepoint"
colnames(d_P_a)[4]="Dose_Site"
colnames(d_P_a)[5]="Site"

#split by combos
d_P_a_split=split(d_P_a, with(d_P_a, interaction(Dose_disease_state,T3_disease_state,Timepoint, Dose_Site, Site)), drop = TRUE)
d_P_a_split=lapply(d_P_a_split, function(x) {x[1:5]=list(NULL);x})
list2env(d_P_a_split, environment())
#calculate means and se
full_mean=data.frame(NA)
full_mean$D_D_Dose_CK4=mean(as.matrix(Diseased.Diseased.dose.CK4.CK4))
full_mean$D_D_Dose_CK14=mean(as.matrix(Diseased.Diseased.dose.CK14.CK14))
full_mean$H_H_T1_CK4=mean(as.matrix(Healthy.Healthy.one.CK4.CK14))
full_mean$H_H_T1_CK14=mean(as.matrix(Healthy.Healthy.one.CK14.CK4))
full_mean$D_D_T3_CK4_CK4=mean(as.matrix(Diseased.Diseased.three.CK4.CK4))
full_mean$D_D_T3_CK14_CK4=mean(as.matrix(Diseased.Diseased.three.CK14.CK4))
full_mean$D_D_T3_CK14_CK14=mean(as.matrix(Diseased.Diseased.three.CK14.CK14))
full_mean$D_D_T3_CK4_CK14=mean(as.matrix(Diseased.Diseased.three.CK4.CK14))
full_mean=full_mean[-1]

se=function(x) sd(x)/sqrt(length(x))
full_se=data.frame(NA)
full_se$D_D_Dose_CK4=se(as.matrix(Diseased.Diseased.dose.CK4.CK4))
full_se$D_D_Dose_CK14=se(as.matrix(Diseased.Diseased.dose.CK14.CK14))
full_se$H_H_T1_CK4=se(as.matrix(Healthy.Healthy.one.CK4.CK14))
full_se$H_H_T1_CK14=se(as.matrix(Healthy.Healthy.one.CK14.CK4))
full_se$D_D_T3_CK4_CK4=se(as.matrix(Diseased.Diseased.three.CK4.CK4))
full_se$D_D_T3_CK14_CK4=se(as.matrix(Diseased.Diseased.three.CK14.CK4))
full_se$D_D_T3_CK14_CK14=se(as.matrix(Diseased.Diseased.three.CK14.CK14))
full_se$D_D_T3_CK4_CK14=se(as.matrix(Diseased.Diseased.three.CK4.CK14))
full_se=full_se[-1]

full_a=rbind(full_se, full_mean)
row.names(full_a)=c("SE", "mean")
full_a=t(full_a)
full_a=as.data.frame(full_a)
full_a$Dose_Site=c("CK4", "CK14", "NA", "NA", "CK4", 
                   "CK14", "CK14", "CK4")
full_a$Site=c("NA", "NA", "CK4", "CK14", "CK4", 
              "CK4", "CK14", "CK14")
full_a$Time=c("Inoculant", "Inoculant", "T1", "T1", "T3", "T3", "T3", "T3")
full_a$Family="Alteromonadaceae"
full_a$Number_OTUs_CK14_CK4="7"
full_a$Number_OTUs_CK4_CK14="135"


##Rhodobacteraceae
P_r=Secondary[Secondary$family_s=="Rhodobacteraceae",]
c2norm_P_r=c2norm[row.names(c2norm) %in% P_r$OTUID,]
b_P_r=t(c2norm_P_r)
b_P_r=as.data.frame(b_P_r)
s1_s=cbind.data.frame(s1$T3_disease_state, s1$Dose_disease_state, s1$Timepoint, s1$Dose_Site, s1$Site)
d_P_r=cbind.data.frame(s1_s, b_P_r)
colnames(d_P_r)[1]="T3_disease_state"
colnames(d_P_r)[2]="Dose_disease_state"
colnames(d_P_r)[3]="Timepoint"
colnames(d_P_r)[4]="Dose_Site"
colnames(d_P_r)[5]="Site"

#split by combos
d_P_r_split=split(d_P_r, with(d_P_r, interaction(Dose_disease_state,T3_disease_state,Timepoint, Dose_Site, Site)), drop = TRUE)
d_P_r_split=lapply(d_P_r_split, function(x) {x[1:5]=list(NULL);x})
list2env(d_P_r_split, environment())
#calculate means and se
full_mean=data.frame(NA)
full_mean$D_D_Dose_CK4=mean(as.matrix(Diseased.Diseased.dose.CK4.CK4))
full_mean$D_D_Dose_CK14=mean(as.matrix(Diseased.Diseased.dose.CK14.CK14))
full_mean$H_H_T1_CK4=mean(as.matrix(Healthy.Healthy.one.CK4.CK14))
full_mean$H_H_T1_CK14=mean(as.matrix(Healthy.Healthy.one.CK14.CK4))
full_mean$D_D_T3_CK4_CK4=mean(as.matrix(Diseased.Diseased.three.CK4.CK4))
full_mean$D_D_T3_CK14_CK4=mean(as.matrix(Diseased.Diseased.three.CK14.CK4))
full_mean$D_D_T3_CK14_CK14=mean(as.matrix(Diseased.Diseased.three.CK14.CK14))
full_mean$D_D_T3_CK4_CK14=mean(as.matrix(Diseased.Diseased.three.CK4.CK14))
full_mean=full_mean[-1]

se=function(x) sd(x)/sqrt(length(x))
full_se=data.frame(NA)
full_se$D_D_Dose_CK4=se(as.matrix(Diseased.Diseased.dose.CK4.CK4))
full_se$D_D_Dose_CK14=se(as.matrix(Diseased.Diseased.dose.CK14.CK14))
full_se$H_H_T1_CK4=se(as.matrix(Healthy.Healthy.one.CK4.CK14))
full_se$H_H_T1_CK14=se(as.matrix(Healthy.Healthy.one.CK14.CK4))
full_se$D_D_T3_CK4_CK4=se(as.matrix(Diseased.Diseased.three.CK4.CK4))
full_se$D_D_T3_CK14_CK4=se(as.matrix(Diseased.Diseased.three.CK14.CK4))
full_se$D_D_T3_CK14_CK14=se(as.matrix(Diseased.Diseased.three.CK14.CK14))
full_se$D_D_T3_CK4_CK14=se(as.matrix(Diseased.Diseased.three.CK4.CK14))
full_se=full_se[-1]

full_r=rbind(full_se, full_mean)
row.names(full_r)=c("SE", "mean")
full_r=t(full_r)
full_r=as.data.frame(full_r)
full_r$Dose_Site=c("CK4", "CK14", "NA", "NA", "CK4", 
                   "CK14", "CK14", "CK4")
full_r$Site=c("NA", "NA", "CK4", "CK14", "CK4", 
              "CK4", "CK14", "CK14")
full_r$Time=c("Inoculant", "Inoculant", "T1", "T1", "T3", "T3", "T3", "T3")
full_r$Family="Rhodobacteraceae"
full_r$Number_OTUs_CK14_CK4="17"
full_r$Number_OTUs_CK4_CK14="68"


##Methylococcaceae
P_m=Secondary[Secondary$family_s=="Methylococcaceae",]
c2norm_P_m=c2norm[row.names(c2norm) %in% P_m$OTUID,]
b_P_m=t(c2norm_P_m)
b_P_m=as.data.frame(b_P_m)
s1_s=cbind.data.frame(s1$T3_disease_state, s1$Dose_disease_state, s1$Timepoint, s1$Dose_Site, s1$Site)
d_P_m=cbind.data.frame(s1_s, b_P_m)
colnames(d_P_m)[1]="T3_disease_state"
colnames(d_P_m)[2]="Dose_disease_state"
colnames(d_P_m)[3]="Timepoint"
colnames(d_P_m)[4]="Dose_Site"
colnames(d_P_m)[5]="Site"

#split by combos
d_P_m_split=split(d_P_m, with(d_P_m, interaction(Dose_disease_state,T3_disease_state,Timepoint, Dose_Site, Site)), drop = TRUE)
d_P_m_split=lapply(d_P_m_split, function(x) {x[1:5]=list(NULL);x})
list2env(d_P_m_split, environment())
#calculate means and se
full_mean=data.frame(NA)
full_mean$D_D_Dose_CK4=mean(as.matrix(Diseased.Diseased.dose.CK4.CK4))
full_mean$D_D_Dose_CK14=mean(as.matrix(Diseased.Diseased.dose.CK14.CK14))
full_mean$H_H_T1_CK4=mean(as.matrix(Healthy.Healthy.one.CK4.CK14))
full_mean$H_H_T1_CK14=mean(as.matrix(Healthy.Healthy.one.CK14.CK4))
full_mean$D_D_T3_CK4_CK4=mean(as.matrix(Diseased.Diseased.three.CK4.CK4))
full_mean$D_D_T3_CK14_CK4=mean(as.matrix(Diseased.Diseased.three.CK14.CK4))
full_mean$D_D_T3_CK14_CK14=mean(as.matrix(Diseased.Diseased.three.CK14.CK14))
full_mean$D_D_T3_CK4_CK14=mean(as.matrix(Diseased.Diseased.three.CK4.CK14))
full_mean=full_mean[-1]

se=function(x) sd(x)/sqrt(length(x))
full_se=data.frame(NA)
full_se$D_D_Dose_CK4=se(as.matrix(Diseased.Diseased.dose.CK4.CK4))
full_se$D_D_Dose_CK14=se(as.matrix(Diseased.Diseased.dose.CK14.CK14))
full_se$H_H_T1_CK4=se(as.matrix(Healthy.Healthy.one.CK4.CK14))
full_se$H_H_T1_CK14=se(as.matrix(Healthy.Healthy.one.CK14.CK4))
full_se$D_D_T3_CK4_CK4=se(as.matrix(Diseased.Diseased.three.CK4.CK4))
full_se$D_D_T3_CK14_CK4=se(as.matrix(Diseased.Diseased.three.CK14.CK4))
full_se$D_D_T3_CK14_CK14=se(as.matrix(Diseased.Diseased.three.CK14.CK14))
full_se$D_D_T3_CK4_CK14=se(as.matrix(Diseased.Diseased.three.CK4.CK14))
full_se=full_se[-1]

full_m=rbind(full_se, full_mean)
row.names(full_m)=c("SE", "mean")
full_m=t(full_m)
full_m=as.data.frame(full_m)
full_m$Dose_Site=c("CK4", "CK14", "NA", "NA", "CK4", 
                   "CK14", "CK14", "CK4")
full_m$Site=c("NA", "NA", "CK4", "CK14", "CK4", 
              "CK4", "CK14", "CK14")
full_m$Time=c("Inoculant", "Inoculant", "T1", "T1", "T3", "T3", "T3", "T3")
full_m$Family="Methylococcaceae"
full_m$Number_OTUs_CK14_CK4="11"
full_m$Number_OTUs_CK4_CK14="0"



##concatenate
full_full=rbind(full_a, full_m, full_f, full_fl, full_r, full_c)
full_full$Number_OTUs_CK14_CK4=as.numeric(full_full$Number_OTUs_CK14_CK4)
full_full$Number_OTUs_CK4_CK14=as.numeric(full_full$Number_OTUs_CK4_CK14)


full_full_dose=full_full[full_full$Time=="Inoculant",]
full_full_one=full_full[full_full$Time=="T1",]
full_full_T3=full_full[full_full$Time=="T3",]

##plot
dodge=position_dodge(width=0.9)
plot=ggplot() + 
  geom_bar(data=full_full_T3, aes(x=Family, y=mean, fill=Family), stat="identity", position="identity")+
  geom_errorbar(data=full_full_T3, aes(x=Family, y=mean, ymin=mean-SE, ymax=mean+SE), width=0.5, colour="black")+
  theme_bw()+ ylim(0,1)+ 
  scale_fill_manual(values=c("#E6A0C4", "#F1BB7B", "#5B1A18", "#3B9AB2", "#00A08A", "#E58601")) +
  facet_grid(Dose_Site~Site, scale="free_y") 
print(plot)
ggsave(plot=plot,height=6,width=6,dpi=1000, filename="Output_files/GLMM_results/plots/Secondary_interaction_T3.pdf", useDingbats=FALSE)

##plot doses
dodge=position_dodge(width=0.9)
plot=ggplot() + 
  geom_bar(data=full_full_dose, aes(x=Family, y=mean, fill=Family), stat="identity", position="identity")+
  geom_errorbar(data=full_full_dose, aes(x=Family, y=mean, ymin=mean-SE, ymax=mean+SE), width=0.5, colour="black")+
  theme_bw()+ ylim(0,1)+ 
  scale_fill_manual(values=c("#E6A0C4", "#F1BB7B", "#5B1A18", "#3B9AB2", "#00A08A", "#E58601")) +
  facet_grid(Dose_Site~., scale="free_y") 
print(plot)
ggsave(plot=plot,height=6,width=6,dpi=1000, filename="Output_files/GLMM_results/plots/Secondary_interaction_doses.pdf", useDingbats=FALSE)

##plotT1
dodge=position_dodge(width=0.9)
plot=ggplot() + 
  geom_bar(data=full_full_one, aes(x=Family, y=mean, fill=Family), stat="identity", position="identity")+
  geom_errorbar(data=full_full_one, aes(x=Family, y=mean, ymin=mean-SE, ymax=mean+SE), width=0.5, colour="black")+
  theme_bw()+ ylim(0,1)+ 
  scale_fill_manual(values=c("#E6A0C4", "#F1BB7B", "#5B1A18", "#3B9AB2", "#00A08A", "#E58601")) +
  facet_grid(Site~., scale="free_y") 
print(plot)
ggsave(plot=plot,height=6,width=6,dpi=1000, filename="Output_files/GLMM_results/plots/Secondary_interaction_T1.pdf", useDingbats=FALSE)
