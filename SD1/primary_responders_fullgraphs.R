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

Primary=read.csv(file="Output_files/GLMM_results/Primary_responders2.csv")




####Cryomorphaceae#######
P_c=Primary[Primary$family_s=="Cryomorphaceae",]
c2norm_P_c=c2norm[row.names(c2norm) %in% P_c$OTUID,]
b_P_c=t(c2norm_P_c)
b_P_c=as.data.frame(b_P_c)
s1_s=cbind.data.frame(s1$T3_disease_state, s1$Dose_disease_state, s1$Timepoint)
d_P_c=cbind.data.frame(s1_s, b_P_c)
colnames(d_P_c)[1]="T3_disease_state"
colnames(d_P_c)[2]="Dose_disease_state"
colnames(d_P_c)[3]="Timepoint"

d_P_c_split=split(d_P_c, with(d_P_c, interaction(Dose_disease_state,T3_disease_state,Timepoint)), drop = TRUE)
d_P_c_split=lapply(d_P_c_split, function(x) {x[1:4]=list(NULL);x})
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
full_c$Family="Cryomorphaceae"
full_c$Number_OTUs="22"

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
full_a$Number_OTUs="24"


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
full_fl$Number_OTUs="26"


####Saprospiraceae
P_s=Primary[Primary$family_s=="Saprospiraceae",]
c2norm_P_s=c2norm[row.names(c2norm) %in% P_s$OTUID,]
b_P_s=t(c2norm_P_s)
b_P_s=as.data.frame(b_P_s)
s1_s=cbind.data.frame(s1$T3_disease_state, s1$Dose_disease_state, s1$Timepoint)
d_P_s=cbind.data.frame(s1_s, b_P_s)
colnames(d_P_s)[1]="T3_disease_state"
colnames(d_P_s)[2]="Dose_disease_state"
colnames(d_P_s)[3]="Timepoint"


d_P_s_split=split(d_P_s, with(d_P_s, interaction(Dose_disease_state,T3_disease_state,Timepoint)), drop = TRUE)
d_P_s_split=lapply(d_P_s_split, function(x) {x[1:3]=list(NULL);x})
list2env(d_P_s_split, environment())
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

full_s=rbind(full_se, full_mean)
row.names(full_s)=c("SE", "mean")
full_s=t(full_s)
full_s=as.data.frame(full_s)
full_s$Result=c("Dose", "Diseased", "Diseased", "Healthy", "Healthy", 
                "Control", "Healthy", "Healthy", "Healthy")
full_s$Treatments=c("Dosed", "Dosed", "Dosed", "Dosed", "Dosed", 
                    "Control", "Control", "Control", "Control")
full_s$Time=c("Inoculant", "T3", "T2", "T3", "T2", "Inoculant", "T1", "T3", "T2")
full_s$Family="Saprospiraceae"
full_s$Number_OTUs="20"

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
full_r$Number_OTUs="18"

##concatenate
full_full=rbind(full_a, full_c, full_s, full_fl, full_r)
full_full$Number_OTUs=as.numeric(full_full$Number_OTUs)

full_full_dose=full_full[full_full$Time=="Inoculant",]
full_full_one=full_full[full_full$Time=="T1",]
full_full_one$Treatments="Dosed"
full_full_T123=full_full[!full_full$Time=="Inoculant",]
full_full_1234=rbind(full_full_T123, full_full_one)
##plot
dodge=position_dodge(width=0.9)
plot=ggplot(data=full_full_1234, aes(x=Time, y=mean, group=Result, shape=Result, colour=Family)) + 
  geom_point(aes(size=Number_OTUs))+
  geom_errorbar(aes(ymin=mean-SE, ymax=mean+SE), width=0.05, colour="black")+
  theme_bw()+ facet_grid(Family+Treatments~., scale="free_y")+ ylim(-1,22)+
  geom_line()+ scale_shape_manual(values=c( 19, 1)) + 
  scale_color_manual(values=c("#E6A0C4", "#CCC591","#5B1A18","#E58601", "#D8A499")) +
  theme(plot.margin=unit(c(5,5,5,5), "cm")) + scale_size_continuous(range=c(4,5)) +
  scale_x_discrete(expand=c(0.5, 0))
plot
ggsave(plot=plot,height=21,width=8,dpi=1000, filename="Output_files/GLMM_results/plots/Primary_responders_Time.pdf", useDingbats=FALSE)

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
