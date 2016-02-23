rm(list=ls())
require("ggplot2")
require("reshape")
library(wesanderson)

setwd("~/Documents/SD1/SD1_analyses/")

ts=read.csv("input_files/rep_set_tax_assignments_silva_clean2.csv")
colnames(ts)=c("OTUID", "kingdom_s", "phylum_s", "class_s", "order_s", "family_s", "genus_s", "species_s", "evalue", "ref")

tgg=read.csv("input_files/rep_set_tax_assignments_split2.csv")

s1=read.csv("input_files/SD1_map_clean.csv", header=TRUE, row.names=1)
s1=as.data.frame(s1)
c2norm=read.csv("input_files/c2norm.csv", header=TRUE, row.names=1)
#c2norm=as.matrix(sapply(c2norm, as.numeric)) 
check=cbind(rownames(s1), colnames(c2norm))

##Import all sig otu files 
Genotype_sigOTUs_alltaxr=read.csv("Output_files/GLMM_results/Genotype_sigOTUs_alltax.csv")
T3_disease_state_sigOTUs_alltaxr=read.csv("Output_files/GLMM_results/T3_disease_state_sigOTUs_alltax.csv", header=TRUE, row.names=1)
Dose_disease_state_sigOTUs_alltaxr=read.csv("Output_files/GLMM_results/Dose_disease_state_sigOTUs_alltax.csv", header=TRUE, row.names=1)
Site_sigOTUs_alltaxr=read.csv("Output_files/GLMM_results/Site_sigOTUs_alltax.csv", header=TRUE, row.names=1)
Timepoint_sigOTUs_alltaxr=read.csv("Output_files/GLMM_results/Timepoint_sigOTUs_alltax.csv", header=TRUE, row.names=1)
Dose_Site_sigOTUs_alltaxr=read.csv("Output_files/GLMM_results/Dose_Site_sigOTUs_alltax.csv", header=TRUE, row.names=1)
Site_Dose_disease_state_sigOTUs_alltaxr=read.csv("Output_files/GLMM_results/Site_Dose_disease_state_sigOTUs_alltax.csv", header=TRUE, row.names=1)
Site_Timepoint_sigOTUs_alltaxr=read.csv("Output_files/GLMM_results/Site_Timepoint_sigOTUs_alltax.csv", header=TRUE, row.names=1)
Dose_disease_state_Timepoint_sigOTUs_alltaxr=read.csv("Output_files/GLMM_results/Dose_disease_state_Timepoint_sigOTUs_alltax.csv", header=TRUE, row.names=1)
Site_Dose_Site_sigOTUs_alltaxr=read.csv("Output_files/GLMM_results/Site_Dose_Site_sigOTUs_alltax.csv", header=TRUE, row.names=1)
Dose_disease_state_Dose_Site_sigOTUs_alltaxr=read.csv("Output_files/GLMM_results/Dose_disease_state_Dose_Site_sigOTUs_alltax.csv", header=TRUE, row.names=1)
Timepoint_Dose_Site_sigOTUs_alltaxr=read.csv("Output_files/GLMM_results/Timepoint_Dose_Site_sigOTUs_alltax.csv", header=TRUE, row.names=1)
Site_Dose_disease_state_Timepoint_sigOTUs_alltaxr=read.csv("Output_files/GLMM_results/Site_Dose_disease_state_Timepoint_sigOTUs_alltax.csv", header=TRUE, row.names=1)
Site_Dose_disease_state_Dose_site_sigOTUs_alltaxr=read.csv("Output_files/GLMM_results/Site_Dose_disease_state_Dose_Site_sigOTUs_alltax.csv", header=TRUE, row.names=1)
Site_Timepoint_Dose_Site_sigOTUs_alltaxr=read.csv("Output_files/GLMM_results/Site_Timepoint_Dose_Site_sigOTUs_alltax.csv", header=TRUE, row.names=1)
Dose_disease_state_Timepoint_Dose_Site_sigOTUs_alltaxr=read.csv("Output_files/GLMM_results/Dose_disease_state_Timepoint_Dose_Site_sigOTUs_alltax.csv", header=TRUE, row.names=1)
Site_Dose_disease_state_Timepoint_Dose_Site_sigOTUs_alltaxr=read.csv("Output_files/GLMM_results/Site_Dose_disease_state_Timepoint_Dose_Site_sigOTUs_alltax.csv", header=TRUE, row.names=1)


#beneficial
#OTUs significant in genotype and dose
Geno_T3=merge(Genotype_sigOTUs_alltaxr, T3_disease_state_sigOTUs_alltaxr, by="OTUID")
Geno_dose=merge(Genotype_sigOTUs_alltaxr, Dose_disease_state_sigOTUs_alltaxr, by="OTUID")
Geno_T3_time=merge(Geno_T3, Timepoint_sigOTUs_alltaxr, by="OTUID")

##Find OTUs up in healthy
c2norm_genoT3=c2norm[row.names(c2norm) %in% Geno_T3$OTUID, ]
b_genoT3=t(c2norm_genoT3)
d_genoT3=cbind.data.frame(s1$Dose_disease_state, s1$T3_disease_state, s1$Timepoint, s1$Genotype, b_genoT3)
colnames(d_genoT3)[1]="Dose_disease_state"
colnames(d_genoT3)[2]="T3_disease_state"
colnames(d_genoT3)[3]="Timepoint"
colnames(d_genoT3)[4]="Genotype"
d_genoT3=d_genoT3[d_genoT3$Timepoint %in% c("one","two", "three"),]

#average diseased and healthy
d_genoT3_dh=aggregate(. ~ Dose_disease_state:T3_disease_state:Timepoint:Genotype, data=d_genoT3, FUN=mean)
row.names(d_genoT3_dh)=paste(d_genoT3_dh$Dose_disease_state, d_genoT3_dh$T3_disease_state, d_genoT3_dh$Timepoint, sep="_")
d_genoT3_dh=d_genoT3_dh[4:NCOL(d_genoT3_dh)]
d_genoT3_dh=t(d_genoT3_dh)
d_genoT3_dh=as.data.frame(d_genoT3_dh)
geno_beneficial=subset(d_genoT3_dh, d_genoT3_dh$Healthy_Healthy_one>0&d_genoT3_dh$Healthy_Healthy_two>0&d_genoT3_dh$Healthy_Healthy_three>0)
geno_beneficial$DDT3vsHHT1=geno_beneficial$Healthy_Healthy_one-geno_beneficial$Diseased_Diseased_three
geno_beneficial$DDT3vsHHT3=geno_beneficial$Healthy_Healthy_three-geno_beneficial$Diseased_Diseased_three
geno_beneficial$DDT2vsHHT2=geno_beneficial$Healthy_Healthy_two-geno_beneficial$Diseased_Diseased_two
geno_beneficial$DHT2vsHHT2=geno_beneficial$Healthy_Healthy_two-geno_beneficial$Diseased_Healthy_two
geno_beneficial$DHT3vsHHT3=geno_beneficial$Healthy_Healthy_three-geno_beneficial$Diseased_Healthy_three
geno_beneficial_tax=merge(geno_beneficial, ts, by.x="row.names", by.y="OTUID")
geno_beneficial2_tax=subset(geno_beneficial_tax, geno_beneficial_tax$DDT3vsHHT3>0&geno_beneficial_tax$DDT2vsHHT2>0&geno_beneficial_tax$DDT3vsHHT1>0)

geno_beneficial2_defensive_tax=subset(geno_beneficial2_tax, geno_beneficial2_tax$DHT2vsHHT2>0&geno_beneficial2_tax$DHT3vsHHT3)

c2norm_genoBen=c2norm[row.names(c2norm) %in% geno_beneficial2_tax$Row.names, ]
b_genoBen=t(c2norm_genoBen)
d_genoBen=cbind.data.frame(s1$Dose_disease_state, s1$T3_disease_state, s1$Timepoint, s1$Genotype, b_genoBen)
colnames(d_genoBen)[1]="Dose_disease_state"
colnames(d_genoBen)[2]="T3_disease_state"
colnames(d_genoBen)[3]="Timepoint"
colnames(d_genoBen)[4]="Genotype"
d_genoBen=d_genoBen[d_genoBen$Timepoint %in% c("one","two", "three"),]
d_genoBen_agg=aggregate(. ~ T3_disease_state:Genotype, data=d_genoBen, FUN=mean)
row.names(d_genoBen_agg)=paste(d_genoBen_agg$T3_disease_state, d_genoBen_agg$Genotype, sep="_")
d_genoBen_agg=d_genoBen_agg[5:NCOL(d_genoBen_agg)]
d_genoBen_agg=t(d_genoBen_agg)
d_genoBen_agg=as.data.frame(d_genoBen_agg)
genoBen_tax=merge(d_genoBen_agg, ts, by.x="row.names", by.y="OTUID")
genoBen_tax$OTUID=row.names(genoBen_tax)
table(genoBen_tax$family_s)
genoBen_tax_endo=genoBen_tax[genoBen_tax$family_s=="Hahellaceae",]
colnames(genoBen_tax_endo)=c("OTUID", "D_B_CK14", "H_B_CK14", "D_B_CK4", "H_B_CK4", "D_G_CK14",
                             "H_G_CK14", "D_G_CK4", "H_G_CK4", "D_O_CK14","H_O_CK14",
                             "D_O_CK4", "H_O_CK4", "D_P_CK14", "H_P_CK14", "D_P_CK4",
                             "H_P_CK4", "D_W_CK14", "H_W_CK14",  "D_W_CK4", "H_W_CK4",
                             "kingdom_s", "phylum_s", "class_s", "order_s", "family_s", "genus_s",
                             "species_s", "evalue", "ref", "no")

otu.long_h = melt(genoBen_tax_endo, id = c("OTUID", "family_s"), measure = c("D_B_CK14", "H_B_CK14", "D_B_CK4", "H_B_CK4", "D_G_CK14",
                                                                        "H_G_CK14", "D_G_CK4", "H_G_CK4", "D_O_CK14","H_O_CK14",
                                                                        "D_O_CK4", "H_O_CK4", "D_P_CK14", "H_P_CK14", "D_P_CK4",
                                                                        "H_P_CK4", "D_W_CK14", "H_W_CK14",  "D_W_CK4", "H_W_CK4" ))
otu.long_h$OTUID = factor(otu.long_h$OTUID, levels = otu.long_h$OTUID[order(otu.long_h$family_s)])
plotgb=ggplot(data=otu.long_h, aes(x=OTUID, y=value )) +
  geom_bar(aes(fill=family_s), stat="identity") + 
  theme_bw()+ theme(legend.position = "left") +
  ggtitle("Genotype specific healthy symbionts") +
  facet_grid(.~variable, scale="free_y") + ylim(0,45)+coord_flip() 

  plotgb



d_genoT3dose_dh$DDT3vsHHT3=d_genoT3dose_dh$Healthy_Healthy_three-d_genoT3dose_dh$Diseased_Diseased_three
d_genoT3dose_dh$DDT2vsHHT2=d_genoT3dose_dh$Healthy_Healthy_two-d_genoT3dose_dh$Diseased_Diseased_two
d_genoT3dose_dh$DDT2vsDHT2=d_genoT3dose_dh$Diseased_Healthy_two-d_genoT3dose_dh$Diseased_Diseased_two
d_genoT3dose_dh$HHT2vsDHT2=d_genoT3dose_dh$Diseased_Healthy_two-d_genoT3dose_dh$Healthy_Healthy_two
d_genoT3dose_dh_tax=merge(d_genoT3dose_dh, ts, by.x="row.names", by.y="OTUID")
geno_beneficial=subset(d_genoT3dose_dh_tax, d_genoT3dose_dh_tax$DDT3vsHHT3>0|d_genoT3dose_dh_tax$DDT2vsHHT2>0|d_genoT3dose_dh_tax$DDT2vsDHT2>0)
row.names(geno_beneficial)=geno_beneficial$Row.names
geno_beneficial=geno_beneficial[,12:23]
geno_beneficial$type=ifelse(geno_beneficial$DDT3vsHHT3>=0&geno_beneficial$DDT2vsHHT2>=0&geno_beneficial$DDT3vsHHT1>0, "Healthy_symbiont", "Defensive" )
#geno_beneficial$Row.names %in% Dose_disease_state_Timepoint_sigOTUs_alltaxr$OTUID
geno_beneficial$OTUID=row.names(geno_beneficial)
Healthy_symbionts=geno_beneficial[geno_beneficial$type=="Healthy_symbiont",]
Defensive=geno_beneficial[geno_beneficial$type=="Defensive",]
Defensive=Defensive[!Defensive$OTUID=="New.ReferenceOTU68613",]

#plot geno_beneficial
#pdf(file="Output_files/Taxonomy_plots/sFamily_CK4.pdf", width=10, height=15)

otu.long_h = melt(Healthy_symbionts, id = c("OTUID", "family_s"), measure = c("DDT3vsHHT1", "DDT3vsHHT3", "DDT2vsHHT2", "DDT2vsDHT2", "HHT2vsDHT2"))
otu.long_h$OTUID = factor(otu.long_h$OTUID, levels = otu.long_h$OTUID[order(otu.long_h$family_s)])
plotgb=ggplot(data=otu.long_h, aes(x=OTUID, y=value )) +
  geom_bar(aes(fill=family_s), stat="identity") + coord_flip() +
  theme_bw()+theme(legend.position = "left") +
  ggtitle("Genotype specific healthy symbionts") +
  facet_grid(.~variable, scale="free_y") + ylim(-1,1) +geom_hline(xintercept=0)
plotgb

otu.long = melt(Defensive, id = c("OTUID", "family_s", "type"), measure = c("DDT3vsHHT1", "DDT3vsHHT3", "DDT2vsHHT2", "DDT2vsDHT2", "HHT2vsDHT2"))
otu.long$OTUID <- factor(otu.long$OTUID, levels = otu.long$OTUID[order(otu.long$family_s)])
plotgb=ggplot(data=otu.long, aes(x=OTUID, y=value )) +
  geom_bar(aes(fill=family_s), stat="identity") + coord_flip() +
  theme_bw()+theme(legend.position = "left") +
  ggtitle("Genotype specific defensive bacteria") +
  facet_grid(.~variable, scale="free_y") + ylim(-1,1) +geom_hline(xintercept=0)
plotgb

##not genotype specific healthy symbionts
#OTUs significant in genotype and dose
T3dose=T3_disease_state_Dose_disease_state_sigOTUs_alltaxr[!T3_disease_state_Dose_disease_state_sigOTUs_alltaxr$OTUID %in% Genotype_sigOTUs_alltaxr$OTUID,]

##Find OTUs up in healthy
c2norm_T3dose=c2norm[row.names(c2norm) %in% T3dose$OTUID, ]
b_T3dose=t(c2norm_T3dose)
d_T3dose=cbind.data.frame(s1$Dose_disease_state, s1$T3_disease_state, s1$Timepoint, b_T3dose)
colnames(d_T3dose)[1]="Dose_disease_state"
colnames(d_T3dose)[2]="T3_disease_state"
colnames(d_T3dose)[3]="Timepoint"
#d_genoT3dose=d_genoT3dose[d_genoT3dose$Timepoint %in% c("two", "three"),]

#average diseased and healthy
d_T3dose_dh=aggregate(. ~ Dose_disease_state:T3_disease_state:Timepoint, data=d_T3dose, FUN=mean)
row.names(d_T3dose_dh)=paste(d_T3dose_dh$Dose_disease_state, d_T3dose_dh$T3_disease_state, d_T3dose_dh$Timepoint, sep="_")
d_T3dose_dh=d_T3dose_dh[4:NCOL(d_T3dose_dh)]
d_T3dose_dh=t(d_T3dose_dh)
d_T3dose_dh=as.data.frame(d_T3dose_dh)
d_T3dose_dh$DDT3vsHHT1=d_T3dose_dh$Healthy_Healthy_one-d_T3dose_dh$Diseased_Diseased_three
d_T3dose_dh$DDT3vsHHT3=d_T3dose_dh$Healthy_Healthy_three-d_T3dose_dh$Diseased_Diseased_three
d_T3dose_dh$DDT2vsHHT2=d_T3dose_dh$Healthy_Healthy_two-d_T3dose_dh$Diseased_Diseased_two
d_T3dose_dh$DDT2vsDHT2=d_T3dose_dh$Diseased_Healthy_two-d_T3dose_dh$Diseased_Diseased_two
d_T3dose_dh$HHT2vsDHT2=d_T3dose_dh$Diseased_Healthy_two-d_T3dose_dh$Healthy_Healthy_two
d_T3dose_dh_tax=merge(d_T3dose_dh, ts, by.x="row.names", by.y="OTUID")
beneficial=subset(d_T3dose_dh_tax, d_T3dose_dh_tax$DDT3vsHHT3>0|d_T3dose_dh_tax$DDT2vsHHT2>0|d_T3dose_dh_tax$DDT2vsDHT2>0)
row.names(beneficial)=beneficial$Row.names
#beneficial$Row.names %in% Dose_disease_state_Timepoint_sigOTUs_alltaxr$OTUID
beneficial=beneficial[,12:23]
beneficial$OTUID=row.names(beneficial)
beneficial$type=ifelse(beneficial$DDT2vsDHT2>=0&beneficial$DDT2vsHHT2>=0, "Defensive", "Healthy_symbiont" )
Healthy_symbionts=beneficial[beneficial$type=="Healthy_symbiont",]
Defensive=beneficial[beneficial$type=="Defensive",]

#plot geno_beneficial
#pdf(file="Output_files/Taxonomy_plots/sFamily_CK4.pdf", width=10, height=15)

otu.long_h = melt(Healthy_symbionts, id = c("OTUID", "family_s"), measure = c("DDT3vsHHT1", "DDT3vsHHT3", "DDT2vsHHT2", "DDT2vsDHT2", "HHT2vsDHT2"))
plotgb=ggplot(data=otu.long_h, aes(x=OTUID, y=value )) +
  geom_bar(aes(fill=family_s), stat="identity") + coord_flip() +
  theme_bw()+theme(legend.position = "left") +
  ggtitle("Not genotype specific healthy symbionts") +
  facet_grid(.~variable, scale="free_y") + ylim(-1,1) +geom_hline(xintercept=0)
plotgb

otu.long = melt(Defensive, id = c("OTUID", "family_s", "type"), measure = c("DDT3vsHHT1", "DDT3vsHHT3", "DDT2vsHHT2", "DDT2vsDHT2", "HHT2vsDHT2"))
otu.long$OTUID <- factor(otu.long$OTUID, levels = otu.long$OTUID[order(otu.long$family_s)])
plotgb=ggplot(data=otu.long, aes(x=OTUID, y=value )) +
  geom_bar(aes(fill=family_s), stat="identity") + coord_flip() +
  theme_bw()+theme(legend.position = "left") +
  ggtitle("Not genotype specific defensive bacteria") +
  facet_grid(.~variable, scale="free_y") + ylim(-1,1) +geom_hline(xintercept=0)
plotgb


##########Disease############################

##primary pathogen
int1=T3_disease_state_Dose_disease_state_sigOTUs_alltaxr[!T3_disease_state_Dose_disease_state_sigOTUs_alltaxr$OTUID %in% Dose_disease_state_sigOTUs_alltaxr,]
int2=T3_disease_state_sigOTUs_alltaxr[!T3_disease_state_sigOTUs_alltaxr$OTUID %in% Dose_disease_state_sigOTUs_alltaxr,]

T3dnoTime=T3_disease_state_Dose_disease_state_sigOTUs_alltaxr[!T3_disease_state_Dose_disease_state_sigOTUs_alltaxr$OTUID %in% Dose_disease_state_Timepoint_sigOTUs_alltaxr$OTUID,]

DoseTime=Dose_disease_state_sigOTUs_alltaxr[Dose_disease_state_sigOTUs_alltaxr$OTUID %in% Dose_disease_state_Timepoint_sigOTUs_alltaxr$OTUID,]

##Create sig otu otu table
c2norm_T3=c2norm[row.names(c2norm) %in% T3dnoTime$OTUID, ]
b_T3=t(c2norm_T3)
d_T3=cbind.data.frame(s1$T3_disease_state, s1$Dose_disease_state, b_T3)
colnames(d_T3)[1]="T3_disease_state"
colnames(d_T3)[2]="Dose_disease_state"

#check order of samples
check=cbind(row.names(d_T3), s1[1:1])

#average diseased and healthy
##T3 disease state
d_T3_dh=aggregate(. ~T3_disease_state:Dose_disease_state, data=d_T3, FUN=mean)
row.names(d_T3_dh)=paste(d_T3_dh$Dose_disease_state, d_T3_dh$T3_disease_state, sep="_")
d_T3_dh=d_T3_dh[,-c(1,2)]
d_T3_dh=t(d_T3_dh)
d_T3_dh=as.data.frame(d_T3_dh)
d_T3_dh$DvH=d_T3_dh$Diseased_Healthy-d_T3_dh$Healthy_Healthy
UpD=subset(d_T3_dh, d_T3_dh$DvH>0)
UpDtax=merge(UpD, ts, by.x="row.names", by.y="OTUID")


##Dose disease state
d_dose_dh=aggregate(. ~Dose_disease_state , data=d_T3, FUN=mean)
row.names(d_dose_dh)=d_dose_dh[,1]
d_dose_dh=d_dose_dh[,-c(1,2)]
d_dose_dh=t(d_dose_dh)
d_dose_dh=as.data.frame(d_dose_dh)
d_dose_dh$DvH=d_dose_dh$Diseased-d_dose_dh$Healthy
UpDdose=subset(d_dose_dh, d_dose_dh$DvH>0)



##Timepoint x Dose disease state
int1=Dose_disease_state_Timepoint_sigOTUs_alltaxr[!Dose_disease_state_Timepoint_sigOTUs_alltaxr$OTUID %in% Dose_disease_state_sigOTUs_alltaxr,]
int2=Dose_disease_state_sigOTUs_alltaxr[!Dose_disease_state_sigOTUs_alltaxr$OTUID %in% Timepoint_sigOTUs_alltaxr,]




#####looking for otus that are more abundant at dose disease time two than healthy and also than time three...

DoseTime=Dose_disease_state_sigOTUs_alltaxr[Dose_disease_state_sigOTUs_alltaxr$OTUID %in% Dose_disease_state_Timepoint_sigOTUs_alltaxr$OTUID,]

DoseT3Time=DoseTime[DoseTime$OTUID %in% T3_disease_state_Dose_disease_state_sigOTUs_alltaxr$OTUID,]

##Create sig otu otu table
c2norm_DT3T=c2norm[row.names(c2norm) %in% DoseT3Time$OTUID, ]
b_DT3T=t(c2norm_DT3T)
d_DT3T=cbind.data.frame(s1$T3_disease_state, s1$Dose_disease_state, s1$Timepoint, b_DT3T)
colnames(d_DT3T)[1]="T3_disease_state"
colnames(d_DT3T)[2]="Dose_disease_state"
colnames(d_DT3T)[3]="Timepoint"

#check order of samples
check=cbind(row.names(d_DT3T), s1[1:1])


#average diseased and healthy
##T3 disease state
d_DT3T_dh=aggregate(. ~T3_disease_state:Dose_disease_state:Timepoint, data=d_DT3T, FUN=mean)
row.names(d_DT3T_dh)=paste(d_DT3T_dh$Dose_disease_state, d_DT3T_dh$T3_disease_state, d_DT3T_dh$Timepoint, sep="_")
d_DT3T_dh=d_DT3T_dh[,-c(1,2,3)]
d_DT3T_dh=t(d_DT3T_dh)
d_DT3T_dh=as.data.frame(d_DT3T_dh)
d_DT3T_dh_T2D$DvH=d_DT3T_dh$Diseased-d_DT3T_dh$Healthy
UpT2D=subset(d_DT3T_dh, d_DT3T_dh$Diseased_Healthy_two>0)
UpT2Dtax=merge(UpT2D, ts, by.x="row.names", by.y="OTUID")

UpT3D=subset(d_DT3T_dh, d_DT3T_dh$Diseased_Diseased_three>0)
UpT3Dtax=merge(UpT3D, ts, by.x="row.names", by.y="OTUID")


#####looking for otus that are more abundant in T3 disease state DD time 2 and time 3. 
T3ex=T3_disease_state_sigOTUs_alltaxr[!T3_disease_state_sigOTUs_alltaxr$OTUID %in% Dose_Site_sigOTUs_alltaxr$OTUID,]
T3ex1=T3ex[!T3ex$OTUID %in% Genotype_sigOTUs_alltaxr$OTUID,]
T3ex2=T3ex1[!T3ex1$OTUID %in% Site_sigOTUs_alltaxr$OTUID,]
T3ex3=T3ex2[!T3ex2$OTUID %in% Site_Dose_Site_sigOTUs_alltaxr$OTUID,]

##Create sig otu otu table
c2norm_T3ex3=c2norm[row.names(c2norm) %in% T3ex3$OTUID, ]
b_T3ex3=t(c2norm_T3ex3)
d_T3ex3=cbind.data.frame(s1$T3_disease_state, s1$Dose_disease_state, s1$Timepoint, b_T3ex3)
colnames(d_T3ex3)[1]="T3_disease_state"
colnames(d_T3ex3)[2]="Dose_disease_state"
colnames(d_T3ex3)[3]="Timepoint"

#check order of samples
check=cbind(row.names(d_T3ex3), s1[1:1])

#average diseased and healthy
##T3 disease state
d_T3ex3_dh=aggregate(. ~T3_disease_state:Dose_disease_state:Timepoint, data=d_T3ex3, FUN=mean)
row.names(d_T3ex3_dh)=paste(d_T3ex3_dh$Dose_disease_state, d_T3ex3_dh$T3_disease_state, d_T3ex3_dh$Timepoint, sep="_")
d_T3ex3_dh=d_T3ex3_dh[,-c(1,2,3)]
d_T3ex3_dh=t(d_T3ex3_dh)
d_T3ex3_dh=as.data.frame(d_T3ex3_dh)
d_T3ex3_dh$DvH=d_T3ex3_dh$Diseased_Diseased_two-d_T3ex3_dh$Healthy_Healthy_two


UpT2D=subset(d_T3ex3_dh, d_T3ex3_dh$Diseased_Diseased_two>0&d_T3ex3_dh$Diseased_Diseased_three>0&d_T3ex3_dh$Diseased_Diseased_dose>0&d_T3ex3_dh$Healthy_Diseased_two==0)
UpT2Dtax=merge(UpT2D, ts, by.x="row.names", by.y="OTUID")
UpT2Dtax$DvHdose=UpT2Dtax$Diseased_Diseased_dose-UpT2Dtax$Healthy_Healthy_dose
UpT2Dtax$DvHT2=UpT2Dtax$Diseased_Diseased_two-UpT2Dtax$Healthy_Healthy_two
UpT2Dtax$DvHT3=UpT2Dtax$Diseased_Diseased_three-UpT2Dtax$Healthy_Healthy_three
UpT2DtaxUp=subset(UpT2Dtax, UpT2Dtax$DvHdose>0&UpT2Dtax$DvHT2>0&UpT2Dtax$DvHT3>0)
colnames(UpT2DtaxUp)= c("OTUID", "D_D_dose", "H_H_dose", "H_H_T1", "D_D_T3",
                        "D_H_T3", "H_H_T3", "D_D_T2", "D_H_T2", "H_D_T2",
                        "H_H_T2", "DvH", "kingdom_s", "phylum_s", "class_s",
                        "order_s", "family_s", "genus_s", "species_s",               
                        "evalue", "ref", "DvHdose", "DvHT2", "DvHT3")  
UpT2DtaxUp$H_D_T2=NULL
UpT2DtaxUpfs=subset(UpT2DtaxUp, UpT2DtaxUp$family_s=="Pasteurellaceae"|UpT2DtaxUp$family_s=="Francisellaceae")

counts=as.data.frame(table(UpT2DtaxUp$family_s))
counts=counts[order(-counts$Freq),]
otu.long_d = melt(UpT2DtaxUpfs, id = c("OTUID","family_s"), 
                  measure = c("H_H_dose", "D_D_dose", "H_H_T1", "H_H_T2", 
                              "D_H_T2", "D_D_T2", "H_H_T3", "D_H_T3", "D_D_T3"
                              ))
otu.long_d$variable = factor(otu.long_d$variable, levels = c("H_H_dose", "D_D_dose", "H_H_T1", "H_H_T2", 
                                                             "D_H_T2", "D_D_T2", "H_H_T3", "D_H_T3", "D_D_T3"))
otu.long_d$OTUID = factor(otu.long_d$OTUID, levels = otu.long_d$OTUID[order(otu.long_d$family_s)])


pdf(file="Output_files/GLMM_results/Primary_pathogens_fp.pdf", width=5, height=7)
plotgb=ggplot(data=otu.long_d, aes(x=OTUID, y=value ))+
  geom_bar(aes(fill=family_s), stat="identity") + theme_bw() +
  theme(legend.position = "left") +
  ggtitle("Primary pathogens") + 
  facet_grid(variable~., scale="free_y")+ ylim(0,160)
plotgb=plotgb + theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank()) +
  scale_x_discrete("OTUID")+scale_fill_manual(values=wes_palette(n=2, name="Moonrise2"))
plotgb

dev.off()




d_T3_dh=t(d_T3_dh)
d_T3_dh$log2fc=(log2(d_T3_dh$Diseased)-log2(d_T3_dh$Healthy))
table(d_T3_dh$log2fc>2)
table(d_T3_dh$log2fc==Inf& d_T3_dh$Healthy==0.00000000&d_T3_dh$Diseased>2)


T3_and_interaction=merge(T3_disease_state_Dose_disease_state_sigOTUs_alltaxr, Site_sigOTUs_alltaxr, by="OTUID")
T3_and_dose_and_interaction=merge(T3_disease_state_Dose_disease_state_sigOTUs_alltaxr, Dose_disease_state_Dose_Site_sigOTUs_alltaxr, by="OTUID")
T3_dose_exclude2=T3_disease_state_Dose_disease_state_sigOTUs_alltaxr[!T3_disease_state_Dose_disease_state_sigOTUs_alltaxr$OTUID %in% c(Dose_disease_state_Timepoint_sigOTUs_alltaxr$OTUID), ]
T3_dose_exclude2=T3_disease_state_Dose_disease_state_sigOTUs_alltaxr[!T3_disease_state_Dose_disease_state_sigOTUs_alltaxr$OTUID %in% Timepoint_sigOTUs_alltaxr$OTUID, ]

T3_dose_exclude2=Dose_disease_state_sigOTUs_alltaxr[!Dose_disease_state_sigOTUs_alltaxr$OTUID %in% Dose_disease_state_Timepoint_sigOTUs_alltaxr$OTUID, ]


#beneficial
Healthy=T3_and_dose_and_interaction[!T3_and_dose_and_interaction$OTUID 
                                                            %in% c(Dose_disease_state_Timepoint_sigOTUs_alltaxr$OTUID,
                                                                   Timepoint_Dose_Site_sigOTUs_alltaxr$OTUID,
                                                                   Timepoint_sigOTUs_alltaxr$OTUID,
                                                                   Site_sigOTUs_alltaxr$OTUID,
                                                                   Dose_Site_sigOTUs_alltaxr$OTUID), ]


