rm(list=ls())
require("ggplot2")
require("reshape")
require("plyr")
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
Genotype_sigOTUs_alltaxr=read.csv("Output_files/GLMM_results/lists/Genotype_sigOTUs_alltax.csv")
T3_disease_state_sigOTUs_alltaxr=read.csv("Output_files/GLMM_results/lists/T3_disease_state_sigOTUs_alltax.csv", header=TRUE, row.names=1)
Dose_disease_state_sigOTUs_alltaxr=read.csv("Output_files/GLMM_results/lists/Dose_disease_state_sigOTUs_alltax.csv", header=TRUE, row.names=1)
Site_sigOTUs_alltaxr=read.csv("Output_files/GLMM_results/lists/Site_sigOTUs_alltax.csv", header=TRUE, row.names=1)
Timepoint_sigOTUs_alltaxr=read.csv("Output_files/GLMM_results/lists/Timepoint_sigOTUs_alltax.csv", header=TRUE, row.names=1)
Dose_Site_sigOTUs_alltaxr=read.csv("Output_files/GLMM_results/lists/Dose_Site_sigOTUs_alltax.csv", header=TRUE, row.names=1)
Site_Dose_disease_state_sigOTUs_alltaxr=read.csv("Output_files/GLMM_results/lists/Site_Dose_disease_state_sigOTUs_alltax.csv", header=TRUE, row.names=1)
Site_Timepoint_sigOTUs_alltaxr=read.csv("Output_files/GLMM_results/lists/Site_Timepoint_sigOTUs_alltax.csv", header=TRUE, row.names=1)
Dose_disease_state_Timepoint_sigOTUs_alltaxr=read.csv("Output_files/GLMM_results/lists/Dose_disease_state_Timepoint_sigOTUs_alltax.csv", header=TRUE, row.names=1)
Site_Dose_Site_sigOTUs_alltaxr=read.csv("Output_files/GLMM_results/lists/Site_Dose_Site_sigOTUs_alltax.csv", header=TRUE, row.names=1)
Dose_disease_state_Dose_Site_sigOTUs_alltaxr=read.csv("Output_files/GLMM_results/lists/Dose_disease_state_Dose_Site_sigOTUs_alltax.csv", header=TRUE, row.names=1)
Timepoint_Dose_Site_sigOTUs_alltaxr=read.csv("Output_files/GLMM_results/lists/Timepoint_Dose_Site_sigOTUs_alltax.csv", header=TRUE, row.names=1)
Site_Dose_disease_state_Timepoint_sigOTUs_alltaxr=read.csv("Output_files/GLMM_results/lists/Site_Dose_disease_state_Timepoint_sigOTUs_alltax.csv", header=TRUE, row.names=1)
Site_Dose_disease_state_Dose_site_sigOTUs_alltaxr=read.csv("Output_files/GLMM_results/lists/Site_Dose_disease_state_Dose_Site_sigOTUs_alltax.csv", header=TRUE, row.names=1)
Site_Timepoint_Dose_Site_sigOTUs_alltaxr=read.csv("Output_files/GLMM_results/lists/Site_Timepoint_Dose_Site_sigOTUs_alltax.csv", header=TRUE, row.names=1)
Dose_disease_state_Timepoint_Dose_Site_sigOTUs_alltaxr=read.csv("Output_files/GLMM_results/lists/Dose_disease_state_Timepoint_Dose_Site_sigOTUs_alltax.csv", header=TRUE, row.names=1)
Site_Dose_disease_state_Timepoint_Dose_Site_sigOTUs_alltaxr=read.csv("Output_files/GLMM_results/lists/Site_Dose_disease_state_Timepoint_Dose_Site_sigOTUs_alltax.csv", header=TRUE, row.names=1)


#beneficial
#OTUs significant in genotype and T3
Geno_T3=merge(Genotype_sigOTUs_alltaxr, T3_disease_state_sigOTUs_alltaxr, by="OTUID")
#Geno_dose=merge(Genotype_sigOTUs_alltaxr, Dose_disease_state_sigOTUs_alltaxr, by="OTUID")
#Geno_T3_time=merge(Geno_T3, Timepoint_sigOTUs_alltaxr, by="OTUID")

##Find OTUs up in healthy
c2norm_genoT3=c2norm[row.names(c2norm) %in% Geno_T3$OTUID, ]
b_genoT3=t(c2norm_genoT3)
high= c("B_CK4", "P_CK14", "P_CK4", "W_CK4")
s1$highlow=ifelse(s1$Genotype %in% high, "high", "low") 
d_genoT3=cbind.data.frame(s1$Dose_disease_state, s1$T3_disease_state, s1$Timepoint, s1$highlow, b_genoT3)
colnames(d_genoT3)[1]="Dose_disease_state"
colnames(d_genoT3)[2]="T3_disease_state"
colnames(d_genoT3)[3]="Timepoint"
colnames(d_genoT3)[4]="Genotype"
d_genoT3=d_genoT3[d_genoT3$Timepoint %in% c("one","two", "three"),]

#average diseased and healthy
se=function(x) sd(x)/sqrt(length(x))
d_genoT3_dh=ddply(d_genoT3, ~ Dose_disease_state:T3_disease_state:Timepoint:Genotype, colwise(mean))
d_genoT3_dh_se=ddply(d_genoT3, ~ Dose_disease_state:T3_disease_state:Timepoint:Genotype, colwise(se))
row.names(d_genoT3_dh)=gsub(":", "_", d_genoT3_dh[,1])
row.names(d_genoT3_dh_se)=gsub(":", "_", d_genoT3_dh_se[,1])
d_genoT3_dh=d_genoT3_dh[,-c(1)]
d_genoT3_dh_se=d_genoT3_dh_se[,-c(1)]
d_genoT3_dh=t(d_genoT3_dh)
d_genoT3_dh=as.data.frame(d_genoT3_dh)
d_genoT3_dh_se=t(d_genoT3_dh_se)
d_genoT3_dh_se=as.data.frame(d_genoT3_dh_se)
geno_beneficial=subset(d_genoT3_dh, d_genoT3_dh$Healthy_Healthy_one_high>0&
                         d_genoT3_dh$Healthy_Healthy_two_high>0&d_genoT3_dh$Healthy_Healthy_three_high>0|
                         d_genoT3_dh$Healthy_Healthy_one_low>0&d_genoT3_dh$Healthy_Healthy_three_low>0&d_genoT3_dh$Healthy_Healthy_two_low>0)
geno_beneficial$DDT3vsHHT1=geno_beneficial$Healthy_Healthy_one_high-geno_beneficial$Diseased_Diseased_three_high
geno_beneficial$DDT3vsHHT3=geno_beneficial$Healthy_Healthy_three_high-geno_beneficial$Diseased_Diseased_three_high
geno_beneficial$DDT2vsHHT2=geno_beneficial$Healthy_Healthy_two_high-geno_beneficial$Diseased_Diseased_two_high
geno_beneficial$DHT2vsHHT2=geno_beneficial$Healthy_Healthy_two_high-geno_beneficial$Diseased_Healthy_two_high
geno_beneficial$DHT3vsHHT3=geno_beneficial$Healthy_Healthy_three_high-geno_beneficial$Diseased_Healthy_three_high
geno_beneficial$DDT3vsHHT1L=geno_beneficial$Healthy_Healthy_one_low-geno_beneficial$Diseased_Diseased_three_low
geno_beneficial$DDT3vsHHT3L=geno_beneficial$Healthy_Healthy_three_low-geno_beneficial$Diseased_Diseased_three_low
geno_beneficial$DDT2vsHHT2L=geno_beneficial$Healthy_Healthy_two_low-geno_beneficial$Diseased_Diseased_two_low
geno_beneficial$DHT2vsHHT2L=geno_beneficial$Healthy_Healthy_two_low-geno_beneficial$Diseased_Healthy_two_low
geno_beneficial$DHT3vsHHT3L=geno_beneficial$Healthy_Healthy_three_low-geno_beneficial$Diseased_Healthy_three_low


geno_beneficial_tax=merge(geno_beneficial, ts, by.x="row.names", by.y="OTUID")
geno_beneficial2_tax=subset(geno_beneficial_tax, geno_beneficial_tax$DDT3vsHHT3>0&geno_beneficial_tax$DDT2vsHHT2>0&
                              geno_beneficial_tax$DDT3vsHHT1>0|geno_beneficial_tax$DDT3vsHHT3L>0&
                              geno_beneficial_tax$DDT2vsHHT2L>0&geno_beneficial_tax$DDT3vsHHT1L>0)

geno_beneficial2_tax$Healthy_Diseased_two_low=NULL
d_genoT3_dh_se$Healthy_Diseased_two_low=NULL

d_genoT3_dh_se=d_genoT3_dh_se[row.names(d_genoT3_dh_se) %in% geno_beneficial2_tax$Row.names,]

d_genoT3_dh_se_filt=d_genoT3_dh_se[apply(d_genoT3_dh_se, 2, function(x) all(x<50)),]
geno_beneficial2_tax_se=geno_beneficial2_tax[geno_beneficial2_tax$Row.names %in% row.names(d_genoT3_dh_se_filt),]

counts=as.data.frame(table(geno_beneficial2_tax_se$family_s))
counts=counts[order(-counts$Freq),]

#defensive otus
geno_beneficial2_tax_def=subset(geno_beneficial2_tax_se, geno_beneficial2_tax_se$DHT2vsHHT2<0& geno_beneficial2_tax_se$DHT3vsHHT3<0)
d_genoT3_def_se=d_genoT3_dh_se[row.names(d_genoT3_dh_se) %in% geno_beneficial2_tax_def$Row.names,]


#save list of beneficial otus
Geno_Bene=cbind(geno_beneficial2_tax_se[,1], geno_beneficial2_tax_se[,16:34])
colnames(Geno_Bene)[1]="OTUID"
write.csv(Geno_Bene, file="Output_files/GLMM_results/Genotype_Beneficial.csv")

#save list of defensive beneficial otus
Geno_Bene_def=cbind(geno_beneficial2_tax_def[,1], geno_beneficial2_tax_def[,16:34])
colnames(Geno_Bene_def)[1]="OTUID"
countsdef=as.data.frame(table(Geno_Bene_def$family_s))
countsdef=countsdef[order(-countsdef$Freq),]
write.csv(Geno_Bene_def, file="Output_files/GLMM_results/Genotype_Beneficial_defensive.csv")

#isolate just endos
geno_beneficial2_tax_se_endo=geno_beneficial2_tax_se[geno_beneficial2_tax_se$family_s=="Hahellaceae",]
geno_beneficial2_tax_se_endo=cbind(geno_beneficial2_tax_se_endo[,1:15], geno_beneficial2_tax_se_endo[26:NCOL(geno_beneficial2_tax_se_endo)])
colnames(geno_beneficial2_tax_se_endo)=c("OTUID", "D_D_3_high", "D_D_3_low", "D_D_2_high", "D_D_2_low",
                                         "D_H_3_high", "D_H_3_low", "D_H_2_high", "D_H_2_low", "H_H_1_high",
                                         "H_H_1_low", "H_H_3_high", "H_H_3_low", "H_H_2_high",    
                                         "H_H_2_low", "kingdom_s", "phylum_s", "class_s", 
                                         "order_s", "family_s", "genus_s", "species_s", "evalue", "ref")

endose=d_genoT3_dh_se_filt[row.names(d_genoT3_dh_se_filt) %in% geno_beneficial2_tax_se_endo$OTUID,]
colnames(endose)=c("D_D_3_high", "D_D_3_low", "D_D_2_high", "D_D_2_low",
                                "D_H_3_high", "D_H_3_low", "D_H_2_high", "D_H_2_low", "H_H_1_high",
                                "H_H_1_low", "H_H_3_high", "H_H_3_low", "H_H_2_high",    
                                "H_H_2_low")

endose$OTUID=row.names(endose)

#isolate defensive endos
endo_def=geno_beneficial2_tax_se_endo[geno_beneficial2_tax_se_endo$OTUID %in% geno_beneficial2_tax_def$Row.names,]
endo_def_se=endose[endose$OTUID %in% row.names(d_genoT3_def_se),]


#plot all beneficial OTUs (endozoicomonas)
otu.long_h = melt(geno_beneficial2_tax_se_endo, id = c("OTUID","family_s"), 
                  measure = c("H_H_1_high", "H_H_2_high", "D_H_2_high", "D_D_2_high",
                              "H_H_3_high", "D_H_3_high", "D_D_3_high", "H_H_1_low",
                              "H_H_2_low", "D_H_2_low", "D_D_2_low", "H_H_3_low",
                              "D_H_3_low", "D_D_3_low"))
otu.long_se = melt(endose, id = c("OTUID"), 
                   measure = c("H_H_1_high", "H_H_2_high", "D_H_2_high", "D_D_2_high",
                                 "H_H_3_high", "D_H_3_high", "D_D_3_high", "H_H_1_low",
                                 "H_H_2_low", "D_H_2_low", "D_D_2_low", "H_H_3_low",
                                 "D_H_3_low", "D_D_3_low"))
colnames(otu.long_se)[3]="SE"
otu.long_se$OTU_group=paste(otu.long_se$OTUID, otu.long_h$variable, sep="_")

otu.long_h$OTU_group=paste(otu.long_h$OTUID, otu.long_h$variable, sep="_")

otu.long_full=merge(otu.long_h, otu.long_se, by="OTU_group")
otu.long_full=otu.long_full[,3:8]
colnames(otu.long_full)=c("family_s", "variable", "value", "OTUID", "variable.y", "SE")

otu.long_full$OTUID = factor(otu.long_full$OTUID, levels = otu.long_full$OTUID[order(otu.long_full$family_s)])
otu.long_full$variable = factor(otu.long_full$variable, levels = c("H_H_1_high", "H_H_2_high", "D_H_2_high", "D_D_2_high",
                                                                   "H_H_3_high", "D_H_3_high", "D_D_3_high", "H_H_1_low",
                                                                   "H_H_2_low", "D_H_2_low", "D_D_2_low", "H_H_3_low",
                                                                   "D_H_3_low", "D_D_3_low"))


pdf(file="Output_files/GLMM_results/Beneficial_endo2.pdf", width=8, height=15)
dodge=position_dodge(width=0.9)
plotgb=ggplot(data=otu.long_full, aes(x=OTUID, y=value ))+
  geom_bar(aes(fill=family_s),position=dodge, stat="identity") + theme_bw() +
  theme(legend.position = "left") +
  ggtitle("Endozoicomonas: genotype specific healthy symbionts") + 
  facet_grid(variable~., scale="free_y") + ylim(-1,70)
plotgb=plotgb + theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank()) +
  scale_x_discrete("OTUID")+scale_fill_manual(values="#9A8822"))
plotgb=plotgb+geom_errorbar(aes(ymin=value-SE, ymax=value+SE), position=dodge, width=0.001)
plotgb
dev.off()


#plot defensive beneficial OTUs (endozoicomonas)
otu.long_h = melt(endo_def, id = c("OTUID","family_s"), 
                  measure = c("H_H_2_high", "D_H_2_high", 
                              "H_H_3_high", "D_H_3_high", 
                              "H_H_2_low", "D_H_2_low", "H_H_3_low",
                              "D_H_3_low"))
otu.long_se = melt(endo_def_se, id = c("OTUID"), 
                   measure = c("H_H_2_high", "D_H_2_high", 
                               "H_H_3_high", "D_H_3_high", 
                               "H_H_2_low", "D_H_2_low", "H_H_3_low",
                               "D_H_3_low"))
colnames(otu.long_se)[3]="SE"
otu.long_se$OTU_group=paste(otu.long_se$OTUID, otu.long_h$variable, sep="_")

otu.long_h$OTU_group=paste(otu.long_h$OTUID, otu.long_h$variable, sep="_")

otu.long_full=merge(otu.long_h, otu.long_se, by="OTU_group")
otu.long_full=otu.long_full[,3:8]
colnames(otu.long_full)=c("family_s", "variable", "value", "OTUID", "variable.y", "SE")

otu.long_full$OTUID = factor(otu.long_full$OTUID, levels = otu.long_full$OTUID[order(otu.long_full$family_s)])
otu.long_full$variable = factor(otu.long_full$variable, levels = c("H_H_2_high", "D_H_2_high", 
                                                                   "H_H_3_high", "D_H_3_high", 
                                                                   "H_H_2_low", "D_H_2_low", "H_H_3_low",
                                                                   "D_H_3_low"))


pdf(file="Output_files/GLMM_results/Beneficial_endo_defensive.pdf", width=8, height=10)
dodge=position_dodge(width=0.9)
plotgb=ggplot(data=otu.long_full, aes(x=OTUID, y=value ))+
  geom_bar(aes(fill=family_s),position=dodge, stat="identity") + theme_bw() +
  theme(legend.position = "left") +
  ggtitle("Endozoicomonas: genotype specific defensive symbionts") + 
  facet_grid(variable~., scale="free_y") + ylim(-1,62)
plotgb=plotgb + theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank()) +
  scale_x_discrete("OTUID")+scale_fill_manual(values=wes_palette(n=2, name="Royal2"))
plotgb=plotgb+geom_errorbar(aes(ymin=value-SE, ymax=value+SE), position=dodge, width=0.001)
plotgb
dev.off()

######divide T1 sig otus by genotype
c2norm_genoT1=c2norm[row.names(c2norm) %in% geno_beneficial2_tax_se$Row.names, ]
b_genoT1=t(c2norm_genoT1)
d_genoT1=cbind.data.frame(s1$Genotype, s1$Timepoint, b_genoT1)
colnames(d_genoT1)[1]="Genotype"
colnames(d_genoT1)[2]="Timepoint"

d_genoT1=d_genoT1[d_genoT1$Timepoint %in% c("one"),]

#average by genotype
se=function(x) sd(x)/sqrt(length(x))
d_genoT1$Timepoint=NULL
d_genoT1_geno=ddply(d_genoT1, ~ Genotype, colwise(mean))
d_genoT1_geno_se=ddply(d_genoT1, ~ Genotype, colwise(se))
row.names(d_genoT1_geno)=d_genoT1_geno[,1]
d_genoT1_geno=d_genoT1_geno[,-c(1)]
row.names(d_genoT1_geno_se)=d_genoT1_geno_se[,1]
d_genoT1_geno_se=d_genoT1_geno_se[,-c(1)]
d_genoT1_geno=t(d_genoT1_geno)
d_genoT1_geno=as.data.frame(d_genoT1_geno)
d_genoT1_geno_se=t(d_genoT1_geno_se)
d_genoT1_geno_se=as.data.frame(d_genoT1_geno_se)

Geno_g_tax=merge(d_genoT1_geno, ts, by.x="row.names", by.y="OTUID")
Geno_g_tax_endo=Geno_g_tax[Geno_g_tax$family_s=="Hahellaceae",]
endo_g_se=d_genoT1_geno_se[row.names(d_genoT1_geno_se) %in% Geno_g_tax_endo$OTUID,]
endo_g_se$OTUID=row.names(endo_g_se)
colnames(Geno_g_tax_endo)[1]="OTUID"

#plot all beneficial OTUs (endozoicomonas)
otu.long_h = melt(Geno_g_tax_endo, id = c("OTUID","family_s"), 
                  measure = c("B_CK14", "B_CK4", "G_CK14", "G_CK4", "O_CK14", "O_CK4", 
                              "P_CK14", "P_CK4", "W_CK14", "W_CK4" ))
otu.long_se = melt(endo_g_se, id = c("OTUID"), 
                   measure = c("B_CK14", "B_CK4", "G_CK14", "G_CK4", "O_CK14", "O_CK4", 
                               "P_CK14", "P_CK4", "W_CK14", "W_CK4" ))
colnames(otu.long_se)[3]="SE"
otu.long_se$OTU_group=paste(otu.long_se$OTUID, otu.long_h$variable, sep="_")

otu.long_h$OTU_group=paste(otu.long_h$OTUID, otu.long_h$variable, sep="_")

otu.long_full=merge(otu.long_h, otu.long_se, by="OTU_group")
otu.long_full=otu.long_full[,3:8]
colnames(otu.long_full)=c("family_s", "variable", "value", "OTUID", "variable.y", "SE")

otu.long_full$OTUID = factor(otu.long_full$OTUID, levels = otu.long_full$OTUID[order(otu.long_full$family_s)])
otu.long_full$variable = factor(otu.long_full$variable, levels = c("H_H_1_high", "H_H_2_high", "D_H_2_high", "D_D_2_high",
                                                                   "H_H_3_high", "D_H_3_high", "D_D_3_high", "H_H_1_low",
                                                                   "H_H_2_low", "D_H_2_low", "D_D_2_low", "H_H_3_low",
                                                                   "D_H_3_low", "D_D_3_low"))


pdf(file="Output_files/GLMM_results/Beneficial_endo_T1geno.pdf", width=15, height=8)
dodge=position_dodge(width=0.9)
plotgb=ggplot(data=otu.long_full, aes(x=OTUID, y=value ))+
  geom_bar(aes(fill=family_s),position=dodge, stat="identity") + theme_bw() +
  theme(legend.position = "left") +
  ggtitle("Endozoicomonas abundance at T1") + 
  facet_grid(.~variable, scale="free_y") + ylim(-1,100)
plotgb=plotgb + theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank()) +
  scale_x_discrete("OTUID")+scale_fill_manual(values=wes_palette(n=2, name="Royal2"))
plotgb=plotgb+geom_errorbar(aes(ymin=value-SE, ymax=value+SE), position=dodge, width=0.001)
plotgb
dev.off()
