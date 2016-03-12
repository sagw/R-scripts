rm(list=ls())
require("ggplot2")
require("reshape")
library(wesanderson)
require("plyr")

setwd("~/Documents/SD1/SD1_analyses/")

ts=read.csv("input_files/rep_set_tax_assignments_silva_clean2.csv")
colnames(ts)=c("OTUID", "kingdom_s", "phylum_s", "class_s", "order_s", "family_s", "genus_s", "species_s", "evalue", "ref")

tgg=read.csv("input_files/rep_set_tax_assignments_split2.csv")

s1=read.csv("input_files/SD1_map_clean.csv", header=TRUE, row.names=1)
s1=as.data.frame(s1)
c2norm=read.csv("input_files/c2norm.csv", header=TRUE, row.names=1)

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


#####looking for otus that are more abundant in T3 disease state DD time 3 not time 2. 
T3ex=T3_disease_state_sigOTUs_alltaxr[T3_disease_state_sigOTUs_alltaxr$OTUID %in% Dose_disease_state_Timepoint_sigOTUs_alltaxr$OTUID,]
T3ex1=T3ex[!T3ex$OTUID %in% Genotype_sigOTUs_alltaxr$OTUID,]
T3ex3=T3ex

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
se=function(x) sd(x)/sqrt(length(x))
d_T3ex3_dh_mean=ddply(d_T3ex3, ~Dose_disease_state:T3_disease_state:Timepoint, colwise(mean))
d_T3ex3_dh_se=ddply(d_T3ex3, ~Dose_disease_state:T3_disease_state:Timepoint, colwise(se))
row.names(d_T3ex3_dh_mean)=gsub(":", "_", d_T3ex3_dh_mean[,1])
row.names(d_T3ex3_dh_se)=gsub(":", "_", d_T3ex3_dh_se[,1])
d_T3ex3_dh_mean=d_T3ex3_dh_mean[,-c(1)]
d_T3ex3_dh_se=d_T3ex3_dh_se[,-c(1)]
d_T3ex3_dh_mean=t(d_T3ex3_dh_mean)
d_T3ex3_dh_se=t(d_T3ex3_dh_se)
d_T3ex3_dh_mean=as.data.frame(d_T3ex3_dh_mean)
d_T3ex3_dh_se=as.data.frame(d_T3ex3_dh_se)
d_T3ex3_dh=d_T3ex3_dh_mean
UpT3D=subset(d_T3ex3_dh, d_T3ex3_dh$Diseased_Diseased_two==0&d_T3ex3_dh$Diseased_Diseased_three>0&d_T3ex3_dh$Diseased_Diseased_dose==0&d_T3ex3_dh$Healthy_Healthy_one>0)
UpT3Dtax=merge(UpT3D, ts, by.x="row.names", by.y="OTUID")
UpT3Dtax$DvHdose=UpT3Dtax$Diseased_Diseased_dose-UpT3Dtax$Healthy_Healthy_dose
UpT3Dtax$DvHT2=UpT3Dtax$Diseased_Diseased_two-UpT3Dtax$Healthy_Healthy_two
UpT3Dtax$DvHT3=UpT3Dtax$Diseased_Diseased_three-UpT3Dtax$Healthy_Healthy_three
UpT3DtaxUp=subset(UpT3Dtax, UpT3Dtax$DvHT3>0)
UpT3DtaxUp$Healthy_Diseased_two=NULL
d_T3ex3_dh_se$Healthy_Diseased_two=NULL
d_T3ex3_dh_se=d_T3ex3_dh_se[row.names(d_T3ex3_dh_se) %in% UpT3DtaxUp$Row.names,]
d_T3ex3_dh_se=as.matrix(d_T3ex3_dh_se)
d_T3ex3_dh_se_filt=d_T3ex3_dh_se[apply(d_T3ex3_dh_se, 2, function(x) all(x<50)),]
UpT3DtaxUpse=UpT3DtaxUp[UpT3DtaxUp$Row.names %in% row.names(d_T3ex3_dh_se_filt),]

Secondary=cbind(UpT3DtaxUpse[,1], UpT3DtaxUpse[,11:22])
colnames(Secondary)[1]="OTUID"
write.csv(Secondary, file="Output_files/GLMM_results/Secondary_responders2.csv")


colnames(UpT3DtaxUpse)= c("OTUID", "D_D_dose", "D_D_T3", "D_D_T2", 
                          "D_H_T3", "D_H_T2", "H_H_dose", "H_H_T1", 
                          "H_H_T3", "H_H_T2", "DvH", "kingdom_s", "phylum_s", "class_s",
                          "order_s", "family_s", "genus_s", "species_s",               
                          "evalue", "ref", "DvHdose", "DvHT2", "DvHT3")  
colnames(d_T3ex3_dh_se_filt)= c("D_D_dose", "D_D_T3", "D_D_T2", 
                                "D_H_T3", "D_H_T2", "H_H_dose", "H_H_T1", 
                                "H_H_T3", "H_H_T2") 



counts=as.data.frame(table(UpT3DtaxUpse$family_s))
counts=counts[order(-counts$Freq),]

##graph campylobacteraceae

####Campylobacter
UpT3DtaxUpc=subset(UpT3DtaxUpse, UpT3DtaxUpse$family_s=="Campylobacteraceae")
sec=d_T3ex3_dh_se_filt[row.names(d_T3ex3_dh_se_filt) %in% UpT3DtaxUpc$OTUID,]
sec=as.data.frame(sec)
sec$OTUID=row.names(sec)

otu.long_d = melt(UpT3DtaxUpc, id = c("OTUID","family_s"), 
                  measure = c("H_H_dose", "D_D_dose", "H_H_T1", "H_H_T2", 
                              "D_H_T2", "D_D_T2", "H_H_T3", "D_H_T3", "D_D_T3"
                  ))
otu.long_se = melt(sec, id = c("OTUID"), 
                   measure = c("H_H_dose", "D_D_dose", "H_H_T1", "H_H_T2", 
                               "D_H_T2", "D_D_T2", "H_H_T3", "D_H_T3", "D_D_T3"
                   ))
colnames(otu.long_se)[3]="SE"
otu.long_se$OTU_group=paste(otu.long_se$OTUID, otu.long_d$variable, sep="_")

otu.long_d$OTU_group=paste(otu.long_d$OTUID, otu.long_d$variable, sep="_")

otu.long_full=merge(otu.long_d, otu.long_se, by="OTU_group")
otu.long_full=otu.long_full[,3:8]
colnames(otu.long_full)=c("family_s", "variable", "value", "OTUID", "variable.y", "SE")

otu.long_full$OTUID = factor(otu.long_full$OTUID, levels = otu.long_full$OTUID[order(otu.long_full$family_s)])
otu.long_full$variable = factor(otu.long_full$variable, levels = c("H_H_dose", "D_D_dose", "H_H_T1", "H_H_T2", 
                                                                   "D_H_T2", "D_D_T2", "H_H_T3", "D_H_T3", "D_D_T3"))

pdf(file="Output_files/GLMM_results/Secondary_responders_c.pdf", width=6, height=7)
dodge=position_dodge(width=0.9)
plotgb=ggplot(data=otu.long_full, aes(x=OTUID, y=value ))+
  geom_bar(aes(fill=family_s),position=dodge, stat="identity") + theme_bw() +
  theme(legend.position = "left") +
  ggtitle("Secondary responders") + 
  facet_grid(variable~., scale="free_y") + ylim(-0.001,3)
plotgb=plotgb + theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank()) +
  scale_x_discrete("OTUID")+scale_fill_manual(values=wes_palette(n=2, name="GrandBudapest"))
plotgb=plotgb+geom_errorbar(aes(ymin=value-SE, ymax=value+SE), position=dodge, width=0.001)
plotgb

dev.off()

####Flavobacteriaceae
UpT2DtaxUpf=subset(UpT3DtaxUpse, UpT3DtaxUpse$family_s=="Flavobacteriaceae")
sef=d_T3ex3_dh_se_filt[row.names(d_T3ex3_dh_se_filt) %in% UpT2DtaxUpf$OTUID,]
sef=as.data.frame(sef)
sef$OTUID=row.names(sef)

otu.long_d = melt(UpT2DtaxUpf, id = c("OTUID","family_s"), 
                  measure = c("H_H_dose", "D_D_dose", "H_H_T1", "H_H_T2", 
                              "D_H_T2", "D_D_T2", "H_H_T3", "D_H_T3", "D_D_T3"
                  ))
otu.long_se = melt(sef, id = c("OTUID"), 
                   measure = c("H_H_dose", "D_D_dose", "H_H_T1", "H_H_T2", 
                               "D_H_T2", "D_D_T2", "H_H_T3", "D_H_T3", "D_D_T3"
                   ))
colnames(otu.long_se)[3]="SE"
otu.long_se$OTU_group=paste(otu.long_se$OTUID, otu.long_d$variable, sep="_")

otu.long_d$OTU_group=paste(otu.long_d$OTUID, otu.long_d$variable, sep="_")

otu.long_full=merge(otu.long_d, otu.long_se, by="OTU_group")
otu.long_full=otu.long_full[,3:8]
colnames(otu.long_full)=c("family_s", "variable", "value", "OTUID", "variable.y", "SE")

otu.long_full$OTUID = factor(otu.long_full$OTUID, levels = otu.long_full$OTUID[order(otu.long_full$family_s)])
otu.long_full$variable = factor(otu.long_full$variable, levels = c("H_H_dose", "D_D_dose", "H_H_T1", "H_H_T2", 
                                                                   "D_H_T2", "D_D_T2", "H_H_T3", "D_H_T3", "D_D_T3"))

pdf(file="Output_files/GLMM_results/Secondary_responders_fl.pdf", width=6, height=7)
dodge=position_dodge(width=0.9)
plotgb=ggplot(data=otu.long_full, aes(x=OTUID, y=value ))+
  geom_bar(aes(fill=family_s),position=dodge, stat="identity") + theme_bw() +
  theme(legend.position = "left") +
  ggtitle("Primary pathogens") + 
  facet_grid(variable~., scale="free_y") + ylim(-0.001,1)
plotgb=plotgb + theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank()) +
  scale_x_discrete("OTUID")+scale_fill_manual(values="#78B7C5")
plotgb=plotgb+geom_errorbar(aes(ymin=value-SE, ymax=value+SE), position=dodge, width=0.001)
plotgb

dev.off()

####Alteromonadaceae
UpT2DtaxUpa=subset(UpT3DtaxUpse, UpT3DtaxUpse$family_s=="Alteromonadaceae")
sea=d_T3ex3_dh_se_filt[row.names(d_T3ex3_dh_se_filt) %in% UpT2DtaxUpa$OTUID,]
sea=as.data.frame(sea)
sea$OTUID=row.names(sea)

otu.long_d = melt(UpT2DtaxUpa, id = c("OTUID","family_s"), 
                  measure = c("H_H_dose", "D_D_dose", "H_H_T1", "H_H_T2", 
                              "D_H_T2", "D_D_T2", "H_H_T3", "D_H_T3", "D_D_T3"
                  ))
otu.long_se = melt(sea, id = c("OTUID"), 
                   measure = c("H_H_dose", "D_D_dose", "H_H_T1", "H_H_T2", 
                               "D_H_T2", "D_D_T2", "H_H_T3", "D_H_T3", "D_D_T3"
                   ))
colnames(otu.long_se)[3]="SE"
otu.long_se$OTU_group=paste(otu.long_se$OTUID, otu.long_d$variable, sep="_")

otu.long_d$OTU_group=paste(otu.long_d$OTUID, otu.long_d$variable, sep="_")

otu.long_full=merge(otu.long_d, otu.long_se, by="OTU_group")
otu.long_full=otu.long_full[,3:8]
colnames(otu.long_full)=c("family_s", "variable", "value", "OTUID", "variable.y", "SE")

otu.long_full$OTUID = factor(otu.long_full$OTUID, levels = otu.long_full$OTUID[order(otu.long_full$family_s)])
otu.long_full$variable = factor(otu.long_full$variable, levels = c("H_H_dose", "D_D_dose", "H_H_T1", "H_H_T2", 
                                                                   "D_H_T2", "D_D_T2", "H_H_T3", "D_H_T3", "D_D_T3"))

pdf(file="Output_files/GLMM_results/Secondary_responders_a.pdf", width=6, height=7)
dodge=position_dodge(width=0.9)
plotgb=ggplot(data=otu.long_full, aes(x=OTUID, y=value ))+
  geom_bar(aes(fill=family_s),position=dodge, stat="identity") + theme_bw() +
  theme(legend.position = "left") +
  ggtitle("Secondary responders") + 
  facet_grid(variable~., scale="free_y") + ylim(-0.001,3)
plotgb=plotgb + theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank()) +
  scale_x_discrete("OTUID")+scale_fill_manual(values=wes_palette(n=2, name="GrandBudapest2"))
plotgb=plotgb+geom_errorbar(aes(ymin=value-SE, ymax=value+SE), position=dodge, width=0.001)
plotgb

dev.off()

