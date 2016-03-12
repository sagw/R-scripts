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


#####looking for otus that are more abundant in T3 disease state DD time 2 and time 3. 
T3ex=T3_disease_state_sigOTUs_alltaxr
T3ex2=T3ex[!T3ex$OTUID %in% Dose_Site_sigOTUs_alltaxr$OTUID,]
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
#d_T3ex3_dh$DvH=d_T3ex3_dh$Diseased_Diseased_two-d_T3ex3_dh$Healthy_Healthy_two
UpT2D=subset(d_T3ex3_dh, d_T3ex3_dh$Diseased_Diseased_two>0&d_T3ex3_dh$Diseased_Diseased_three>0&d_T3ex3_dh$Diseased_Diseased_dose==0&d_T3ex3_dh$Healthy_Healthy_one>0)
UpT2Dtax=merge(UpT2D, ts, by.x="row.names", by.y="OTUID")
UpT2Dtax$DvHdose=UpT2Dtax$Diseased_Diseased_dose-UpT2Dtax$Healthy_Healthy_dose
UpT2Dtax$DvHT2=UpT2Dtax$Diseased_Diseased_two-UpT2Dtax$Healthy_Healthy_two
UpT2Dtax$DvHT3=UpT2Dtax$Diseased_Diseased_three-UpT2Dtax$Healthy_Healthy_three
UpT2DtaxUp=subset(UpT2Dtax, UpT2Dtax$DvHT2>0&UpT2Dtax$DvHT3>0)
UpT2DtaxUp$Healthy_Diseased_two=NULL
d_T3ex3_dh_se$Healthy_Diseased_two=NULL
T2D_se=d_T3ex3_dh_se[row.names(d_T3ex3_dh_se) %in% UpT2DtaxUp$Row.names,]
T2D_se=as.matrix(T2D_se)
T2D_se_filt=T2D_se[apply(T2D_se, 2, function(x) all(x<500)),]
UpT2DtaxUpse=UpT2DtaxUp[UpT2DtaxUp$Row.names %in% row.names(T2D_se_filt),]
UpT2DtaxUpse$DvH=NULL

colnames(UpT2DtaxUpse)= c("OTUID", "D_D_dose", "D_D_T3", "D_D_T2", 
                          "D_H_T3", "D_H_T2", "H_H_dose", "H_H_T1", 
                          "H_H_T3", "H_H_T2", "kingdom_s", "phylum_s", "class_s",
                          "order_s", "family_s", "genus_s", "species_s",               
                          "evalue", "ref", "DvHdose", "DvHT2", "DvHT3")  
colnames(T2D_se_filt)= c("D_D_dose", "D_D_T3", "D_D_T2", 
                         "D_H_T3", "D_H_T2", "H_H_dose", "H_H_T1", 
                         "H_H_T3", "H_H_T2") 


Primary=cbind(UpT2DtaxUpse[,1], UpT2DtaxUpse[,11:22])
colnames(Primary)[1]="OTUID"
write.csv(Primary, file="Output_files/GLMM_results/Primary_responders2.csv")
countso=as.data.frame(table(UpT2DtaxUpse$order_s))
countso=countso[order(-countso$Freq),]
counts=as.data.frame(table(UpT2DtaxUpse$family_s))
counts=counts[order(-counts$Freq),]

###pull out families and graph
#Rhodobacteraceae
UpT2DtaxUpr=subset(UpT2DtaxUpse, UpT2DtaxUpse$family_s=="Rhodobacteraceae")
ser=T2D_se_filt[row.names(T2D_se_filt) %in% UpT2DtaxUpr$OTUID,]
ser=as.data.frame(ser)
ser$OTUID=row.names(ser)

otu.long_d = melt(UpT2DtaxUpr, id = c("OTUID","family_s"), 
                  measure = c("H_H_dose", "D_D_dose", "H_H_T1", "H_H_T2", 
                              "D_H_T2", "D_D_T2", "H_H_T3", "D_H_T3", "D_D_T3"
                  ))
otu.long_se = melt(ser, id = c("OTUID"), 
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

pdf(file="Output_files/GLMM_results/plots/Primary_responders_r.pdf", width=6, height=10)
dodge=position_dodge(width=0.9)
plotgb=ggplot(data=otu.long_full, aes(x=OTUID, y=value ))+
  geom_bar(aes(fill=family_s),position=dodge, stat="identity") + theme_bw() +
  theme(legend.position = "left") +
  ggtitle("Primary pathogens") + 
  facet_grid(variable~., scale="free_y") + ylim(-1,15)
plotgb=plotgb + theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank()) +
  scale_x_discrete("OTUID")+scale_fill_manual(values="#E58601")
plotgb=plotgb+geom_errorbar(aes(ymin=value-SE, ymax=value+SE), position=dodge, width=0.001)
plotgb

dev.off()

####Plot cryomorphaceae##################
#separate Cryomorphaceae
UpT2DtaxUpp=subset(UpT2DtaxUpse, UpT2DtaxUpse$family_s=="Cryomorphaceae")
seps=T2D_se_filt[row.names(T2D_se_filt) %in% UpT2DtaxUpp$OTUID,]
seps=as.data.frame(seps)
seps$OTUID=row.names(seps)

#melt dataframes
otu.long_d = melt(UpT2DtaxUpp, id = c("OTUID","family_s"), 
                  measure = c("H_H_dose", "D_D_dose", "H_H_T1", "H_H_T2", 
                              "D_H_T2", "D_D_T2", "H_H_T3", "D_H_T3", "D_D_T3"
                  ))
otu.long_se = melt(seps, id = c("OTUID"), 
                   measure = c("H_H_dose", "D_D_dose", "H_H_T1", "H_H_T2", 
                               "D_H_T2", "D_D_T2", "H_H_T3", "D_H_T3", "D_D_T3"
                   ))
#combine means and ses
colnames(otu.long_se)[3]="SE"
otu.long_se$OTU_group=paste(otu.long_se$OTUID, otu.long_d$variable, sep="_")
otu.long_d$OTU_group=paste(otu.long_d$OTUID, otu.long_d$variable, sep="_")
otu.long_full=merge(otu.long_d, otu.long_se, by="OTU_group")
otu.long_full=otu.long_full[,3:8]
colnames(otu.long_full)=c("family_s", "variable", "value", "OTUID", "variable.y", "SE")

otu.long_full$OTUID = factor(otu.long_full$OTUID, levels = otu.long_full$OTUID[order(otu.long_full$family_s)])
otu.long_full$variable = factor(otu.long_full$variable, levels = c("H_H_dose", "D_D_dose", "H_H_T1", "H_H_T2", 
                                                                   "D_H_T2", "D_D_T2", "H_H_T3", "D_H_T3", "D_D_T3"))

pdf(file="Output_files/GLMM_results/plots/Primary_responders_c.pdf", width=6, height=10)
dodge=position_dodge(width=0.9)
plotgb=ggplot(data=otu.long_full, aes(x=OTUID, y=value ))+
  geom_bar(aes(fill=family_s),position=dodge, stat="identity") + theme_bw() +
  theme(legend.position = "left") +
  ggtitle("Primary pathogens") + 
  facet_grid(variable~., scale="free_y") + ylim(-1,150)
plotgb=plotgb + theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank()) +
  scale_x_discrete("OTUID")+scale_fill_manual(values="#5B1A18")
plotgb=plotgb+geom_errorbar(aes(ymin=value-SE, ymax=value+SE), position=dodge, width=0.001)
plotgb

dev.off()


####Alteromonadaceae
UpT2DtaxUpa=subset(UpT2DtaxUpse, UpT2DtaxUpse$family_s=="Alteromonadaceae")
sea=T2D_se_filt[row.names(T2D_se_filt) %in% UpT2DtaxUpa$OTUID,]
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

pdf(file="Output_files/GLMM_results/plots/Primary_responders_a.pdf", width=6, height=10)
dodge=position_dodge(width=0.9)
plotgb=ggplot(data=otu.long_full, aes(x=OTUID, y=value ))+
  geom_bar(aes(fill=family_s),position=dodge, stat="identity") + theme_bw() +
  theme(legend.position = "left") +
  ggtitle("Primary pathogens") + 
  facet_grid(variable~., scale="free_y") + ylim(-1,75)
plotgb=plotgb + theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank()) +
  scale_x_discrete("OTUID")+scale_fill_manual(values="#E6A0C4")
plotgb=plotgb+geom_errorbar(aes(ymin=value-SE, ymax=value+SE), position=dodge, width=0.001)
plotgb

dev.off()



####Flavobacteriaceae
UpT2DtaxUpf=subset(UpT2DtaxUpse, UpT2DtaxUpse$family_s=="Flavobacteriaceae")
sef=T2D_se_filt[row.names(T2D_se_filt) %in% UpT2DtaxUpf$OTUID,]
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

pdf(file="Output_files/GLMM_results/Primary_responders_fl.pdf", width=6, height=10)
dodge=position_dodge(width=0.9)
plotgb=ggplot(data=otu.long_full, aes(x=OTUID, y=value ))+
  geom_bar(aes(fill=family_s),position=dodge, stat="identity") + theme_bw() +
  theme(legend.position = "left") +
  ggtitle("Primary pathogens") + 
  facet_grid(variable~., scale="free_y") + ylim(-1,70)
plotgb=plotgb + theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank()) +
  scale_x_discrete("OTUID")+scale_fill_manual(values="#78B7C5")
plotgb=plotgb+geom_errorbar(aes(ymin=value-SE, ymax=value+SE), position=dodge, width=0.001)
plotgb

dev.off()


###Primary colonizers not in T3
UpT2Donly=subset(d_T3ex3_dh, d_T3ex3_dh$Diseased_Diseased_two>0&d_T3ex3_dh$Diseased_Diseased_dose>0&d_T3ex3_dh$Healthy_Diseased_two==0)
UpT2Donlytax=merge(UpT2Donly, ts, by.x="row.names", by.y="OTUID")
UpT2Donlytax$DvHdose=UpT2Donlytax$Diseased_Diseased_dose-UpT2Donlytax$Healthy_Healthy_dose
UpT2Donlytax$DvHT2=UpT2Donlytax$Diseased_Diseased_two-UpT2Donlytax$Healthy_Healthy_two
UpT2Donlytax$DvHT3=UpT2Donlytax$Diseased_Diseased_three-UpT2Donlytax$Healthy_Healthy_three
UpT2DonlytaxUp=subset(UpT2Donlytax, UpT2Donlytax$DvHdose>0&UpT2Donlytax$DvHT2>0&UpT2Donlytax$DvHT3<=0&
                        UpT2Donlytax$Healthy_Healthy_dose==0)
UpT2DonlytaxUp$Healthy_Diseased_two=NULL
d_T3ex3_dh_se$Healthy_Diseased_two=NULL
T2D_only_se=d_T3ex3_dh_se[row.names(d_T3ex3_dh_se) %in% UpT2DonlytaxUp$Row.names,]
T2D_only_se=as.matrix(T2D_only_se)
T2D_only_se_filt=T2D_only_se[apply(T2D_only_se, 2, function(x) all(x<500)),]
UpT2DonlytaxUpse=UpT2DonlytaxUp[UpT2DonlytaxUp$Row.names %in% row.names(T2D_only_se_filt),]
UpT2DonlytaxUpse$DvH=NULL

Primary2=cbind(UpT2DonlytaxUpse[,1], UpT2DonlytaxUpse[,11:22])
colnames(Primary2)[1]="OTUID"
write.csv(Primary2, file="Output_files/GLMM_results/Primary_pathogens_2.csv")
countso=as.data.frame(table(UpT2DonlytaxUpse$order_s))
countso=countso[order(-countso$Freq),]
counts=as.data.frame(table(UpT2DonlytaxUpse$family_s))
counts=counts[order(-counts$Freq),]

colnames(UpT2DonlytaxUpse)= c("OTUID", "D_D_dose", "D_D_T3", "D_D_T2", 
                              "D_H_T3", "D_H_T2", "H_H_dose", "H_H_T1", 
                              "H_H_T3", "H_H_T2", "kingdom_s", "phylum_s", "class_s",
                              "order_s", "family_s", "genus_s", "species_s",               
                              "evalue", "ref", "DvHdose", "DvHT2", "DvHT3")  
colnames(T2D_only_se_filt)= c("D_D_dose", "D_D_T3", "D_D_T2", 
                              "D_H_T3", "D_H_T2", "H_H_dose", "H_H_T1", 
                              "H_H_T3", "H_H_T2")
T2D_only_se_filt=as.data.frame(T2D_only_se_filt)
T2D_only_se_filt$OTUID=row.names(T2D_only_se_filt)
###graph
otu.long_d = melt(UpT2DonlytaxUpse, id = c("OTUID","order_s"), 
                  measure = c("H_H_dose", "D_D_dose", "H_H_T1", "H_H_T2", 
                              "D_H_T2", "D_D_T2", "H_H_T3", "D_H_T3", "D_D_T3"
                  ))
otu.long_se = melt(T2D_only_se_filt, id = c("OTUID"), 
                   measure = c("H_H_dose", "D_D_dose", "H_H_T1", "H_H_T2", 
                               "D_H_T2", "D_D_T2", "H_H_T3", "D_H_T3", "D_D_T3"
                   ))
colnames(otu.long_se)[3]="SE"
otu.long_se$OTU_group=paste(otu.long_se$OTUID, otu.long_d$variable, sep="_")
otu.long_d$OTU_group=paste(otu.long_d$OTUID, otu.long_d$variable, sep="_")

otu.long_full=merge(otu.long_d, otu.long_se, by="OTU_group")
otu.long_full=otu.long_full[,3:8]
colnames(otu.long_full)=c("order", "variable", "value", "OTUID", "variable.y", "SE")

otu.long_full$OTUID = factor(otu.long_full$OTUID, levels = otu.long_full$OTUID[order(otu.long_full$order)])
otu.long_full$variable = factor(otu.long_full$variable, levels = c("H_H_dose", "D_D_dose", "H_H_T1", "H_H_T2", 
                                                                   "D_H_T2", "D_D_T2", "H_H_T3", "D_H_T3", "D_D_T3"))

pdf(file="Output_files/GLMM_results/Temporary_primary_colonizers.pdf", width=6, height=7)

dodge=position_dodge(width=0.9)
plotgb=ggplot(data=otu.long_full, aes(x=OTUID, y=value ))+
  geom_bar(aes(fill=order),position=dodge, stat="identity") + theme_bw() +
  theme(legend.position = "left") +
  ggtitle("Temporary primary colonizers") + 
  facet_grid(variable~., scale="free_y") + ylim(-1,40)
plotgb=plotgb + theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank()) +
  scale_x_discrete("OTUID")
plotgb=plotgb+geom_errorbar(aes(ymin=value-SE, ymax=value+SE), position=dodge, width=0.001)
plotgb

dev.off()
