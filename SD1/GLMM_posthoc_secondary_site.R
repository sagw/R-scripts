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
T3ex=T3_disease_state_sigOTUs_alltaxr[T3_disease_state_sigOTUs_alltaxr$OTUID %in% Dose_disease_state_Timepoint_Dose_Site_sigOTUs_alltaxr$OTUID,]
T3ex1=T3ex[!T3ex$OTUID %in% Genotype_sigOTUs_alltaxr$OTUID,]
T3ex3=T3ex1


##Create sig otu otu table
c2norm_T3ex3=c2norm[row.names(c2norm) %in% T3ex3$OTUID, ]
b_T3ex3=t(c2norm_T3ex3)
d_T3ex3=cbind.data.frame(s1$T3_disease_state, s1$Dose_disease_state, s1$Timepoint, s1$Dose_Site, b_T3ex3)
colnames(d_T3ex3)[1]="T3_disease_state"
colnames(d_T3ex3)[2]="Dose_disease_state"
colnames(d_T3ex3)[3]="Timepoint"
colnames(d_T3ex3)[4]="Dose_site"

se=function(x) sd(x)/sqrt(length(x))
d_T3ex3_dh_mean=ddply(d_T3ex3, ~Dose_disease_state:T3_disease_state:Timepoint:Dose_site, colwise(mean))
d_T3ex3_dh_se=ddply(d_T3ex3, ~Dose_disease_state:T3_disease_state:Timepoint:Dose_site, colwise(se))
row.names(d_T3ex3_dh_mean)=gsub(":", "_", d_T3ex3_dh_mean[,1])
row.names(d_T3ex3_dh_se)=gsub(":", "_", d_T3ex3_dh_se[,1])
d_T3ex3_dh_mean=d_T3ex3_dh_mean[,-c(1)]
d_T3ex3_dh_se=d_T3ex3_dh_se[,-c(1)]
d_T3ex3_dh_mean=t(d_T3ex3_dh_mean)
d_T3ex3_dh_se=t(d_T3ex3_dh_se)
d_T3ex3_dh_mean=as.data.frame(d_T3ex3_dh_mean)
d_T3ex3_dh_se=as.data.frame(d_T3ex3_dh_se)
d_T3ex3_dh_mean$Healthy_Diseased_two_CK14=NULL
UpT2Dtax=merge(d_T3ex3_dh_mean, ts, by.x="row.names", by.y="OTUID")
d_T3ex3_dh_se$Healthy_Diseased_two_CK14=NULL

colnames(UpT2Dtax)= c("OTUID", "D_D_dose_CK14", "D_D_dose_CK4", 
                      "D_D_T3_CK14", "D_D_T3_CK4", "D_D_T2_CK14", 
                      "D_D_T2_CK4", "D_H_T3_CK14", "D_H_T3_CK4",
                      "D_H_T2_CK14", "D_H_T2_CK4", "H_H_dose_CK14", "H_H_dose_CK4", 
                      "H_H_T1_CK14", "H_H_T1_CK4", "H_H_T3_CK14", "H_H_T3_CK4",
                      "H_H_T2_CK14", "H_H_T2_CK4", "kingdom_s", "phylum_s", "class_s",
                      "order_s", "family_s", "genus_s", "species_s",               
                      "evalue", "ref")  
colnames(d_T3ex3_dh_se)= c("D_D_dose_CK14", "D_D_dose_CK4", 
                           "D_D_T3_CK14", "D_D_T3_CK4", "D_D_T2_CK14", 
                           "D_D_T2_CK4", "D_H_T3_CK14", "D_H_T3_CK4",
                           "D_H_T2_CK14", "D_H_T2_CK4", "H_H_dose_CK14", "H_H_dose_CK4", 
                           "H_H_T1_CK14", "H_H_T1_CK4", "H_H_T3_CK14", "H_H_T3_CK4",
                           "H_H_T2_CK14", "H_H_T2_CK4") 
d_T3ex3_dh_se$OTUID=row.names(d_T3ex3_dh_se)

Secondary=cbind(UpT2Dtax[,1], UpT2Dtax[,18:28])
colnames(Secondary)[1]="OTUID"
write.csv(Primary, file="Output_files/GLMM_results/Primary_pathogens_site.csv")
countso=as.data.frame(table(UpT2Dtax$order_s))
countso=countso[order(-countso$Freq),]
counts=as.data.frame(table(UpT2Dtax$family_s))
counts=counts[order(-counts$Freq),]



counts=as.data.frame(table(UpT3DtaxUpse$family_s))
counts=counts[order(-counts$Freq),]

write.csv(counts, file="Output_files/GLMM_results/secondary_list.csv")

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

pdf(file="Output_files/GLMM_results/Secondary_colonizers_c.pdf", width=6, height=7)
dodge=position_dodge(width=0.9)
plotgb=ggplot(data=otu.long_full, aes(x=OTUID, y=value ))+
  geom_bar(aes(fill=family_s),position=dodge, stat="identity") + theme_bw() +
  theme(legend.position = "left") +
  ggtitle("Secondary colonizers") + 
  facet_grid(variable~., scale="free_y") + ylim(-1,10)
plotgb=plotgb + theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank()) +
  scale_x_discrete("OTUID")+scale_fill_manual(values=wes_palette(n=2, name="GrandBudapest"))
plotgb=plotgb+geom_errorbar(aes(ymin=value-SE, ymax=value+SE), position=dodge, width=0.001)
plotgb

dev.off()
