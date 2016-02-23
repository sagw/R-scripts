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
geno_beneficial=subset(d_genoT3_dh, d_genoT3_dh$Healthy_Healthy_one_high>0&d_genoT3_dh$Healthy_Healthy_two_high>0&d_genoT3_dh$Healthy_Healthy_three_high>0)
geno_beneficial$DDT3vsHHT1=geno_beneficial$Healthy_Healthy_one_high-geno_beneficial$Diseased_Diseased_three_high
geno_beneficial$DDT3vsHHT3=geno_beneficial$Healthy_Healthy_three_high-geno_beneficial$Diseased_Diseased_three_high
geno_beneficial$DDT2vsHHT2=geno_beneficial$Healthy_Healthy_two_high-geno_beneficial$Diseased_Diseased_two_high
geno_beneficial$DHT2vsHHT2=geno_beneficial$Healthy_Healthy_two_high-geno_beneficial$Diseased_Healthy_two_high
geno_beneficial$DHT3vsHHT3=geno_beneficial$Healthy_Healthy_three_high-geno_beneficial$Diseased_Healthy_three_high
geno_beneficial_tax=merge(geno_beneficial, ts, by.x="row.names", by.y="OTUID")
geno_beneficial2_tax=subset(geno_beneficial_tax, geno_beneficial_tax$DDT3vsHHT3>0&geno_beneficial_tax$DDT2vsHHT2>0&geno_beneficial_tax$DDT3vsHHT1>0)
geno_beneficial2_tax$Healthy_Diseased_two_low=NULL
d_genoT3_dh_se$Healthy_Diseased_two_low=NULL

d_genoT3_dh_se=d_genoT3_dh_se[row.names(d_genoT3_dh_se) %in% geno_beneficial2_tax$Row.names,]

d_genoT3_dh_se_filt=d_genoT3_dh_se[apply(d_genoT3_dh_se, 2, function(x) all(x<50)),]
geno_beneficial2_tax_se=geno_beneficial2_tax[geno_beneficial2_tax$Row.names %in% row.names(d_genoT3_dh_se_filt),]

Geno_Bene=cbind(geno_beneficial2_tax_se[,1], geno_beneficial2_tax_se[,16:29])
colnames(Geno_Bene)[1]="OTUID"
write.csv(Geno_Bene, file="Output_files/GLMM_results/Genotype_Beneficial.csv")


#geno_beneficial2_defensive_tax=subset(geno_beneficial2_tax, geno_beneficial2_tax$DHT2vsHHT2>0&geno_beneficial2_tax$DHT3vsHHT3)

geno_beneficial2_tax_se_endo=geno_beneficial2_tax_se[geno_beneficial2_tax_se$family_s=="Hahellaceae",]
colnames(geno_beneficial2_tax_se_endo)=c("OTUID", "D_D_3_high", "D_D_3_low", "D_D_2_high", "D_D_2_low",
                                         "D_H_3_high", "D_H_3_low", "D_H_2_high", "D_H_2_low", "H_H_1_high",
                                         "H_H_1_low", "H_H_3_high", "H_H_3_low", "H_H_2_high",    
                                         "H_H_2_low", "DDT3vsHHT1", "DDT3vsHHT3", "DDT2vsHHT2",                  
                                         "DHT2vsHHT2", "DHT3vsHHT3", "kingdom_s", "phylum_s", "class_s", 
                                         "order_s", "family_s", "genus_s", "species_s", "evalue", "ref")

endose=d_genoT3_dh_se_filt[row.names(d_genoT3_dh_se_filt) %in% geno_beneficial2_tax_se_endo$OTUID,]
colnames(endose)=c("D_D_3_high", "D_D_3_low", "D_D_2_high", "D_D_2_low",
                                "D_H_3_high", "D_H_3_low", "D_H_2_high", "D_H_2_low", "H_H_1_high",
                                "H_H_1_low", "H_H_3_high", "H_H_3_low", "H_H_2_high",    
                                "H_H_2_low")

endose$OTUID=row.names(endose)

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
  facet_grid(variable~., scale="free_y") + ylim(-1,60)
plotgb=plotgb + theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank()) +
  scale_x_discrete("OTUID")+scale_fill_manual(values=wes_palette(n=2, name="Royal2"))
plotgb=plotgb+geom_errorbar(aes(ymin=value-SE, ymax=value+SE), position=dodge, width=0.001)
plotgb
dev.off()

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

