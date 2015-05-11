rm(list=ls())
require(vegan)
require(DESeq2)
setwd("~/Documents/SD1/SD1_analyses/")
c=read.csv("input_files//SD1_OTU_table_97.csv", h=T, row.names=1)

#information about samples
s=read.csv("input_files//SD1_map.csv", h=T, as.is=T)
b=t(c)
d=cbind.data.frame(s,b)

#check order of samples
check<-cbind(row.names(d), d[1:1])

#separate controls
control<-subset(d, Dose=="control")

#only OTUs in controls
cont<-control[10:NCOL(control)]
cont<-t(cont)
rs = rowSums (cont)
use = (rs > 0)
cont = cont [ use, ]

cont<-t(cont)
contd<-cbind(control[0:9], cont)

#add taxonomy
t=read.csv("input_files//rep_set_tax_assignments.csv", header=T)
contd<-t(contd)
conttax<-merge(t, contd, by.x="OTUID", by.y="row.names", all.y=TRUE)
write.csv(as.data.frame(conttax), file="neg_control_otus.csv")
