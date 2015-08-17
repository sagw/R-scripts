rm(list=ls())
setwd("~/Documents/SD1/SD1_analyses/")
c=read.csv("input_files//SD1_OTU_table_97.csv", h=T, row.names=1)
chimera=read.csv("input_files/uchime_chimeras3.csv", h=T)
neg_control=read.csv("Output_files/unfiltered//neg_control_otus.csv", h=T, row.names=1)
c$OTUID<-row.names(c)
c2 <- c[!c$OTUID %in% chimera$OTUID, ]
c3 <- c[!c$OTUID %in% neg_control$OTUID, ]
c3<-c2
c4<-c3[0:(NCOL(c3)-1)]
c5<-t(c4)
rs = rowSums (c4)
use = (rs > 0)
c6 = c4[ use, ]


s=read.csv("input_files//SD1_map.csv", h=T, as.is=T)
b=t(c6)
d=cbind.data.frame(s,b)



t=read.csv("input_files//rep_set_tax_assignments.csv", header=T)


#filter out chloroplast DNA
dtax<-merge(t, c6, by.x="OTUID", by.y="row.names", all.y=T)

dtax1<-dtax[- grep("Chloroplast", dtax$Tax),]
row.names(dtax1)<-dtax1$OTUID
c<-dtax1[5:NCOL(dtax1)]
write.csv(dtax1, file="filtered_otutable.csv")

m<-matrix(1, ncol=263, nrow=1)

filteredc<-rbind(c, m)