rm(list=ls())
setwd("~/Documents/Miami_WBD/Miami_97_otus/")
c=read.csv("otu_table_mc2_2.csv", h=T, row.names=1)

chimera=read.table("chimeras3.txt", h=T)
neg_control=read.csv("~/Documents/SD1/SD1_analyses/Output_files/unfiltered//neg_control_otus.csv", h=T, row.names=1)
c$OTUID=row.names(c)
c2 = c[!c$OTUID %in% chimera$X101015, ]
c3 = c[!c2$OTUID %in% neg_control$OTUID, ]
c4=c3[0:(NCOL(c3)-1)]

t=read.csv("tax/rep_set_tax_assignments.csv", header=F)
colnames(t)=c("OTUID", "tax", "evalue", "ref")

write.csv(t, file="rep_set_tax_assignments.csv")


#filter out chloroplast DNA
dtax<-merge(t, c4, by.x="OTUID", by.y="row.names")

dtax1=dtax[- grep("Chloroplast", dtax$tax),]
row.names(dtax1)<-dtax1$OTUID
cfilt=dtax1[5:NCOL(dtax1)]
write.csv(cfilt, file="filtered_otutable.csv")
