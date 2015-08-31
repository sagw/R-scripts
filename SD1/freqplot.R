rm(list=ls())
setwd("~/Documents/SD1/SD1_analyses/")
c=read.csv("input_files//filtered_otutable.csv", h=T, row.names=1)
c=c[5:NCOL(c)]

#information about samples
s=read.csv("input_files//SD1_map.csv", h=T, as.is=T)
b=t(c)
d=cbind.data.frame(s,b)
d=d[!d$Sample_name=="SD1.CK14D2.CK14.O",]
d=d[d$Timepoint %in% c("one","two","three"),]

check=cbind(row.names(d), d[1:1])

c=d[14:NCOL(d)]
c=t(c)
s=d[0:13]

data=as.data.frame(rowSums(c))
colnames(data)="abund"

data$OTUID=rownames(c)

data=data[with (data, order(-abund)),]
data$rank=1:nrow(data)

freq=apply(c!=0, 1, sum)
freq=as.data.frame(freq)
freq$OTUID=row.names(c)

dataf=merge(data, freq, by.x="OTUID", by.y="OTUID")

dataf=dataf[with (dataf, order(rank)),]

dataf$cumfreq=cumsum(dataf$freq)

plot(dataf$rank, dataf$cumfreq, lwd=0.25, main="Cumulative frequency by otu rank abundance")

abline(a=NULL, b=NULL, v=20000)

filtered=dataf[dataf$rank<20000,]

write.csv(filtered, file=("otus_filtered20000.csv"))
