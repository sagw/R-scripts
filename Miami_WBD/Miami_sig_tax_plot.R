rm(list=ls())
require(ggplot2)
require(plyr)
setwd("~/Documents/Miami_WBD/Miami_analyses/")
tx=read.csv("output/Miami_sig_taxOnly.csv", h=T)
kingdom=as.data.frame(table(tx$kingdom))
subset(tx$phylum, tx$kingdom)

phylum=as.data.frame(table(tx$phylum))
class=as.data.frame(table(tx$class))
order=as.data.frame(table(tx$order))

full = data.frame(Samples= numeric(0), Counts= integer(0), Taxonomy = character(0))

for (i in 1:ncol(b1tcagg10)) {
  fams=data.frame("Samples"=row.names(b1tcagg10), 
                  "Counts"=b1tcagg10[,i], "Taxonomy"=colnames(b1tcagg10)[i])
  full=rbind(full,fams)
}


plot=ggplot(data=tx ) +
  geom_bar(aes(size=Counts,color=Taxonomy)) +       
  scale_colour_discrete(guide=FALSE)+
  scale_size_continuous(range=c(0,16)) + theme_bw()+theme(legend.position = "left") +
  ggtitle("CK4 corals")