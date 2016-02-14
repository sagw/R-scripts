rm(list=ls())
require(ggplot2)
setwd("~/Documents/Miami_WBD/Miami_analyses/")
c=read.csv("input/filtered_otutable.csv", h=T, row.names=1)

#information about samples
s=read.csv("input/Miami_map.csv", h=T, as.is=T)
b=t(c)
d=cbind.data.frame(s,b)
#check order of samples
check=cbind(row.names(d), d[1:1])

t=read.csv("input/rep_set_tax_assignments_clean.csv", header=F)
colnames(t)=c("OTUID", "kingdom", "phylum", "class", "order", "family", "genus", "species", "evalue", "ref")

#convert to percent relative abundance
c1=t(c)
c1per=(c1/rowSums(c1))*100
b1=t(c1per)

#simplify taxonomy
tcof=t[,c("OTUID","class", "order", "family")]

#merge tax with counts
b1t=merge(b1, tcof, by.x="row.names", by.y="OTUID")

#Divide taxonomy into class, order, and family
b1tc=b1t
b1tc$family=NULL
b1tc$order=NULL
b1tc$Row.names=NULL

b1to=b1t
b1to$family=NULL
b1to$class=NULL
b1to$Row.names=NULL

b1tf=b1t
b1tf$class=NULL
b1tf$order=NULL
b1tf$Row.names=NULL
#aggregate by class family or order
b1tcagg=aggregate(. ~ class, data=b1tc, FUN=sum)
b1tfagg=aggregate(. ~ family, data=b1tf, FUN=sum)
b1toagg=aggregate(. ~ order, data=b1to, FUN=sum)

##subset only those families or orders that are found in at least 10% of one coral
b1tcagg10=b1tcagg[apply(b1tcagg[,2:NCOL(b1tcagg)], MARGIN=1, function(x) any(x>10)),]
b1tfagg10=b1tfagg[apply(b1tfagg[,2:NCOL(b1tfagg)], MARGIN=1, function(x) any(x>1)),]
b1toagg10=b1toagg[apply(b1toagg[,2:NCOL(b1toagg)], MARGIN=1, function(x) any(x>10)),]

b1toagg10$order=as.character(b1toagg10$order)
b1toagg10[1,1]="None"

#plot class########################
row.names(b1tcagg10)=b1tcagg10$class
b1tcagg10=b1tcagg10[2:NCOL(b1tcagg10)]
b1tcagg10=t(b1tcagg10)
b1tcagg10=as.data.frame(b1tcagg10)

cfull = data.frame(Samples= numeric(0), Counts= integer(0), Taxonomy = character(0))

for (i in 1:ncol(b1tcagg10)) {
  fams=data.frame("Samples"=row.names(b1tcagg10), 
                  "Counts"=b1tcagg10[,i], "Taxonomy"=colnames(b1tcagg10)[i])
  cfull=rbind(cfull,fams)
}

#plot order#####################################
row.names(b1toagg10)=b1toagg10$order
b1toagg10=b1toagg10[2:NCOL(b1toagg10)]
b1toagg10=t(b1toagg10)

ofull = data.frame(Samples= numeric(0), Counts= integer(0), Taxonomy = character(0))

for (i in 1:ncol(b1toagg10)) {
  fams=data.frame("Samples"=row.names(b1toagg10), 
                  "Counts"=b1toagg10[,i], "Taxonomy"=colnames(b1toagg10)[i])
  ofull=rbind(ofull,fams)
}


##plot family#####################
row.names(b1tfagg10)=b1tfagg10$family
b1tfagg10=b1tfagg10[2:NCOL(b1tfagg10)]
b1tfagg10=t(b1tfagg10)

ffull = data.frame(Samples= numeric(0), Counts= integer(0), Taxonomy = character(0))

for (i in 1:ncol(b1tfagg10)) {
  fams=data.frame("Samples"=row.names(b1tfagg10), 
                  "Counts"=b1tfagg10[,i], "Taxonomy"=colnames(b1tfagg10)[i])
  ffull=rbind(ffull,fams)
}



##plot###############################
fulldat=merge(ofull, s, by.x="Samples", by.y="Sample_name", all.x=TRUE)

ggplotColours =function(n=6, h=c(0, 360) +15){
  if ((diff(h)%%360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}



##plot all samples 
pdf(file="output/Order_taxplots.pdf", width=8, height=5)
plotorderalltime=ggplot(data=fulldat, aes(x=rep, y=Taxonomy )) +
  geom_point(aes(size=Counts,color=Taxonomy)) +
  scale_colour_discrete(guide=FALSE)+
  scale_size_continuous(range=c(0,16), name="Percent") + theme_bw()+theme(legend.position = "left") +
  ggtitle("Microbiomes of Healthy and Diseased Corals from Miami")
plotorderalltime=plotorderalltime  + xlab("Pair") + ylab("Order") +
  theme(axis.text.y = element_text(color=ggplotColours(n=5)))+ facet_grid (Disease_state~., scales="free_x") 
  
plotorderalltime
dev.off()

