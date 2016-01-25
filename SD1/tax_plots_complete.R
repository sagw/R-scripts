rm(list=ls())
require(ggplot2)
require(grid)
setwd("~/Documents/SD1/SD1_analyses/")
c=read.csv("input_files/filtered_otutable.csv", h=T, row.names=1)
c=c[7:NCOL(c)]

#information about samples
s=read.csv("input_files//SD1_map.csv", h=T, as.is=T)
s=s[!s$Timepoint=="control",]
b=t(c)
d=cbind.data.frame(s,b)
row.names(d)=row.names(b)

d=d[!d$Sample_name=="SD1.CK14D2.CK14.O",]
d=d[!d$Sample_name=="SD1.CK14.D2.CK4.W.T57.D",]

d$Dose_disease_state=ifelse(d$Timepoint=="one", "Healthy", d$Dose_disease_state)

#check order of samples
check=cbind(row.names(d), d[1:1])

#remove past data
d[13]=NULL

experimental=d[d$Timepoint %in% c("one", "two", "three"),]

d1=experimental
c1=d1[13:NCOL(d1)]
s1=d1[0:12]

c1=t(c1)
rs = rowSums (c1)
use = (rs > 0)
c1 = c1 [ use, ]

c1=t(c1)

#convert to percent relative abundance
c1per=(c1/rowSums(c1))*100
b1=t(c1per)

##use gg taxonomy#################################
#simplify taxonomy
t=read.csv("input_files//rep_set_tax_assignments_split2.csv", header=T)
tcof=t[,c("OTUID","class", "order", "family")]

##use silva taxonomy##########################
ts=read.csv("input_files/rep_set_tax_assignments_silva_clean2.csv", header=F)
colnames(ts)=c("OTUID", "kingdom", "phylum", "class", "order", "family", "genus", "species", "evalue", "ref")
tcof=ts[,c("OTUID","class", "order", "family")]

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

##subset only those families or orders that are found in at least 15% of one coral
b1tcagg10=b1tcagg[apply(b1tcagg[,2:NCOL(b1tcagg)], MARGIN=1, function(x) any(x>15)),]
b1tfagg10=b1tfagg[apply(b1tfagg[,2:NCOL(b1tfagg)], MARGIN=1, function(x) any(x>15)),]
b1toagg10=b1toagg[apply(b1toagg[,2:NCOL(b1toagg)], MARGIN=1, function(x) any(x>15)),]

b1tcagg10$class=as.character(b1tcagg10$class)
b1tcagg10[1,1]="None"

b1toagg10$order=as.character(b1toagg10$order)
b1toagg10[1,1]="None"

b1tfagg10$family=as.character(b1tfagg10$family)
b1tfagg10[1,1]="None"

#plot class########################
row.names(b1tcagg10)=b1tcagg10$class
b1tcagg10=b1tcagg10[2:NCOL(b1tcagg10)]
b1tcagg10=t(b1tcagg10)
b1tcagg10=as.data.frame(b1tcagg10)

full = data.frame(Samples= numeric(0), Counts= integer(0), Taxonomy = character(0))

for (i in 1:ncol(b1tcagg10)) {
  fams=data.frame("Samples"=row.names(b1tcagg10), 
                  "Counts"=b1tcagg10[,i], "Taxonomy"=colnames(b1tcagg10)[i])
  full=rbind(full,fams)
}

#plot order#####################################
row.names(b1toagg10)=b1toagg10$order
b1toagg10=b1toagg10[2:NCOL(b1toagg10)]
b1toagg10=t(b1toagg10)

full = data.frame(Samples= numeric(0), Counts= integer(0), Taxonomy = character(0))

for (i in 1:ncol(b1toagg10)) {
  fams=data.frame("Samples"=row.names(b1toagg10), 
                  "Counts"=b1toagg10[,i], "Taxonomy"=colnames(b1toagg10)[i])
  full=rbind(full,fams)
}


##plot family#####################
row.names(b1tfagg10)=b1tfagg10$family
b1tfagg10=b1tfagg10[2:NCOL(b1tfagg10)]
b1tfagg10=t(b1tfagg10)

full = data.frame(Samples= numeric(0), Counts= integer(0), Taxonomy = character(0))

for (i in 1:ncol(b1tfagg10)) {
  fams=data.frame("Samples"=row.names(b1tfagg10), 
                  "Counts"=b1tfagg10[,i], "Taxonomy"=colnames(b1tfagg10)[i])
  full=rbind(full,fams)
}



##plot###############################
fulldat=merge(full, s1, by.x="Samples", by.y="Sample_name", all.x=TRUE)

CK4alltime=fulldat[fulldat$Site=="CK4",]
CK14alltime=fulldat[fulldat$Site=="CK14",]

ggplotColours =function(n=6, h=c(0, 360) +15){
  if ((diff(h)%%360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

#order levels
CK4alltime$Timepoint=factor(CK4alltime$Timepoint, levels=c("one","two","three"))
CK4alltime$Dose_disease_state=factor(CK4alltime$Dose_disease_state, levels=c("Healthy","Diseased"))
CK4alltime$T3_disease_state=factor(CK4alltime$T3_disease_state, levels=c("Healthy","Diseased"))

CK14alltime$Timepoint=factor(CK14alltime$Timepoint, levels=c("one","two","three"))
CK14alltime$Dose_disease_state=factor(CK14alltime$Dose_disease_state, levels=c("Healthy","Diseased"))
CK14alltime$T3_disease_state=factor(CK14alltime$T3_disease_state, levels=c("Healthy","Diseased"))

##plot all timepoints CK4 
####DONT FORGET TO CHANGE NAMES OF FILES AND N FOR COLOR!!!
pdf(file="Output_files/Taxonomy_plots/sFamily_CK4.pdf", width=10, height=15)
plotCK4alltime=ggplot(data=CK4alltime, aes(x=Coral, y=Taxonomy )) +
  geom_point(aes(size=Counts,color=Taxonomy)) +       
  scale_colour_discrete(guide=FALSE)+
  scale_size_continuous(range=c(0,16)) + theme_bw()+theme(legend.position = "left") +
  ggtitle("CK4 corals")
plotCK4alltime=plotCK4alltime+facet_grid (Timepoint+T3_disease_state+Dose_disease_state~Genotype, scales="free_x") +
  theme(axis.text.y = element_text(color=ggplotColours(n=16)))+ 
  theme(axis.text.x = element_blank())

plotCK4alltime
dev.off()

##plot all timepoints CK14
pdf(file="Output_files/Taxonomy_plots/sFamily_CK14.pdf", width=10, height=16)
plotCK14alltime=ggplot(data=CK14alltime, aes(x=Coral, y=Taxonomy )) +
  geom_point(aes(size=Counts,color=Taxonomy)) +
  scale_colour_discrete(guide=FALSE)+
  scale_size_continuous(range=c(0,16)) + theme_bw()+theme(legend.position = "left") +
  ggtitle("CK14 corals")
plotCK14alltime=plotCK14alltime + facet_grid (Timepoint+T3_disease_state+Dose_disease_state~Genotype, scales="free_x") +
  theme(axis.text.y = element_text(color=ggplotColours(n=16)))+ 
  theme(axis.text.x = element_blank())

plotCK14alltime
dev.off()
