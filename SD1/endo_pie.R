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
b1=t(c1)




##use silva taxonomy##########################
ts=read.csv("input_files/rep_set_tax_assignments_silva_clean2.csv", header=F)
colnames(ts)=c("OTUID", "kingdom", "phylum", "class", "order", "family", "genus", "species", "evalue", "ref")
tcof=ts[,c("OTUID","class", "order", "family")]

#merge tax with counts
b1t=merge(b1, tcof, by.x="row.names", by.y="OTUID")

#Divide taxonomy into class, order, and family

b1tf=b1t
b1tf$class=NULL
b1tf$order=NULL
row.names(b1tf)=b1tf$Row.names
b1tf$Row.names=NULL
b1tf=t(b1tf)

b1tfdat=merge(b1tf, s1, by.x="row.names", by.y="Sample_name", all.x=TRUE)
fulldatone=fulldat[fulldat$Timepoint=="one",]

d_genoT1_geno=ddply(d_genoT1, ~ Genotype, colwise(mean))




#aggregate by family or order
b1tfagg=aggregate(. ~ family, data=b1tf, FUN=sum)

##subset Endozoicomonas
b1tfagg10=b1tfagg[b1tfagg$family=="Hahellaceae",]


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
fulldatone=fulldat[fulldat$Timepoint=="one",]

fulldatoneagg=ddply(fulldatone, ~ Genotype, colwise(mean))

endo=cbind(fulldatoneagg$Genotype, fulldatoneagg$Counts)
colnames(endo)=c("Colony", "Endozoicomonas")
endo=as.data.frame(endo)
endo$Endozoicomonas=as.numeric(as.character(endo$Endozoicomonas))
endo$Other=(100-endo$Endozoicomonas)
row.names(endo)=endo$Colony
endo$Colony=NULL

##reshape
full = data.frame(Colony= numeric(0), Counts= integer(0), Taxonomy = character(0))

for (i in 1:ncol(endo)) {
  fams=data.frame("Colony"=row.names(endo), 
                  "Counts"=endo[,i], "Taxonomy"=colnames(endo)[i])
  full=rbind(full,fams)
}

##plot pies
bp= ggplot(full, aes(x="", y=Counts, fill=Taxonomy))+
  geom_bar(width = 1, stat = "identity")+facet_grid (.~Colony, scales="free_x") +
  scale_fill_manual(values=c("#9A8822", "#D5D5D3")) +
  theme_bw()
bp=bp+coord_polar(theta = "y")
bp
ggsave(plot=bp,height=6,width=48,dpi=1000, filename="Output_files/GLMM_results/plots/Endo_piecharts.pdf", useDingbats=FALSE)

  
