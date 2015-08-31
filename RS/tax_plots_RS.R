rm(list=ls())
require(ggplot2)
require(grid)

setwd("~/Documents/cbbs4_te/cbbs4_te/")

c=read.csv("cbbs4_te_allotus.csv", h=T, row.names=1)

s=read.csv("cbbs4_te_data.csv", h=T, as.is=T)

b=t(c)

d=cbind.data.frame(s,b)

check=cbind(row.names(d), d[1:1])

#remove ambiguous disease state and mid samples
d=d[d$Disease_state!="?",]

d=d[d$base_mid!="M",]

# Remove first column
d=d[2:(NCOL(d))]

c=d[8:NCOL(d)]
c=t(c)

s=d[0:7]

t=read.csv("rep_set_tax_assignments3.csv", h=T)

c1=c

rs = rowSums (c1)
use = (rs > 0)
c1 = c1 [ use, ]

#remove chloroplast
ct=merge(c1, t, by.x="row.names", by.y="OTUID")
row.names(ct)=ct$Row.names
ctnoc=ct[ct$class!="Chloroplast",]

c1=ctnoc[,2:(NCOL(ctnoc)-7)]

#convert to percentages
c1=t(c1)
c1per=(c1/rowSums(c1))*100
b1=t(c1per)

#use just order and family
tof=t[,c("OTUID", "order", "family")] 

b1t=merge(b1, tof, by.x="row.names", by.y="OTUID")

b1to=b1t[(b1t$family==""),]
b1tf=b1t[!(b1t$family==""),]
b1tf=b1tf[2:NCOL(b1tf)]
b1tf$order=NULL
b1to=b1to[2:NCOL(b1to)]
b1to$family=NULL
b1tfagg=aggregate(. ~ family, data=b1tf, FUN=sum)

##subset only those families that are found in at least 10% of one coral
b1tfagg10=b1tfagg[apply(b1tfagg[,2:NCOL(b1tfagg)], MARGIN=1, function(x) any(x>10)),]

b1toagg=aggregate(. ~ order, data=b1to, FUN=sum)
b1toagg10=b1toagg[apply(b1toagg[,2:NCOL(b1toagg)], MARGIN=1, function(x) any(x>10)),]
#b1toagg10=b1toagg10[-c(1),]

## rename columns so they are the same
colnames(b1toagg10)[1]="taxonomy"
colnames(b1tfagg10)[1]="taxonomy"

b1t10=rbind(b1tfagg10, b1toagg10)
row.names(b1t10)=b1t10$taxonomy
b1t10=b1t10[2:NCOL(b1t10)]
b1t10=t(b1t10)


fullfams2 = data.frame(Samples= numeric(0), Counts= integer(0), Taxonomy = character(0))

  for (i in 1:ncol(b1t10)) {
    fams=data.frame("Samples"=row.names(b1t10), 
    "Counts"=b1t10[,i], "Taxonomy"=colnames(b1t10)[i])
    fullfams2=rbind(fullfams2,fams)
}


fulldat=merge(fullfams2, s, by.x="Samples", by.y="row.names", all.x=TRUE)

T0=fulldat[fulldat$Day=="zero",]
T1=fulldat[fulldat$Day=="one",]
T3H=fulldat[fulldat$Day=="three"&fulldat$Disease_state=="Healthy"& fulldat$Exposure=="Healthy",]
T3DH=fulldat[fulldat$Day=="three"&fulldat$Disease_state=="Healthy"& fulldat$Exposure=="Disease",]
T3DD=fulldat[fulldat$Day=="three"&fulldat$Disease_state=="Diseased"& fulldat$Exposure=="Disease",]

ggplotColours =function(n=6, h=c(0, 360) +15){
  if ((diff(h)%%360) < 1) h[2] = h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

plot=ggplot(data=T3DD, aes(x=Samples, y=Taxonomy )) +
  geom_point(aes(size=Counts,color=Taxonomy)) +          
  scale_size_continuous(range=c(0,15)) + theme_bw()+theme(legend.position = "none")  

plot + facet_grid (.~resistant_susceptible, scales="free_x")+theme(axis.text.y = element_text(color=ggplotColours(n=5)))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
