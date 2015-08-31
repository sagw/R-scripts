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
#check order of samples
check=cbind(row.names(d), d[1:1])

t=read.csv("input_files//rep_set_tax_assignments_split2.csv", header=T)

d[14]=NULL
d[13]=NULL

experimental=d[d$Timepoint %in% c("one", "two", "three", "dose"),]

timeone=d[d$Timepoint=="one",]
timetwo=d[d$Timepoint=="two",]
timethree=d[d$Timepoint=="three",]


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
#b1=t(c1)
##b125=b1[apply(b1, MARGIN=1, function(x) (sum(x>0)>(0.25*length(x)))), ]

#simplify taxonomy
tof=t[,c("OTUID", "order", "family")]

#merge tax with counts
b1t=merge(b1, tof, by.x="row.names", by.y="OTUID")

#subdivide into otus with taxonomy for family and order
b1to=b1t[(b1t$family==""),]
b1to=b1to[2:NCOL(b1to)]
b1to$family=NULL

b1tf=b1t[!(b1t$family==""),]
b1tf=b1tf[2:NCOL(b1tf)]
b1tf$order=NULL

#aggregate by family or order
b1tfagg=aggregate(. ~ family, data=b1tf, FUN=sum)
b1toagg=aggregate(. ~ order, data=b1to, FUN=sum)

##subset only those families or orders that are found in at least 20% of one coral
b1tfagg10=b1tfagg[apply(b1tfagg[,2:NCOL(b1tfagg)], MARGIN=1, function(x) any(x>20)),]
b1toagg10=b1toagg[apply(b1toagg[,2:NCOL(b1toagg)], MARGIN=1, function(x) any(x>20)),]

##subset samples with greater than 50% Endozoicomonas
b1tendo=b1tfagg[b1tfagg$family=="Endozoicimonaceae",]
b1tendo=t(b1tendo)
b1tendo=as.data.frame(b1tendo)
b1tendo$nada="20"
b1tendo=b1tendo[-1,]
colnames(b1tendo)=c("count", "Endo")
b1tendo1=apply(b1tendo, 2, as.numeric)
row.names(b1tendo1)=row.names(b1tendo)
endo20=b1tendo[b1tendo1[,1]>20,]

row.names(endo50)=gsub("4D", "4.D", row.names(endo50))
row.names(endo50)=gsub("4H", "4.H", row.names(endo50))

row.names(endo20)=gsub("4D", "4.D", row.names(endo20))
row.names(endo20)=gsub("4H", "4.H", row.names(endo20))

##remove no taxonomy b1toagg10=b1toagg10[-c(1),]

## rename columns so they are the same
colnames(b1toagg10)[1]="taxonomy"
colnames(b1tfagg10)[1]="taxonomy"

b1t10=rbind(b1tfagg10, b1toagg10)
row.names(b1t10)=b1t10$taxonomy
b1t10=b1t10[2:NCOL(b1t10)]
b1t10=t(b1t10)

full = data.frame(Samples= numeric(0), Counts= integer(0), Taxonomy = character(0))

for (i in 1:ncol(b1t10)) {
  fams=data.frame("Samples"=row.names(b1t10), 
                  "Counts"=b1t10[,i], "Taxonomy"=colnames(b1t10)[i])
  full=rbind(full,fams)
}

fulldat=merge(full, s1, by.x="Samples", by.y="Sample_name", all.x=TRUE)


T1CK4=fulldat[fulldat$Site=="CK4"& fulldat$Timepoint=="one",]

T2CK4H=fulldat[fulldat$Site=="CK4"& fulldat$Timepoint=="two"& 
                 fulldat$Dose_disease_state=="Healthy" & fulldat$T3_disease_state=="Healthy",]
T2CK4D=fulldat[fulldat$Site=="CK4"& fulldat$Timepoint=="two"& 
                 fulldat$Dose_disease_state=="Diseased" & fulldat$T3_disease_state=="Diseased",]
T2CK4DH=fulldat[fulldat$Site=="CK4"& fulldat$Timepoint=="two"& 
                 fulldat$Dose_disease_state=="Diseased" & fulldat$T3_disease_state=="Healthy",]
T3CK4H=fulldat[fulldat$Site=="CK4"& fulldat$Timepoint=="three"& 
                 fulldat$Dose_disease_state=="Healthy" & fulldat$T3_disease_state=="Healthy",]
T3CK4D=fulldat[fulldat$Site=="CK4"& fulldat$Timepoint=="three"& 
                 fulldat$Dose_disease_state=="Diseased" & fulldat$T3_disease_state=="Diseased",]
T3CK4DH=fulldat[fulldat$Site=="CK4"& fulldat$Timepoint=="three"& 
                  fulldat$Dose_disease_state=="Diseased" & fulldat$T3_disease_state=="Healthy",]

T1CK14=fulldat[fulldat$Site=="CK14"& fulldat$Timepoint=="one",]

T2CK14H=fulldat[fulldat$Site=="CK14"& fulldat$Timepoint=="two"& 
                 fulldat$Dose_disease_state=="Healthy" & fulldat$T3_disease_state=="Healthy",]
T2CK14HD=fulldat[fulldat$Site=="CK14"& fulldat$Timepoint=="two"& 
                  fulldat$Dose_disease_state=="Healthy" & fulldat$T3_disease_state=="Diseased",]
T2CK14DH=fulldat[fulldat$Site=="CK14"& fulldat$Timepoint=="two"& 
                  fulldat$Dose_disease_state=="Diseased" & fulldat$T3_disease_state=="Healthy",]
T2CK14D=fulldat[fulldat$Site=="CK14"& fulldat$Timepoint=="two"& 
                   fulldat$Dose_disease_state=="Diseased" & fulldat$T3_disease_state=="Diseased",]
T3CK14H=fulldat[fulldat$Site=="CK14"& fulldat$Timepoint=="three"& 
                   fulldat$Dose_disease_state=="Healthy" & fulldat$T3_disease_state=="Healthy",]
T3CK14HD=fulldat[fulldat$Site=="CK14"& fulldat$Timepoint=="three"& 
                   fulldat$Dose_disease_state=="Healthy" & fulldat$T3_disease_state=="Diseased",]
T3CK14DH=fulldat[fulldat$Site=="CK14"& fulldat$Timepoint=="three"& 
                   fulldat$Dose_disease_state=="Diseased" & fulldat$T3_disease_state=="Healthy",]
T3CK14D=fulldat[fulldat$Site=="CK14"& fulldat$Timepoint=="three"& 
                  fulldat$Dose_disease_state=="Diseased" & fulldat$T3_disease_state=="Diseased",]


ggplotColours =function(n=6, h=c(0, 360) +15){
  if ((diff(h)%%360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
#plot=plot+guides(colour = guide_legend(override.aes = list(size=5)))

plotCK4=ggplot(data=T1CK4, aes(x=Coral, y=Taxonomy )) +
  geom_point(aes(size=Counts,color=Taxonomy)) +          
  scale_size_continuous(range=c(0,16)) + theme_bw()+theme(legend.position = "none") +
  ggtitle("CK4 Time 1")
plotCK4=plotCK4 + facet_grid (.~Genotype, scales="free_x") +
  theme(axis.text.y = element_text(color=ggplotColours(n=12)))+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

plotCK4

plot2CK4H=ggplot(data=T2CK4H, aes(x=Coral, y=Taxonomy )) +
  geom_point(aes(size=Counts,color=Taxonomy)) +          
  scale_size_continuous(range=c(0,8)) + theme_bw()+theme(legend.position = "none")  
plot2CK4H=plot2CK4H + facet_grid (.~Genotype, scales="free_x") +
  theme(axis.text.y = element_text(color=ggplotColours(n=12)))+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("CK4 Time 2 Healthy")

plot2CK4H


plot2CK4DH=ggplot(data=T2CK4DH, aes(x=Coral, y=Taxonomy )) +
  geom_point(aes(size=Counts,color=Taxonomy)) +          
  scale_size_continuous(range=c(0,9)) + theme_bw()+theme(legend.position = "none")  
plot2CK4DH=plot2CK4DH + facet_grid (.~Genotype, scales="free_x") +
  theme(axis.text.y = element_text(color=ggplotColours(n=12)))+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("CK4 Time 2 Disease-exposed Healthy")

plot2CK4DH


plot2CK4D=ggplot(data=T2CK4D, aes(x=Coral, y=Taxonomy )) +
  geom_point(aes(size=Counts,color=Taxonomy)) +          
  scale_size_continuous(range=c(0,10)) + theme_bw()+theme(legend.position = "none")  
plot2CK4D=plot2CK4D + facet_grid (.~Genotype, scales="free_x") +
  theme(axis.text.y = element_text(color=ggplotColours(n=12)))+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("CK4 Time 2 Diseased")

plot2CK4D

plot3CK4H=ggplot(data=T3CK4H, aes(x=Coral, y=Taxonomy )) +
  geom_point(aes(size=Counts,color=Taxonomy)) +          
  scale_size_continuous(range=c(0,10)) + theme_bw()+theme(legend.position = "none")  
plot3CK4H=plot3CK4H + facet_grid (.~Genotype, scales="free_x") +
  theme(axis.text.y = element_text(color=ggplotColours(n=12)))+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("CK4 Time 3 Healthy") 

plot3CK4H


plot3CK4DH=ggplot(data=T3CK4DH, aes(x=Coral, y=Taxonomy )) +
  geom_point(aes(size=Counts,color=Taxonomy)) +          
  scale_size_continuous(range=c(0,5.5)) + theme_bw()+theme(legend.position = "none")  
plot3CK4DH=plot3CK4DH + facet_grid (.~Genotype, scales="free_x") +
  theme(axis.text.y = element_text(color=ggplotColours(n=12)))+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("CK4 Time 3 Disease-exposed Healthy")

plot3CK4DH

plot3CK4D=ggplot(data=T3CK4D, aes(x=Coral, y=Taxonomy )) +
  geom_point(aes(size=Counts,color=Taxonomy)) +        
  scale_size_continuous(range=c(0,12)) + theme_bw()+theme(legend.position = "none")  
plot3CK4D=plot3CK4D + facet_grid (.~Genotype, scales="free_x") +
  theme(axis.text.y = element_text(color=ggplotColours(n=12)))+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("CK4 Time 3 Diseased")

plot3CK4D

plotCK14=ggplot(data=T1CK14, aes(x=Coral, y=Taxonomy )) +
  geom_point(aes(size=Counts,color=Taxonomy)) +          
  scale_size_continuous(range=c(0,16)) + theme_bw()+theme(legend.position = "left")  
plotCK14=plotCK14 + facet_grid (.~Genotype, scales="free_x") +
  theme(axis.text.y = element_text(color=ggplotColours(n=12)))+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("CK14 Time 1")

plotCK14

plot2CK14H=ggplot(data=T2CK14H, aes(x=Coral, y=Taxonomy )) +
  geom_point(aes(size=Counts,color=Taxonomy)) +          
  scale_size_continuous(range=c(0,10)) + theme_bw()+theme(legend.position = "none")  
plot2CK14H=plot2CK14H + facet_grid (.~Genotype, scales="free_x") +
  theme(axis.text.y = element_text(color=ggplotColours(n=12)))+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("CK14 Time 2 Healthy")

plot2CK14H

plot2CK14DH=ggplot(data=T2CK14DH, aes(x=Coral, y=Taxonomy )) +
  geom_point(aes(size=Counts,color=Taxonomy)) +          
  scale_size_continuous(range=c(0,10)) + theme_bw()+theme(legend.position = "none")  
plot2CK14DH=plot2CK14DH + facet_grid (.~Genotype, scales="free_x") +
  theme(axis.text.y = element_text(color=ggplotColours(n=12)))+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("CK14 Time 2 Disease-exposed Healthy")

plot2CK14DH

plot2CK14D=ggplot(data=T2CK14D, aes(x=Coral, y=Taxonomy )) +
  geom_point(aes(size=Counts,color=Taxonomy)) +          
  scale_size_continuous(range=c(0,7)) + theme_bw()+theme(legend.position = "none")  
plot2CK14D=plot2CK14D + facet_grid (.~Genotype, scales="free_x") +
  theme(axis.text.y = element_text(color=ggplotColours(n=12)))+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("CK14 Time 2 Diseased")

plot2CK14D

plot3CK14H=ggplot(data=T3CK14H, aes(x=Coral, y=Taxonomy )) +
  geom_point(aes(size=Counts,color=Taxonomy)) +          
  scale_size_continuous(range=c(0,15)) + theme_bw()+theme(legend.position = "none")  
plot3CK14H=plot3CK14H + facet_grid (.~Genotype, scales="free_x") +
  theme(axis.text.y = element_text(color=ggplotColours(n=13)))+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("CK14 Time 3 Healthy")
plot3CK14H

plot3CK14DH=ggplot(data=T3CK14DH, aes(x=Coral, y=Taxonomy )) +
  geom_point(aes(size=Counts,color=Taxonomy)) +          
  scale_size_continuous(range=c(0,9)) + theme_bw()+theme(legend.position = "none")  
plot3CK14DH=plot3CK14DH + facet_grid (.~Genotype, scales="free_x") +
  theme(axis.text.y = element_text(color=ggplotColours(n=12)))+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("CK14 Time 3 Disease-exposed Healthy")
plot3CK14DH  

plot3CK14D=ggplot(data=T3CK14D, aes(x=Coral, y=Taxonomy )) +
  geom_point(aes(size=Counts,color=Taxonomy)) +          
  scale_size_continuous(range=c(0,13)) + theme_bw()+theme(legend.position = "none")  
plot3CK14D=plot3CK14D + facet_grid (.~Genotype, scales="free_x") +
  theme(axis.text.y = element_text(color=ggplotColours(n=12)))+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("CK14 Time 3 Diseased")
plot3CK14D

multiplot(plotCK4, plot3CK14D)

multiplot(plot2CK4D,plot2CK4H, plot2CK4DH)






# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


