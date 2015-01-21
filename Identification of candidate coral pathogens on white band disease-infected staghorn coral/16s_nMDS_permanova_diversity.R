rm(list=ls())
require(vegan)

setwd("~/Desktop/cbbs4/cbbs4_97open_otus_coldose2/")

#this is the otu table from qiime converted to .csv
c=read.csv("complete_otu_table_frombiom.csv", h=T, row.names=1)

#information about samples
s=read.csv("dosecol_data.csv", h=T, as.is=T)
b=t(c)
d=cbind.data.frame(s,b)

#check order of samples
check<-cbind(row.names(d), d[1:1])

# Remove first column
d=d[2:(NCOL(d))]

#subset dosage and collection data
dosage<-subset(d, Experiment=="Dosage")
collection<-subset(d, Experiment=="Collection" & Sample_name!="bad")

#analyze collection  data
d=collection
c=d[6:NCOL(d)]
c=t(c)
s=d[0:5]

#filter out OTUS present only in other dataset
rs = rowSums (c)
use = (rs > 0)
c = c [ use, ]

#import normalized counts from deseq2
comm.normalized=read.csv("counts_norm_collection.csv", h=T)

factors.num=5
names.factors=colnames(d)[1:(factors.num-1)]
column.num=factors.num+1

# Run nMDS: because nMDS is partially stochastic, metaMDS actually performs a 
# number of nMDS based on different initial conditions to determine whether the
# final ordination is identical (or at the very least very similar) across the
# different nMDS. If they are, then the algorithm is said to have converged.
comm.mds <- metaMDS(comm.normalized, autotransform=F, k=2, trymax=1000)

# Compute Rsq
(Rsq=1-(comm.mds$stress)^2)

# Fit explanatory factors onto ordination space
env=d[, 1:(factors.num-1)]

# Remove spaces in label names
env=sapply(env, FUN=function (x) gsub(x, pattern=" ", replace=""))

# Convert table to data frame and characters to factors
env=as.data.frame(env)
ef=envfit(comm.mds, env)

# Shorten labels by removing factor name and only using level names
for (i in 1:(factors.num-1)) {
  labs=paste0(colnames(d)[i], levels(env[, i]))
  rownames(ef$factors$centroids)[which(rownames(ef$factors$centroids) %in% labs)]=levels(env[, i])  
}

#plot nmds
# Prepare an empty plot
pdf(file="nmds_collection.pdf", width=8, height=8)
par(tck=0.02, mgp=c(1.5, 0.2, 0))
plot (comm.mds$points[,1], comm.mds$points[,2], t="n", xlab="nMDS 1", xaxs="i", yaxs="i",
      ylab="nMDS 2",
      main=paste("nMDS of OTU counts \n field-collected (2-D stress = ", 
                 round(comm.mds$stress, digit=2), ")", sep=""),
      xlim=range(species.scores[,1]+1, comm.mds$points[,1]-1, na.rm=T),
      ylim=range(species.scores[,2], comm.mds$points[,2]-1, na.rm=T))

## Draw grid with dashed lines about origin
abline(h=0, lty=2)
abline(v=0, lty=2)

legend(x="bottomright", pch=c(1, 2, 0, 5, 1, 2, 0, 5), 
       col=c("blue", "blue", "blue", "blue", "red", "red", "red", "red"),
       legend=c("CK4, Healthy", "CK5, Healthy", "CK6, Healthy", "POPA, Healthy", "CK4, 
                Diseased", "CK5", "Diseased", "CK6, Diseased", "POPA", "Diseased"))

#Use different symbols for disease and site
p=subset(nmds.coords, Disease_state=="Healthy" & Site=="CK4")
points (p[,1], p[,2], pch=1, col="blue", bg="blue")

p=subset(nmds.coords, Disease_state=="Healthy"  & Site=="CK5")
points (p[,1], p[,2], pch=2, col="blue", bg="blue")

p=subset(nmds.coords, Disease_state=="Healthy" & Site=="CK6")
points (p[,1], p[,2], pch=0, col="blue", bg="blue")

p=subset(nmds.coords, Disease_state=="Healthy" & Site=="POPA")
points (p[,1], p[,2], pch=5, col="blue", bg="blue")

p=subset(nmds.coords, Disease_state=="Diseased" & Site=="CK4")
points (p[,1], p[,2], pch=1, col="red", bg="red")

p=subset(nmds.coords, Disease_state=="Diseased"  & Site=="CK5")
points (p[,1], p[,2], pch=2, col="red", bg="red")

p=subset(nmds.coords, Disease_state=="Diseased" & Site=="CK6")
points (p[,1], p[,2], pch=0, col="red", bg="red")

p=subset(nmds.coords, Disease_state=="Diseased"  & Site=="POPA")
points (p[,1], p[,2], pch=5, col="red", bg="red")

# Plot diseased centroid
loc=which(rownames(ef$factors$centroids)=="Diseased")
text (x=ef$factors$centroids[loc, "NMDS1"], 
      y=ef$factors$centroids[loc, "NMDS2"],
      rownames(ef$factors$centroids)[loc],
      col="black")
# Plot healthy centroid
loc=which(rownames(ef$factors$centroids)=="Healthy")
text (x=ef$factors$centroids[loc, "NMDS1"], 
      y=ef$factors$centroids[loc, "NMDS2"],
      rownames(ef$factors$centroids)[loc],
      col="black")

# Plot CK4 centroid
loc=which(rownames(ef$factors$centroids)=="CK4")
text (x=ef$factors$centroids[loc, "NMDS1"], 
      y=ef$factors$centroids[loc, "NMDS2"],
      rownames(ef$factors$centroids)[loc],
      col="darkgray")
# Plot CK5 centroid
loc=which(rownames(ef$factors$centroids)=="CK5")
text (x=ef$factors$centroids[loc, "NMDS1"], 
      y=ef$factors$centroids[loc, "NMDS2"],
      rownames(ef$factors$centroids)[loc],
      col="darkgray")

# Plot CK6 centroid
loc=which(rownames(ef$factors$centroids)=="CK6")
text (x=ef$factors$centroids[loc, "NMDS1"], 
      y=ef$factors$centroids[loc, "NMDS2"],
      rownames(ef$factors$centroids)[loc],
      col="darkgray")

# Plot POPA centroid
loc=which(rownames(ef$factors$centroids)=="POPA")
text (x=ef$factors$centroids[loc, "NMDS1"], 
      y=ef$factors$centroids[loc, "NMDS2"],
      rownames(ef$factors$centroids)[loc],
      col="darkgray")

# Plot 95% CI ellipses based on healthy vs. diseased samples
pl <- ordiellipse(comm.mds, env$Disease_state, kind="se", conf=0.95, lwd=2, 
                  draw = "lines", 
                  col="black")
dev.off()

# run permanova
(perm=adonis(comm.normalized ~ Disease_state*Site*Year, data=env, permutations=999, method="bray"))
write.csv(file=paste0("permanova_dosedisease_proposal1.csv"), perm$aov.tab)

#Look at dispersion
#distance matrix
dm<-vegdist(comm.normalized, method="bray", binary=FALSE, diag=FALSE, upper=FALSE,
            na.rm = FALSE)
#set factors
site.factors <- env$Site
disease.factors <- env$Disease_state
year.factors <- env$Year
dose.factors <- env$Dose
interaction<-interaction(disease.factors, site.factors, year.factors)

#create models
mod<-betadisper(dm, interaction, bias.adjust = FALSE)
mod1<-betadisper(dm, disease.factors, bias.adjust = FALSE)
mod2<-betadisper(dm, year.factors, bias.adjust = FALSE)

#plot and run anovas
plot(mod)
anova(mod)
permutest(mod)
mod.hsd<-TukeyHSD(mod)
plot(mod.hsd)
boxplot(mod.hsd)

plot(mod1)
anova(mod1)
permutest(mod1)
mod.hsd<-TukeyHSD(mod1)
plot(mod.hsd1)
boxplot(mod.hsd1)

plot(mod2)
anova(mod2)
permutest(mod2)
mod.hsd<-TukeyHSD(mod2)
plot(mod.hsd)
boxplot(mod.hsd)

#Diversity
bs=t(c)

#calculate shannon diversity
shannon<-diversity(bs, index = "shannon", MARGIN = 1, base = exp(1))
as.data.frame(shannon)
dshannon=cbind.data.frame(s,shannon)

#subset d and h
ColD_shannon<-subset(dshannon, Disease_state=="Disease")
ColH_shannon<-subset(dshannon, Disease_state=="Healthy")

#run T test 
D<-ColD_shannon$shannon
H<-ColH_shannon$shannon
t.test(D,H)

#calculate simpson diversity
simpson<-diversity(bs, index = "simpson", MARGIN = 1, base = exp(1))
as.data.frame(simpson)
dsimpson=cbind.data.frame(s,simpson)

#subset d and h
ColD_simpson<-subset(dsimpson, Disease_state=="Diseased")
ColH_simpson<-subset(dsimpson, Disease_state=="Healthy")

#run t test
D<-ColD_simpson$simpson
H<-ColH_simpson$simpson
t.test(D,H)

#rarefied richness
rich<-rarefy(bs, sample=589,857)
rich<-as.data.frame(rich)
rich=t(rich)
drich=cbind.data.frame(s,rich)
check<-cbind(row.names(drich), drich[1:1])
D_rich<-subset(drich, Disease_state=="Diseased")
H_rich<-subset(drich, Disease_state=="Healthy")

#run t test
D<-D_rich$S
H<-H_rich$S
t.test(D,H)

#analyze dose  data
d=dosage
c=d[6:NCOL(d)]
c=t(c)
s=d[0:5]

#filter out OTUS present only in other dataset
rs = rowSums (c)
use = (rs > 0)
c = c [ use, ]

#import normalized counts from deseq2
comm.normalized=read.csv("counts_norm_dose.csv", h=T)

factors.num=5
names.factors=colnames(d)[1:(factors.num-1)]
column.num=factors.num+1

# Run nMDS: because nMDS is partially stochastic, metaMDS actually performs a 
# number of nMDS based on different initial conditions to determine whether the
# final ordination is identical (or at the very least very similar) across the
# different nMDS. If they are, then the algorithm is said to have converged.
comm.mds <- metaMDS(comm.normalized, autotransform=F, k=2, trymax=1000)
# Compute Rsq
(Rsq=1-(comm.mds$stress)^2)
# Flip axes
comm.mds$points=comm.mds$points*-1
nmds.coords <- cbind(comm.mds$points, d[, 1:(factors.num-1)])
species.scores=comm.mds$species

# Fit explanatory factors onto ordination space
env=d[, 1:(factors.num-1)]
# Remove spaces in label names
env=sapply(env, FUN=function (x) gsub(x, pattern=" ", replace=""))
# Convert table to data frame and characters to factors
env=as.data.frame(env)

ef=envfit(comm.mds, env)
# Shorten labels by removing factor name and only using level names
for (i in 1:(factors.num-1)) {
  labs=paste0(colnames(d)[i], levels(env[, i]))
  rownames(ef$factors$centroids)[which(rownames(ef$factors$centroids) %in% labs)]=levels(env[, i])  
}

# Prepare an empty plot
pdf(file="nmds_dose.pdf", width=8, height=8)
par(tck=0.02, mgp=c(1.5, 0.2, 0))
plot (comm.mds$points[,1], comm.mds$points[,2], t="n", xlab="nMDS 1", xaxs="i", yaxs="i",
      ylab="nMDS 2",
      main=paste("nMDS of OTU counts \n tank-exposed (2-D stress = ", 
                 round(comm.mds$stress, digit=2), ")", sep=""),
      xlim=range(species.scores[,1]+1, comm.mds$points[,1]-1, na.rm=T),
      ylim=range(species.scores[,2], comm.mds$points[,2]-1, na.rm=T))
## Draw grid with dashed lines about origin
abline(h=0, lty=2)
abline(v=0, lty=2)

legend(x="bottomright", pch=c(1, 16, 1), col=c("blue", "blue", "red"),
       legend=c("Dose=Disease, State=Healthy", "Dose=Healthy, State=Healthy", "Dose=Disease, State=Diseased"))

#plot samples by dose and disease state
p=subset(nmds.coords, Dose=="D" & Disease_state=="Healthy")
points (p[,1], p[,2], pch=1, col="blue", bg="blue")

p=subset(nmds.coords, Dose=="H" & Disease_state=="Healthy")
points (p[,1], p[,2], pch=16, col="blue", bg="blue")

p=subset(nmds.coords, Dose=="D" & Disease_state=="Diseased")
points (p[,1], p[,2], pch=1, col="red", bg="red")  

# Plot diseased centroid
loc=which(rownames(ef$factors$centroids)=="Diseased")
text (x=ef$factors$centroids[loc, "NMDS1"], 
      y=ef$factors$centroids[loc, "NMDS2"],
      rownames(ef$factors$centroids)[loc],
      col="red")
# Plot healthy centroid
loc=which(rownames(ef$factors$centroids)=="Healthy")
text (x=ef$factors$centroids[loc, "NMDS1"], 
      y=ef$factors$centroids[loc, "NMDS2"],
      rownames(ef$factors$centroids)[loc],
      col="blue")

# Plot 95% CI ellipses based on healthy vs. diseased samples
pl <- ordiellipse(comm.mds, env$Disease_state, kind="se", conf=0.95, lwd=2, 
                  draw = "lines", 
                  col="black")

dev.off()

# run permanova
(perm=adonis(comm.normalized ~ Disease_state*Site*Year, data=env, permutations=999, method="bray"))
write.csv(file=paste0("permanova_dosedisease_proposal1.csv"), perm$aov.tab)

#Look at dispersion
#distance matrix
dm<-vegdist(comm.normalized, method="bray", binary=FALSE, diag=FALSE, upper=FALSE,
            na.rm = FALSE)
#set factors
disease.factors <- env$Disease_state

#create model
mod1<-betadisper(dm, disease.factors, bias.adjust = FALSE)

#plot and run anovas
plot(mod)
anova(mod)
permutest(mod)
mod.hsd<-TukeyHSD(mod)
plot(mod.hsd)
boxplot(mod.hsd)

#Diversity
bs=t(c)

#calculate shannon diversity
shannon<-diversity(bs, index = "shannon", MARGIN = 1, base = exp(1))
as.data.frame(shannon)
dshannon=cbind.data.frame(s,shannon)

#subset d and h
DoseDD_simpson<-subset(dsimpson, Disease_state=="Diseased")
DoseHH_simpson<-subset(dsimpson, Disease_state=="Healthy" & Dose=="H")
DoseDH_simpson<-subset(dsimpson, Disease_state=="Healthy" & Dose=="D")

#run T tests 
D<-DoseDD_shannon$shannon
H<-DoseHH_shannon$shannon
t.test(D,H)

D<-DoseDD_shannon$shannon
H<-DoseDH_shannon$shannon
t.test(D,H)

D<-DoseDH_shannon$shannon
H<-DoseHH_shannon$shannon
t.test(D,H)

#calculate simpson diversity
simpson<-diversity(bs, index = "simpson", MARGIN = 1, base = exp(1))
as.data.frame(simpson)
dsimpson=cbind.data.frame(s,simpson)

#subset d and h
DoseDD_simpson<-subset(dsimpson, Disease_state=="Diseased")
DoseHH_simpson<-subset(dsimpson, Disease_state=="Healthy" & Dose=="H")
DoseDH_simpson<-subset(dsimpson, Disease_state=="Diseased" & Dose=="H")

#run t tests
D<-DoseDD_simpson$simpson
H<-DoseHH_simpson$simpson
t.test(D,H)

D<-DoseDD_simpson$simpson
H<-DoseDH_simpson$simpson
t.test(D,H)

D<-DoseDH_simpson$simpson
H<-DoseHH_simpson$simpson
t.test(D,H)

#rarefied richness
rich<-rarefy(bs, sample=589,857)
rich<-as.data.frame(rich)
rich=t(rich)
drich=cbind.data.frame(s,rich)
check<-cbind(row.names(drich), drich[1:1])
DD_rich<-subset(drich, Disease_state=="Diseased")
HH_rich<-subset(drich, Disease_state=="Healthy", & Dose=="H")
DH_rich<-subset(drich, Disease_state=="Diseased", & Dose=="H")

#run t tests
D<-DD_rich$S
H<-HH_rich$S
t.test(D,H)

D<-DD_rich$S
H<-DH_rich$S
t.test(D,H)

D<-DH_rich$S
H<-HH_rich$S
t.test(D,H)

