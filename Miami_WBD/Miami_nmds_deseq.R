rm(list=ls())
require(vegan)
require(DESeq2)
require(genefilter)
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

# Remove first column
d=d[2:(NCOL(d))]

## exclude outlier for DESEQ...
d=d[!d$Sample_name=="D10",]
c=d[4:NCOL(d)]
c=t(c)
s=d[0:4]

dds <-DESeqDataSetFromMatrix(countData = c, colData=s, 
                             design= ~ Disease_state)
colData(dds)$Disease_state <- factor(colData(dds)$Disease_state, levels=c("Healthy","Diseased"))
colData(dds)$rep <- factor(colData(dds)$rep, levels=c("six","seven", "eight", "nine", "ten"))

#calculate different geometric means because of zeros.
gm_mean=function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0])))
geoMeans = apply(c, 1, gm_mean) 

#use locfunc=shorth instead of median because of low counts.
dds = estimateSizeFactors(dds, geoMeans=geoMeans, locfunc=shorth)

rs = rowSums ( counts ( dds))
use = (rs > 0)

dds = dds [ use, ]

comm.normalized=t(counts(dds, normalized=TRUE))

braycomm=vegdist(comm.normalized, method="bray")

clust=hclust(braycomm, method="average")
plot(clust, labels=s$rep, cex=0.7)

factors.num=3
names.factors=colnames(s)[1:(factors.num-1)]
column.num=factors.num+1

comm.mds= metaMDS(comm.normalized, autotransform=F, k=2, trymax=1000, distance="bray")
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

##plot nmds
pdf(file="output/miami_nmds_diseaseState.pdf", width=8, height=5)
plot (comm.mds$points[,1], comm.mds$points[,2], t="n", xlab="nMDS 1", xaxs="i", yaxs="i",
      ylab="nMDS 2",
      main=paste("nMDS of OTU counts  (2-D stress = ", 
                 round(comm.mds$stress, digit=2), ")", sep=""),

xlim=range(species.scores[,1]+1, comm.mds$points[,1]-1, na.rm=T),
ylim=range(species.scores[,2]+1, comm.mds$points[,2]-1, na.rm=T))
abline(h=0, lty=2)
abline(v=0, lty=2)


p=subset(nmds.coords, Disease_state=="Diseased")
points (p[,1], p[,2], pch=16, col="red", bg="red")

p=subset(nmds.coords, Disease_state=="Healthy")
points (p[,1], p[,2], pch=16, col="blue", bg="blue")

# Plot diseased centroid
loc=which(rownames(ef$factors$centroids)=="Diseased")
text (x=ef$factors$centroids[loc, "NMDS1"], 
      y=ef$factors$centroids[loc, "NMDS2"],
      rownames(ef$factors$centroids)[loc],
      col="black", cex=0.75)
# Plot healthy centroid
loc=which(rownames(ef$factors$centroids)=="Healthy")
text (x=ef$factors$centroids[loc, "NMDS1"], 
      y=ef$factors$centroids[loc, "NMDS2"],
      rownames(ef$factors$centroids)[loc],
      col="black", cex=0.75)

# Plot 95% CI ellipses based on healthy vs. diseased samples
pl <- ordiellipse(comm.mds, env$Disease_state, kind="se", conf=0.95, lwd=2, 
                  draw = "lines", 
                  col="darkgray")

#add legend
legend(x="bottomright", pch=c(16, 16), col=c("red", "blue"),
       legend=c("Diseased", "Healthy"), cex=0.8)

#add labels to points
#text(nmds.coords, labels=env$rep, cex=0.75, adj=1.5)
dev.off()

##PERMANOVA
(perm=adonis(comm.normalized ~  Disease_state, data=env, permutations=999, method="bray"))

##ENVFIT
envfit(comm.normalized, env)


#calculate shannon diversity
shannon=diversity(b, index = "shannon", MARGIN = 1, base = exp(1))
as.data.frame(shannon)
dshannon=cbind.data.frame(s,shannon)

#subset d and h
D_shannon<-subset(dshannon, Disease_state=="Diseased")
H_shannon<-subset(dshannon, Disease_state=="Healthy")

#run T tests 
D<-D_shannon$shannon
H<-H_shannon$shannon
t.test(D,H)

shannon=qplot(Disease_state, shannon, colour=Disease_state, data=dshannon)
shannon+ theme_bw()

#calculate simpson diversity
simpson<-diversity(b, index = "simpson", MARGIN = 1, base = exp(1))
dsimpson=cbind.data.frame(s,simpson)

#subset d and h
D_simpson<-subset(dsimpson, Disease_state=="Diseased")
H_simpson<-subset(dsimpson, Disease_state=="Healthy")

#run t tests
D<-D_simpson$simpson
H<-H_simpson$simpson
t.test(D,H)

simpson=qplot(rep, simpson, colour=Disease_state, data=dsimpson)
simpson+ theme_bw()

rrarefy(D)

##run rest of deseq
#filter out low aabundance counts
rs = rowSums ( counts ( dds))
use = (rs > 0)
dds = dds [ use, ]

dds=estimateDispersions(dds, fitType='local')
dds=nbinomLRT(dds, reduced=~1)

pdf(file="output/deseq_lrt_volcano.pdf", width=8, height=5)
plotMA(dds, ylim=c(-25,25), alpha=0.05)
dev.off()

res<-results(dds)
ressigs <- subset(res, padj<0.05)
ressigs<-as.data.frame(ressigs)
ressigstax<-merge(ressigs, t, by.x="row.names", by.y="OTUID", all.x=T)

write.csv(ressigstax, file="output/Miami_sigotusall_lrt.csv")

write.csv(ressigsNOD10tax, file="output/Miami_sigotusnoD10.csv")

ressigsNOD10tax<-merge(ressigs, t, by.x="row.names", by.y="OTUID", all.x=T)
