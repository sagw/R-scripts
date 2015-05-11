rm(list=ls())
require(vegan)
require(DESeq2)
require(phyloseq)
require(genefilter)
require(ape)
require(ggplot2)
setwd("~/Documents/SD1/SD1_analyses/")
c=read.csv("input_files//filtered_otutable.csv", h=T, row.names=1)
c=c[5:NCOL(c)]

#information about samples
s=read.csv("input_files//SD1_map.csv", h=T, as.is=T)
b=t(c)
d=cbind.data.frame(s,b)
check=cbind(row.names(d), d[1:1])
experimental=d[d$Timepoint %in% c("one", "two", "three", "dose"),]
d=experimental
c=d[11:NCOL(d)]
c=t(c)
rs = rowSums (c)
use = (rs > 0)
c = c [ use, ]

s=d[0:10]
row.names(s)=s$Sample.name
s=s[2:NCOL(s)]


#Use DESEQ2 to normalize
dds =DESeqDataSetFromMatrix(countData = c, colData=s, 
                             design= ~ final_disease_state)

colData(dds)$Disease_state = factor(colData(dds)$final_disease_state, levels=c("Healthy","Diseased"))

#calculate different geometric means because of zeros.
gm_mean=function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0])))
geoMeans = apply(c, 1, gm_mean) 

#use locfunc=shorth instead of median because of low counts.
dds = estimateSizeFactors(dds, geoMeans=geoMeans, locfunc=shorth)

comm.normalized=(counts(dds, normalized=TRUE))
row.names(comm.normalized)=row.names(c)
otus=(comm.normalized)
otus=as.data.frame(otus)


t=read.csv("input_files//rep_set_tax_assignments_split2.csv", header=T, row.names=1)
t=as.matrix(t)

tree=read.tree("input_files/rep_set_filtered3_phylogeny.tre")

#make into phyloseq object
TREE= tree
OTU = otu_table(comm.normalized, taxa_are_rows = TRUE)
TAX = tax_table(t)
DAT= sample_data(s)
physeq = phyloseq(OTU, TAX, DAT, TREE)


#diversity metrics
#use raw counts
#make into phyloseq object
TREE= tree
OTU = otu_table(c, taxa_are_rows = TRUE)
TAX = tax_table(t)
DAT= sample_data(s)
physeq = phyloseq(OTU, TAX, DAT, TREE)

plot_richness(physeq, x="Timepoint", measures=c("Chao1", "Shannon", "Simpson"), color="final_disease_state")
physeq$Timepoint2 = factor(physeq$Timepoint, c("dose", "one", "two", "three"))
sample_data(physeq)$Timepoint = factor(sample_data(physeq)$Timepoint, levels = c("dose", "one", "two", "three"))

filt = filter_taxa(physeq, function(x) sum(x > 0) > (0.1*length(x)), TRUE)

#distance methods
#plot all distance methods
dist_methods <- unlist(distance("list"))
print(dist_methods)

# Remove them from the vector
dist_methods <- dist_methods[-(1:2)]
dist_methods = dist_methods[-which(dist_methods == "ANY")]
dist_methods = dist_methods[-which(dist_methods == "manhattan")]
plist <- vector("list", length(dist_methods))
names(plist) = dist_methods

for (i in dist_methods) {
  # Calculate distance matrix
  iDist <- distance(filt, method = i)
  # Calculate ordination
  iMDS <- ordinate(filt, "MDS", distance = iDist)
  ## Make plot Don't carry over previous plot (if error, p will be blank)
  p <- NULL
  # Create plot, store as temp variable, p
  p <- plot_ordination(filt, iMDS, color = "T3_disease_state", shape = "Timepoint")
  # Add title to each plot
  p <- p + ggtitle(paste("MDS using distance method ", i, sep = ""))
  # Save the graphic to file.
  plist[[i]] = p
}


df = ldply(plist, function(x) x$data)
names(df)[1] <- "distance"
p = ggplot(df, aes(Axis.1, Axis.2, color = T3_disease_state, shape = Timepoint))
p = p + geom_point(size = 3, alpha = 0.5)
p = p + facet_wrap(~distance, scales = "free")
p = p + ggtitle("MDS on various distance metrics for SD1 dataset")
p

print(plist[["jsd"]])
print(plist[["bray"]])
print(plist[["jaccard"]])
print(plist[["binomial"]])
print(plist[["-1"]])
print(plist[["wb"]])
print(plist[["gl"]])


#bar_plots
filtt1 = subset_samples(filt, Timepoint=="one")


filtt2 = subset_samples(filt, Timepoint=="two")
filtt2ck4d = subset_samples(filtt2, Dose_disease_state=="Diseased" & Dose.Site=="CK4")

filtt2ck14d = subset_samples(filtt2, Dose_disease_state=="Diseased" & Dose.Site=="CK14")
filtt2ck14h = subset_samples(filtt2, Dose_disease_state=="Healthy" & Dose.Site=="CK14")

filtt3 = subset_samples(filt, Timepoint=="three")
filtt3ck14d = subset_samples(filtt3, final_disease_state=="Diseased" & Dose.Site=="CK14")
filtt3ck4d = subset_samples(filtt3, final_disease_state=="Diseased" & Dose.Site=="CK4")
filtt3ck4h = subset_samples(filtt3, final_disease_state=="Healthy" & Dose.Site=="CK4")
filtt3ck14h = subset_samples(filtt3, final_disease_state=="Healthy" & Dose.Site=="CK14")

filtproteo = subset_taxa(filt, phylum=="Proteobacteria")

filtt1p = subset_samples(filtproteo, Timepoint=="one")
filtt1ck14p = subset_samples(filtt1p, Dose.Site=="CK14")
filtt1ck4p = subset_samples(filtt1p, Dose.Site=="CK4")

filtt2p = subset_samples(filtproteo, Timepoint=="two")
filtt2ck4dp = subset_samples(filtt2p, Dose_disease_state=="Diseased" & Dose.Site=="CK4")
filtt2ck4hp = subset_samples(filtt2p, Dose_disease_state=="Healthy" & Dose.Site=="CK4")
filtt2ck14dp = subset_samples(filtt2p, Dose_disease_state=="Diseased" & Dose.Site=="CK14")
filtt2ck14hp = subset_samples(filtt2p, Dose_disease_state=="Healthy" & Dose.Site=="CK14")

filtt3p = subset_samples(filtproteo, Timepoint=="three")
filtt3ck14dp = subset_samples(filtt3p, final_disease_state=="Diseased" & Dose.Site=="CK14")
filtt3ck4dp = subset_samples(filtt3p, final_disease_state=="Diseased" & Dose.Site=="CK4")
filtt3ck4hp = subset_samples(filtt3p, final_disease_state=="Healthy" & Dose_disease_state=="Healthy" & Dose.Site=="CK4")
filtt3ck14hp = subset_samples(filtt3p, final_disease_state=="Healthy" & Dose_disease_state=="Healthy" & Dose.Site=="CK14")

filtgamma = subset_taxa(filt, class=="Gammaproteobacteria")
filtt1g = subset_samples(filtgamma, Timepoint=="one")
filtt1ck14g = subset_samples(filtt1g, Dose.Site=="CK14")
filtt1ck4g = subset_samples(filtt1g, Dose.Site=="CK4")

filtt2g = subset_samples(filtgamma, Timepoint=="two")
filtt2ck4dg = subset_samples(filtt2g, T3_disease_state=="Diseased" & Dose.Site=="CK4")
filtt2ck4hg = subset_samples(filtt2g, Dose_disease_state=="Healthy" & Dose.Site=="CK4")
filtt2ck14dg = subset_samples(filtt2g, Dose_disease_state=="Diseased" & Dose.Site=="CK14")
filtt2ck14hg = subset_samples(filtt2g, Dose_disease_state=="Healthy" & Dose.Site=="CK14")

filtt3g = subset_samples(filtgamma, Timepoint=="three")
filtt3ck14dg = subset_samples(filtt3g, final_disease_state=="Diseased" & Dose.Site=="CK14")
filtt3ck4dg = subset_samples(filtt3g, final_disease_state=="Diseased" & Dose.Site=="CK4")
filtt3ck4hg = subset_samples(filtt3g, final_disease_state=="Healthy"  & Dose.Site=="CK4")
filtt3ck14hg = subset_samples(filtt3g, final_disease_state=="Healthy" & Dose.Site=="CK14")

#just plot oceano, legion, pseudo, pasteurellales
filtlppo=subset_taxa(filt, order==c("Legionellales", "Pseudomonadales", "Pasteurellales", "Oceanospirillales"))

#just oceano
filto=subset_taxa(filt, order== "Oceanospirillales")
filtoe=subset_taxa(filto, family== "Endozoicimonaceae")
filtoo=subset_taxa(filto, family=="Oceanospirillaceae")

filtto1=subset_samples(filto, Timepoint=="one")
filto2d=subset_samples(filto, Timepoint=="two" & T3_disease_state=="Diseased"& Dose_disease_state=="Diseased")
filto2h=subset_samples(filto, Timepoint=="two" & T3_disease_state=="Healthy")
filto3d=subset_samples(filto, Timepoint=="three" & final_disease_state=="Diseased"& Dose_disease_state=="Diseased")
filto3h=subset_samples(filto, Timepoint=="three" & final_disease_state=="Healthy")
filto2=subset_samples(filto, Timepoint=="two")

#just pseudomonadales
filtp=subset_taxa(filt, order== "Pasteurellales")
filtp1=subset_samples(filtp, Timepoint=="one")
filtp2d=subset_samples(filtp, Timepoint=="two" & T3_disease_state=="Diseased"& Dose_disease_state=="Diseased")
filtp2h=subset_samples(filtp, Timepoint=="two" & T3_disease_state=="Healthy")
filtp3d=subset_samples(filtp, Timepoint=="three" & final_disease_state=="Diseased"& Dose_disease_state=="Diseased")
filtp3h=subset_samples(filtp, Timepoint=="three" & final_disease_state=="Healthy")

#just legionellales
filtl=subset_taxa(filt, order== "Legionellales")
filtl1=subset_samples(filtl, Timepoint=="one")
filtl2d=subset_samples(filtl, Timepoint=="two" & T3_disease_state=="Diseased"& Dose_disease_state=="Diseased")
filtl2h=subset_samples(filtl, Timepoint=="two" & T3_disease_state=="Healthy")
filtl3d=subset_samples(filtl, Timepoint=="three" & final_disease_state=="Diseased"& Dose_disease_state=="Diseased")
filtl3h=subset_samples(filtl, Timepoint=="three" & final_disease_state=="Healthy")
filtl3=subset_samples(filtl, Timepoint=="three")




filtlppot2d=subset_samples(filtlppo, Timepoint=="two" & T3_disease_state=="Diseased"& Dose_disease_state=="Diseased")
filtlppot2h=subset_samples(filtlppo, Timepoint=="two" & T3_disease_state=="Healthy")


filtlppot3d=subset_samples(filtlppo, Timepoint=="three"& final_disease_state=="Diseased"& Dose_disease_state=="Diseased")
filtlppot3h=subset_samples(filtlppo, Timepoint=="three"& final_disease_state=="Healthy")

filto2dh=merge_samples(filto2, "T3_disease_state")
sample_data(filto2dh)$T3_disease_state <- levels(sample_data(filto2dh)$T3_disease_state)
filto2dh = transform_sample_counts(filto2dh, function(x) 100*x/sum(x))


p=plot_bar(filto2dh, y="Abundance", x="family", fill="family", facet_grid=~T3_disease_state)

p=plot_bar(filto2dh, y="Abundance", fill="family", facet_grid=~family)
p + geom_bar(aes(color=family, fill=family), stat="identity", position="stack")

filtdf=as.data.frame(otu_table(filt))
write.csv(as.data.frame(filtdf), file="ten_percent_filtered_otutable.csv")

#tree
top100=prune_taxa(names(sort(taxa_sums(filt), TRUE))[1:100], filt)

plot_tree(filtoe, color="family")

plot_tree(filtoe, size="abundance", color = "T3_disease_state", shape = "Timepoint")

plot_tree(, size="abundance", color = "T3_disease_state", 
          shape = "Timepoint", ladderize = "left") + coord_polar(theta = "y")


plot_heatmap(filto)
(p <- plot_heatmap(filto, method="NMDS", distance="bray", sample.label="T3_disease_state", taxa.label="family", sample.order="T3_disease_state"))

filtl3m=tax_glom(filtl3, taxrank="family")
df <- psmelt(filtl3m)

# group by Treatment and Family, calculate mean abundance and standard deviation
avgs <- ddply(df, c("Dose", "family"), summarise, mean=mean(Abundance))

# plot bar graph with standard deviation as error bars
means.barplot=qplot(x=Dose, y=mean, fill=family,
                       data=avgs, geom="bar", stat="identity",
                       position="dodge")
means.barplot


means.sem <- ddply(df, c("Dose", "family"), summarise,
                   mean=mean(Abundance), sem=sd(Abundance)/sqrt(length(Abundance)))
means.sem <- transform(means.sem, lower=mean-sem, upper=mean+sem)
means.barplot + geom_errorbar(aes(ymax=upper,
                                  ymin=lower),
                              position=position_dodge(0.9),
                              data=means.sem)

ggplot(avgs, aes(fill=treatment, x=Family, y=mean)) + 
  geom_bar(position='dodge') +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd))



save(filt, file="filtered_phyloseq.RData")
filtimp=load(file="filtered_phyloseq.RData")