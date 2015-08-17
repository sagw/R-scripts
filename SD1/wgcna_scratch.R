rm(list=ls())
require(WGCNA)
require(afex)
require(graphics)
setwd("~/Documents/SD1/SD1_analyses/")
c=read.csv("input_files//ten_percent_filtered_otutable.csv", h=T, row.names=1)

#information about samples
s=read.csv("input_files//WGCNA_traits.csv", h=T, as.is=T)
s=s[s$Timepoint %in% c("1", "2", "3", "0"),]
b=t(c)
d=cbind.data.frame(s,b)
row.names(d)=row.names(b)

chec=-cbind(row.names(d), d[1:1])

c=t(c)
rs = rowSums (c)
use = (rs > 0)
c = c [ use, ]


options(stringsAsFactors = FALSE) 

gsg = goodSamplesGenes(c, verbose = 3);
gsg$allOK 

c = c[gsg$goodSamples, gsg$goodGenes]

sampleTree = hclust(dist(c), method ="average")

par(cex = 0.6);
par(mar =c(0,4,2,0))
plot(sampleTree, main ="Sampleclusteringtodetectoutliers",sub=" ", xlab=
" ", cex.lab = 1.5,cex.axis= 1.5, cex.main = 2)

abline(h = 9000, col="red")

clust = cutreeStatic(sampleTree, cutHeight = 9000, minSize = 10)
table(clust)
keepSamples = (clust==1)
c = c[keepSamples, ]
nGenes =ncol(c)
nSamples =nrow(c)

s=s[!s$Sample.name=="SD1.CK14H3.CK14.P.T10", ]

traitColors = labels2colors(s)

sampleTree = hclust(dist(c), method ="average")

plotDendroAndColors(sampleTree, traitColors, groupLabels = names(s), 
                    main ="Sampledendrogramandtraitheatmap")


plotDendroAndColors(sampleTree, traitColors, groupLabels = names(s), dendroLabels=FALSE, 
                    main ="Sampledendrogramandtraitheatmap", guideAll=TRUE)

#determine soft threshold
powers =c(c(1:10), seq(from=12, to=20,by=2))

sft = pickSoftThreshold(c, powerVector = powers, networkType="signed", verbose = 5)

sizeGrWindow(9, 5)
par(mfrow =c(1,2))
cex1 = 0.9

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="SoftThreshold(power)",ylab="ScaleFreeTopologyModelFit,signedR^2",type="n",
main =paste("Scaleindependence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1, col="red")
abline(h=0.90,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="SoftThreshold(power)",
ylab="MeanConnectivity", type="n",main =paste("Meanconnectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5],labels=powers, cex=cex1,col="red")

#calculate module eigengenes create TOM
net = blockwiseModules(c, power = 16, TOMType ="signed", minModuleSize = 30,maxBlockSize=6000,
reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE,
saveTOMs = TRUE,saveTOMFileBase ="SD1_TOM" ,verbose = 3, networkType="signed")

table(net$colors)

#plot dendrogram with colors
sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
"Modulecolors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)


nGenes =ncol(c)
nSamples =nrow(c)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]

MEs0 = moduleEigengenes(c, moduleColors)$eigengenes
row.names(MEs0)=row.names(c)
MEs = orderMEs(MEs0)

row.names(s)=s$Sample.name
s=s[2:(NCOL(s))]
names(s)=colnames(s)
s1=apply(s, 2, as.numeric)
row.names(s1)=row.names(s)

## Find modules significantly correlated to traits
s1=as.data.frame(s1)
s2=read.csv("input_files/SD1_map.csv")
s2=s2[s2$Timepoint %in% c("one", "two", "three", "dose"),]
s2=s2[!s2$Sample.name=="SD1.CK14H3.CK14.P.T10", ]
row.names(s2)=s2$Sample.name
s2=s2[2:NCOL(s2)]

MESigs =matrix(NA,nrow=10,ncol=ncol(MEs))
colnames(MESigs)=colnames(MEs)

#mixed effects model for each module 
#repeated measures (including time as factor)
options(scipen=999)
for (i in 1:ncol(MESigs)) {
  exprvec = as.numeric(MEs[,i])
  combined=cbind.data.frame(s2, exprvec)
  combined=combined[!combined$Timepoint=="dose",]
 ## mixed=mixed(exprvec ~ Dose.Site+Site+Dose_disease_state+Dose.Site:Site+Dose.Site:Dose_disease_state+
                ##Site:Dose_disease_state+(1|Coral)+(1|Timepoint)+(1|Genotype),data = combined, method="LRT")
  mixed=mixed(exprvec ~ Dose.Site+Site+Dose_disease_state+Timepoint+Dose.Site:Site+Dose.Site:Dose_disease_state+
                Site:Dose_disease_state+Timepoint:Dose.Site+Timepoint:Site+Timepoint:Dose_disease_state
              +(1|Coral)+(1|Timepoint)+(1|Genotype),data = combined, method="LRT")
  mixed=as.data.frame(summary(mixed))
  pvalues=mixed$p.value
  pvalues=as.data.frame(pvalues)
  pvalues=apply(pvalues, 2, as.numeric)
  row.names(pvalues)=mixed$Effect
  pvalues=t(pvalues)
  DoseSitep=pvalues[,"Dose.Site"]
  Sitep=pvalues[,"Site"]
  Dose_disease_statep=pvalues[,"Dose_disease_state"]
  Timepointp=pvalues[,"Timepoint"]
  DoseSitexSitep=pvalues[,"Dose.Site:Site"]
  DoseSitexDose_disease_statep=pvalues[,"Dose.Site:Dose_disease_state"]
  SitexDose_disease_statep=pvalues[,"Site:Dose_disease_state"]
  TimepointxSitep=pvalues[,"Site:Timepoint"]
  TimepointxDose_disease_statep=pvalues[,"Dose_disease_state:Timepoint"]
  TimepointxDoseSitep=pvalues[,"Dose.Site:Timepoint"]
  MESigs[,i] = c(DoseSitep, Sitep, Dose_disease_statep, Timepointp, DoseSitexSitep, 
                  DoseSitexDose_disease_statep, SitexDose_disease_statep, TimepointxSitep, 
                 TimepointxDose_disease_statep, TimepointxDoseSitep)
}

rownames(MESigs) = c("DoseSite", "Site", "DoseDiseaseState", "Time", "DoseSitexSite", 
                      "DoseSitexDoseDiseaseState", "SitexDoseDiseaseState", "TimexSite",
                     "TimexDoseDiseaseState", "TimexDoseSite")
#get rid of nas (make <0.0001 numeric)
MESigs[is.na(MESigs)] = .00001


#just T2
MESigs =matrix(NA,nrow=6,ncol=ncol(MEs))
colnames(MESigs)=colnames(MEs)
options(scipen=999)
for (i in 1:ncol(MESigs)) {
  exprvec = as.numeric(MEs[,i])
  combined=cbind.data.frame(s2, exprvec)
  combined=combined[combined$Timepoint=="three",]
  mixed=mixed(exprvec ~ Dose.Site+Site+T3_disease_state+Dose.Site:Site+Dose.Site:T3_disease_state+
                Site:T3_disease_state+(1|Genotype),data = combined, method="LRT")
  mixed=as.data.frame(summary(mixed))
  pvalues=mixed$p.value
  pvalues=as.data.frame(pvalues)
  pvalues=apply(pvalues, 2, as.numeric)
  row.names(pvalues)=mixed$Effect
  pvalues=t(pvalues)
  DoseSitep=pvalues[,"Dose.Site"]
  Sitep=pvalues[,"Site"]
  T3_disease_statep=pvalues[,"T3_disease_state"]
  DoseSitexSitep=pvalues[,"Dose.Site:Site"]
  DoseSitexT3_disease_statep=pvalues[,"Dose.Site:T3_disease_state"]
  SitexT3_disease_statep=pvalues[,"Site:T3_disease_state"]
  MESigs[,i] = c(DoseSitep, Sitep, T3_disease_statep, DoseSitexSitep, 
                 DoseSitexT3_disease_statep, SitexT3_disease_statep)
}

rownames(MESigs) = c("DoseSite", "Site", "T3_disease_state", "DoseSitexSite", 
                     "DoseSitexT3_disease_state", "SitexT3_disease_state")



MESigs[is.na(MESigs)] = .0001


#make table with significant values colored.
require("plotrix")
MESigs=apply(MESigs, 2, as.numeric)

mm = as.matrix(MESigs, ncol = 3, row.names=TRUE)
mm2=t(mm)
col=ifelse(mm2>0.05, "white", "light pink")


par(mar = c(0.5, 0.5, 0.5, 0.5))
# create empty plot
plot(mm2, axes = FALSE, xlab = "", ylab = "", type = "n")

addtable2plot(x=0.1, y=0.25, table=mm2,lwd=par("lwd"),bty="n",bg=col,
              cex=0.6,box.col=col,text.col=par("fg"), 
              display.colnames=TRUE,display.rownames=TRUE,hlines=FALSE, vlines=FALSE, title=NULL)



for (i in 1:ncol(MESigs)) {
  exprvec = as.numeric(MEs[,i])
  glmsum=summary(glm(exprvec ~ Timepoint+Site+Dose_Site+Dose_disease_state+final_disease_state+T3_disease_state, data = s1))
  pvalues=coef(glmsum)[,4]
  pvalues=as.data.frame(pvalues)
  pvalues=t(pvalues)
  Timepointp=pvalues[,"Timepoint"]
  Sitep=pvalues[,"Site"]
  Dose_Sitep=pvalues[,"Dose_Site"]
  Dose_disease_statep=pvalues[,"Dose_disease_state"]
  final_disease_statep=pvalues[,"final_disease_state"]
  T3_disease_statep=pvalues[,"T3_disease_state"]
  MESigs[,i] = c(Timepointp, Sitep, Dose_Sitep, Dose_disease_statep, final_disease_statep, T3_disease_statep)
}

rownames(MESigs) = c("Timepoint_p","Site_p","Dose_Site_p","Dose_disease_state_p","final_disease_state_p","T3_disease_state_p")

#make boxplots of eigenvalues
traitdata=cbind(s2, MEs)
traitdata2=cbind(s1, MEs)
traitdata[traitdata == "Diseased"] = "D"
traitdata[traitdata == "Healthy"] = "H"
MEs=apply(MEs, 2, as.numeric)



sizeGrWindow(10,6)
T2traitdata=subset(traitdata, Timepoint=="two")
T3traitdata=subset(traitdata, Timepoint=="three")


boxplot(MEorange~ Dose*T3_disease_state, data=T2traitdata, main="Orange", ylab="Module Eigenvalues", cex.axis=0.75, las=3, subset=Timepoint=="two")


boxplot(MElightcyan~ Dose*T3_disease_state, data=T3traitdata, main="Light Cyan", ylab="Module Eigenvalues", cex.axis=0.5, las=3)

boxplot(MEdarkturquoise~ Dose_disease_state*Site*T3_disease_state, data=T3traitdata, main="Dark Turquoise", ylab="Module Eigenvalues", cex.axis=0.75, las=3)

boxplot(MElightgreen~ T3_disease_state*Timepoint, data=traitdata, main="Light Green", ylab="Module Eigenvalues", cex.axis=0.75, las=3)

boxplot(MEroyalblue~ T3_disease_state*Dose.Site*Timepoint, data=traitdata, main="Royal Blue", ylab="Module Eigenvalues", cex.axis=0.75, las=3)

boxplot(MEturquoise~ Site*T3_disease_state, data=T2traitdata, main="Turquoise", ylab="Module Eigenvalues", cex.axis=0.75, las=3)


boxplot(MEmidnightblue~ Dose*T3_disease_state, data=T3traitdata, main="Midnight Blue", ylab="Module Eigenvalues", cex.axis=0.75, las=3)

boxplot(MEyellow~ Dose*T3_disease_state, data=T3traitdata, main="Yellow", ylab="Module Eigenvalues", cex.axis=0.75, las=3)

boxplot(MEdarkgreen~ Site*T3_disease_state, data=T3traitdata, main="Dark green", ylab="Module Eigenvalues", cex.axis=0.75, las=3)

boxplot(MEred ~ Dose*T3_disease_state, data=T2traitdata, main="Red", ylab="Module Eigenvalues", cex.axis=0.75, las=3)

boxplot(MEcyan~ Site*T3_disease_state, data=T3traitdata, main="Cyan", ylab="Module Eigenvalues", cex.axis=0.75, las=3)

boxplot(MEdarkorange~ Dose*T3_disease_state, data=T3traitdata, main="Dark Orange", ylab="Module Eigenvalues", cex.axis=0.75, las=3)

boxplot(MEmagenta~ Dose*T3_disease_state, data=T3traitdata, main="Magenta", ylab="Module Eigenvalues", cex.axis=0.75, las=3)

boxplot(MEbrown~ Dose*T3_disease_state, data=T3traitdata, main="Brown", ylab="Module Eigenvalues", cex.axis=0.75, las=3)

boxplot(MEdarkgrey~ Dose*T3_disease_state, data=T3traitdata, main="Dark Grey", ylab="Module Eigenvalues", cex.axis=0.75, las=3)

boxplot(MEgreen~ T3_disease_state*Dose.Site*Timepoint, data=traitdata, main="Green", ylab="Module Eigenvalues", cex.axis=0.75, las=3)

boxplot(MEgrey60~ Dose*T3_disease_state, data=T3traitdata, main="Grey 60", ylab="Module Eigenvalues", cex.axis=0.75, las=3)

boxplot(MEblack~ Dose*T3_disease_state, data=T3traitdata, main="Black", ylab="Module Eigenvalues", cex.axis=0.75, las=3)

boxplot(MEpurple~ Dose*T3_disease_state, data=T3traitdata, main="Purple", ylab="Module Eigenvalues", cex.axis=0.75, las=3)

boxplot(MEblue~ Dose*T3_disease_state, data=T3traitdata, main="Blue", ylab="Module Eigenvalues", cex.axis=0.75, las=3)

boxplot(MEgreenyellow~ Dose*T3_disease_state, data=T3traitdata, main="Green yellow", ylab="Module Eigenvalues", cex.axis=0.75, las=3)
abline(h=0)

boxplot(MElightyellow~ Site*T3_disease_state, data=T3traitdata, main="Light yellow", ylab="Module Eigenvalues", cex.axis=0.75, las=3)
boxplot(MEdarkred~ Dose_disease_state*Dose.Site, data=traitdata, main="Dark Red", ylab="Module Eigenvalues", cex.axis=0.75, las=3)

boxplot(MEpink~ Site*T3_disease_state, data=T3traitdata, main="Pink", ylab="Module Eigenvalues", cex.axis=0.75, las=3)

#make heatmap.
sizeGrWindow(10,6)
textMatrix =paste(signif(MESigs, 2), sep ="")
dim(textMatrix)=dim(MESigs)
par(mar=c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = t(MESigs),xLabels =row.names(MESigs),yLabels =names(MEs),
ySymbols =names(MEs), colorLabels = FALSE,colors= greenWhiteRed(50),textMatrix = t(textMatrix),
setStdMargins = FALSE,cex.text= 0.5,zlim =c(-1,1),main =paste("Module-trait relationships"))



#plot relationships between modules
datME=moduleEigengenes(c, moduleColors)$eigengenes
row.names(datME)=row.names(c)
signif(cor(datME, use="p"), 2)
dissimME=(1-t(cor(datME, method="p")))/2
hclustdatME=hclust(as.dist(dissimME), method="average")

sizeGrWindow(8,7)
par(mfrow=c(1,1))
plot(hclustdatME, main="Clustering tree based off the module eigengenes")




#calculate Intramodular connectivity 
ADJ1=abs(cor(c, use="p"))^6
Alldegrees1=intramodularConnectivity(ADJ1, moduleColors)
head(Alldegrees1)

#calculate module eigengene based connectivity for each gene.
datKME=signedKME(c, datME, outputColumnName="MM.")
row.names(datKME)=gsub("X", "", row.names(datKME))
               
#plot module membership and intramodular connectivity
which.color="turquoise"
restrictGenes=moduleColors==which.color
verboseScatterplot(Alldegrees1$kWithin[ restrictGenes],
                   (datKME[restrictGenes, paste("MM.", which.color, sep="")])^6,
                   col=which.color,
                   xlab="Intramodular Connectivity",
                   ylab="(Module Membership)^6")


OTUsinmodule=geneInfo0$moduleColor=="magenta"
FilterOTUs=abs(datKME$MM.magenta)>.8
table(OTUsinmodule)
table(FilterOTUs)

dimnames(data.frame(c))[[2]][FilterGenes]


#names(colors) of the module
modNames=substring(names(MEs), 3)
               
geneModuleMembership=as.data.frame(cor(c, MEs, use="p"))
               
MMPvalue=as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))               
               
names(geneModuleMembership)=paste("MM", modNames, sep="")               
names(MMPvalue)=paste("p.MM", modNames, sep="")               

            
          
               
t=read.csv("input_files/rep_set_tax_assignments.csv", row.names=1)


dim(t)  
probes=colnames(c)               
probes2annot=match(probes, row.names(t))               

sum(is.na(probes2annot))               

#createdf
geneInfo0=data.frame(TAX=t$Tax[probes2annot],
                     moduleColor=moduleColors, row.names=probes)    
               
geneInfo1=merge(geneInfo0, datKME, by="row.names")
               
write.csv(geneInfo1, file="Output_files/WGCNA_geneInfo_klms_signedglm.csv")         


#export to cytoscape
TOM=TOMsimilarityFromExpr(c, power=7)
modules=c("lightcyan")                
               
probes=colnames(c)               
inModule=is.finite(match(moduleColors, modules))
modProbes=probes[inModule]
modOtus=t$Tax[match(modProbes, row.names(t))]               
modTOM=TOM[inModule, inModule]               

dimnames(modTOM)=list(modProbes, modProbes)

cyt=exportNetworkToCytoscape(modTOM, 
edgeFile=paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".csv", sep=","),
nodeFile=paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".csv", sep=","),              
weighted=TRUE,
threshold=0.02,
nodeNames=modProbes,
altNodeNames=modOtus,
#nodeAttr=moduleColors[inModule],
nodeAttr=datKME$MM.lightcyan[inModule])                          
       