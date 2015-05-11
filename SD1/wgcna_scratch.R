rm(list=ls())
require(WGCNA)
setwd("~/Documents/SD1/SD1_analyses/")
c=read.csv("input_files//ten_percent_filtered_otutable.csv", h=T, row.names=1)

#information about samples
s=read.csv("input_files//WGCNA_traits.csv", h=T, as.is=T)
s=s[s$Timepoint %in% c("1", "2", "3", "0"),]
b=t(c)
d=cbind.data.frame(s,b)
row.names(d)=row.names(b)

check<-cbind(row.names(d), d[1:1])



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

plotDendroAndColors(sampleTree, traitColors, groupLabels = names(s), 
                    main ="Sampledendrogramandtraitheatmap")

powers =c(c(1:10),seq(from = 12, to=20,by=2))

sft = pickSoftThreshold(c, powerVector = powers, verbose = 5)

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

net = blockwiseModules(c, power = 7, TOMType ="unsigned", minModuleSize = 30,maxBlockSize=6000,
reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE,
saveTOMs = TRUE,saveTOMFileBase ="SD1_TOM" ,verbose = 3)

table(net$colors)

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
MEs = orderMEs(MEs0)
s=s[2:(NCOL(s))]
names(s)=colnames(s)
s1=apply(s, 2, as.numeric)
moduleTraitCor=cor(MEs, s1, use="p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)
textMatrix =paste(signif(moduleTraitCor, 2),"\n(",signif(moduleTraitPvalue, 1),")", sep ="")
dim(textMatrix)=dim(moduleTraitCor)
par(mar=c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,xLabels =names(s),yLabels =names(MEs),
ySymbols =names(MEs), colorLabels = FALSE,colors= greenWhiteRed(50),textMatrix = textMatrix,
setStdMargins = FALSE,cex.text= 0.5,zlim =c(-1,1),main =paste("Module-trait relationships"))

#define variable T3_disease_State
T3_disease_state=as.data.frame(s$T3_disease_state)
names(T3_disease_state)="T3_disease_state"
               
#names(colors) of the module
modNames=substring(names(MEs), 3)
               
geneModuleMembership=as.data.frame(cor(c, MEs, use="p"))
               
MMPvalue=as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))               
               
names(geneModuleMembership)=paste("MM", modNames, sep="")               
names(MMPvalue)=paste("p.MM", modNames, sep="")               

geneTraitSignificance=as.data.frame(cor(c, T3_disease_state, use="p"))               
GSPvalue=as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))               
               
names(geneTraitSignificance)=paste("GS.", names(T3_disease_state), sep="")              
names(GSPvalue)=paste("p.GS", names(T3_disease_state), sep="")               
               
module="greenyellow"
column=match(module, modNames)
moduleGenes=moduleColors==module               

sizeGrWindow(7, 7)
par(mfrow=c(1,1))               
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab=paste("Module Membership in", module, "module"),
                   ylab="Gene significance for T3 disease state",
                   main=paste("Module membership vs. gene significance\n"),
                   cex.main=1.2, cex.lab=1.2, cex.axis=1.2, col=module)

colnames(c)[moduleColors=="greenyellow"]
               
               
t=read.csv("input_files/rep_set_tax_assignments.csv", row.names=1)
               
dim(t)  
probes=colnames(c)               
probes2annot=match(probes, row.names(t))               

sum(is.na(probes2annot))               

#createdf
geneInfo0=data.frame(TAX=t$Tax[probes2annot],
                     moduleColor=moduleColors,
                     geneTraitSignificance,
                     GSPvalue)    
               
#order modules by significance for t3 disease state               
modOrder=order(-abs(cor(MEs, T3_disease_state, use="p")))

#add module membership info in order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames=names(geneInfo0)
  geneInfo0=data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                       MMPvalue[,modOrder[mod]])
  names(geneInfo0)=c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                     paste("p.MM.", modNames[modOrder[mod]], sep=""))                     
}

#order genes in geneInfo variable by module color then significance               
geneOrder=order(geneInfo0$moduleColor, -abs(geneInfo0$GS.T3_disease_state))
geneInfo=geneInfo0[geneOrder,]
               
write.csv(geneInfo, file="Output_files/WGCNA_geneInfo.csv")         

#export to cytoscape
TOM=TOMsimilarityFromExpr(c, power=7)
modules=c("greenyellow", "magenta")               
               
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
nodeAttr=moduleColors[inModule])                          
                             
datME=moduleEigengenes(c, moduleColors)$eigengenes
row.names(datME)=row.names(c)
sizeGrWindow(8,7)
which.module="turquoise"
ME=datME[, paste("ME", which.module, sep="")]
ME=as.data.frame(ME)
row.names(ME)=row.names(datME)
MEdat=cbind(s, ME)

signif(cor(datME, use="p"), 2)

dissimME=(1-t(cor(datME, method="p")))/2
hclustdatME=hclust(as.dist(dissimME), method="average")

par(mfrow=c(1,1))
plot(hclustdatME, main="Clustering tree based off the module eigengenes")