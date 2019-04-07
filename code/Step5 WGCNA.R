## the weighted gene co-expression network analysis (WGCNA) for the paper "Integrating Genomic and Transcriptomic Data to Reveal Genetic Mechanisms Underlying Piao Chicken Rumpless Trait"

##the path that stores the code folder 
setwd("IntegratingGenomicAndTranscriptomicDataToRevealGeneticMechanismsUnderlyingPiaoChickenRumplessTrait/code")

rm(list=ls())
options(stringsAsFactors=FALSE)
library(WGCNA)

load("../data/AdjustedExpressionData.Rdata")

tExpr = t(datExpr.reg) 
gsg = goodSamplesGenes(tExpr, verbose = 3);
gsg$allOK

##using multiple parallel execution
enableWGCNAThreads()

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(tExpr, powerVector = powers, verbose = 5, networkType="signed", corFnc="bicor")

pdf(file = "../output/automatic_soft_thresholding_power.pdf", width = 9, height = 5)
par(mfrow = c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit,signed R^2", type="n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, col="red")
abline(h=0.80, col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, col="red")
dev.off()

softPower = 12

## Run an automated network analysis; it takes about one hour, and need large memory
net = blockwiseModules(tExpr, power = softPower, minModuleSize = 30, mergeCutHeight = 0.25, numericLabels = TRUE, corType="bicor", networkType="signed", pamStage=FALSE, pamRespectsDendro = FALSE, saveTOMs = FALSE, verbose = 3, maxBlockSize=20000)

save(net, file = "../output/netblockwiseModules_auto.RData")

nGenes = ncol(tExpr)
nSamples = nrow(tExpr)
condition <- as.numeric(as.factor(covar[rownames(tExpr), ]$condition))-1
moduleTraitCor <- cor(net$MEs, condition, use = "p")
moduleTraitCor <- moduleTraitCor[paste("ME", 0:11, sep = ""), ]
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
rowColors = labels2colors(0:11)

# Display the correlation and p-values within a heatmap plot
pdf(file = "../output/Module_trait_relationships.pdf", width = 3, height = 7);
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
labeledHeatmap(Matrix = as.matrix(moduleTraitCor), xLabels = c("Piao"), yLabels = paste("M", rowColors, sep=""), colorLabels = FALSE, ySymbols = names(moduleTraitCor), colors = greenWhiteRed(50), textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 1, zlim = c(-1, 1), main = paste("Module-trait relationships"), xLabelsAngle = 0)
dev.off()

geneTraitSignificance = as.data.frame(cor(tExpr, condition, use = "p"));
names(geneTraitSignificance) ="GS"

ADJ=abs(cor(tExpr, use="p"))^12
Alldegrees=intramodularConnectivity(ADJ, net$colors)
head(Alldegrees)

datKME=signedKME(tExpr, net$MEs, outputColumnName = "kME", corFnc = "bicor")
absKME= abs(datKME)

kmecol = colnames(absKME)

# Calculate topological overlap anew
dissTOM = 1-TOMsimilarityFromExpr(tExpr, power = softPower, networkType = "signed", corType = "bicor");

TOM = 1-dissTOM

sym = read.delim(file = "../data/Galgal4.79.EnsId.SymbolName.txt", sep="\t", header=F)

## produce Cytoscape edges and nodes input files
for (i in c(1:length(kmecol)))
{
which.col = kmecol[i]
FilterGenes= abs(geneTraitSignificance)> 0.2 & absKME[, which.col] >0.8 & moduleColors == substring(which.col, 4)
hubGenes = dimnames(data.frame(tExpr))[[2]][FilterGenes]
hubvalue = cbind(hubGenes, geneTraitSignificance[FilterGenes, ], absKME[hubGenes, which.col], Alldegrees$kWithin[FilterGenes])

topgene = hubvalue[order(as.numeric(hubvalue[,4]), decreasing = TRUE),1][1:50]
topfilter = match(topgene, rownames(FilterGenes))
modTOM = TOM[topfilter, topfilter]
dimnames(modTOM) = list(topgene, topgene)
modGenes = sym[match(topgene, sym[,1]),2]

cyt = exportNetworkToCytoscape(modTOM,
  edgeFile = paste("CytoscapeInput-edges-", substring(which.col,2), ".top50.txt", sep=""),
  nodeFile = paste("CytoscapeInput-nodes-", substring(which.col,2), ".top50.txt", sep=""),
  weighted = TRUE,
  threshold = 0.1,
  nodeNames = topgene,
  altNodeNames = modGenes,
  nodeAttr = moduleColors[topfilter])
}
