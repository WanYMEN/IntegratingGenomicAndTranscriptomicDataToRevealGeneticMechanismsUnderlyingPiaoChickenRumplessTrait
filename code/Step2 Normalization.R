## Gene filtration and expression normalization based on the three expression datasets, i.e. HTSeq union exon model, HTSeq whole gene model and Cufflinks, for the paper "Integrating Genomic and Transcriptomic Data to Reveal Genetic Mechanisms Underlying Piao Chicken Rumpless Trait".

##the path that stores the code folder 
setwd("IntegratingGenomicAndTranscriptomicDataToRevealGeneticMechanismsUnderlyingPiaoChickenRumplessTrait/code")

rm(list=ls())
options(stringsAsFactors=FALSE)

##load the lower bound expression matrix from Cufflinks
dataLB <- read.table("../data/LowerBoundFPKM.txt", sep="\t", header=TRUE)
row.names(dataLB) <- dataLB$gene
ExpLB <- dataLB[, -1]
dim(ExpLB)

keep.LB.CA <- grepl("Piao", colnames(ExpLB))
ExpLB.CA <- ExpLB[, keep.LB.CA]
ExpLB.CO <- ExpLB[, !keep.LB.CA]
dim(ExpLB.CA)
dim(ExpLB.CO)

## Filter genes based on the lower bound expression matrix from Cufflinks. If the lower bound of expression is > 0 in 80% of samples, we keep the genes. We filter in Piao chicken and control chickens separately.
passvec.LB.CA <- apply(ExpLB.CA>0, 1, sum)>(0.8*sum(keep.LB.CA))
passvec.LB.CO <- apply(ExpLB.CO>0, 1, sum)>(0.8*sum(!keep.LB.CA))

ExpLB.filter <- ExpLB[union(rownames(ExpLB[passvec.LB.CA, ]), rownames(ExpLB[passvec.LB.CO, ])), ]

## For the HTSeq counts data of whole gene model and exon union model, we both used a threshold of observing at least 10 fragments in 80% of the samples in Piao chicken and control chickens separately.

##HTSeq.gene
htgene <- read.table("../data/HTSeqGene.txt", sep="\t", header=TRUE)
row.names(htgene) <- htgene$gene
preExp.HTSeq.gene <- htgene[, -1]

keep.htgene.CA <- grepl("Piao", colnames(preExp.HTSeq.gene))
preExp.HTSeq.gene.CA <- preExp.HTSeq.gene[, keep.htgene.CA]
preExp.HTSeq.gene.CO <- preExp.HTSeq.gene[, !keep.htgene.CA]

passvec.CA.gene <- apply(preExp.HTSeq.gene.CA, 1, quantile, 0.8) > 10
passvec.CO.gene <- apply(preExp.HTSeq.gene.CO, 1, quantile, 0.8) > 10
preExp.HTSeq.filter.gene <- preExp.HTSeq.gene[union(rownames(preExp.HTSeq.gene[passvec.CA.gene, ]), rownames(preExp.HTSeq.gene[passvec.CO.gene, ])), ]

##HTSeq.exon
htexon <- read.table("../data/HTSeqExon.txt", sep="\t", header=TRUE)
row.names(htexon) <- htexon$gene
preExp.HTSeq.exon <- htexon[, -1]

keep.htexon.CA <- grepl("Piao", colnames(preExp.HTSeq.exon))
preExp.HTSeq.exon.CA <- preExp.HTSeq.exon[, keep.htexon.CA]
preExp.HTSeq.exon.CO <- preExp.HTSeq.exon[, !keep.htexon.CA]

passvec.CA.exon <- apply(preExp.HTSeq.exon.CA, 1, quantile, 0.8) > 10
passvec.CO.exon <- apply(preExp.HTSeq.exon.CO, 1, quantile, 0.8) > 10
preExp.HTSeq.filter.exon <- preExp.HTSeq.exon[union(rownames(preExp.HTSeq.exon[passvec.CA.exon, ]), rownames(preExp.HTSeq.exon[passvec.CO.exon, ])), ]

##load the "cqn" package
library(cqn)

## gene lengths and GC content file from chicken reference genome
gc <- read.table("../data/GC_lengths.tsv", sep="\t", header=TRUE)
## Keep only the genes with length > 200; the others are hard to assess and would be difficult to map to anyway
gc.len <- gc[gc[, 1]>200, ]

keepvec.htexon <- intersect(rownames(gc.len), rownames(preExp.HTSeq.filter.exon))
geneAnno.htexon <- gc.len[keepvec.htexon, ]
Exp.HTSeq.filter.exon <- preExp.HTSeq.filter.exon[keepvec.htexon, ]

cqn.htexon <- cqn(Exp.HTSeq.filter.exon, lengths = as.numeric(geneAnno.htexon[, 1]), x = as.numeric(geneAnno.htexon[, 2]), lengthMethod=c("smooth"), sqn=FALSE)
Exp.HTSeq.exon.cqn <- cqn.htexon$y + cqn.htexon$offset

## Intersect to get normalized FPKM matrix with expression values that pass filtering criteria
htscint <- intersect(rownames(preExp.HTSeq.filter.gene), intersect(rownames(Exp.HTSeq.exon.cqn), rownames(ExpLB.filter)))

datExpr <- Exp.HTSeq.exon.cqn[match(htscint, rownames(Exp.HTSeq.exon.cqn)), ]
save(datExpr, file="../output/NormalizedExpressionData.Rdata")
