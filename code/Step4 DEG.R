##the path that stores the code folder
setwd("IntegratingGenomicAndTranscriptomicDataToRevealGeneticMechanismsUnderlyingPiaoChickenRumplessTrait/code")

rm(list=ls())
options(stringsAsFactors=FALSE)

##identify DEGs by a linear model
load("../data/NormalizedExpressionData.Rdata")

covar <- read.table("../data/Covariants.txt", sep="\t", header=TRUE)

condition <- as.numeric(as.factor(covar$condition))-1
age <- as.numeric(as.factor(covar$age))
breed <- as.numeric(as.factor(covar$breed))
lane <- as.numeric(as.factor(covar$lane))
regvars <- as.data.frame(cbind(condition, age, breed, lane))

result <- c()
filt <- c()
final <- c()

for (i in 1:nrow(datExpr)) {
  expr <- as.numeric(as.character(datExpr[i,]))
  lmmodl <- lm(expr ~ condition + age + breed + lane)
  sm <- summary(lmmodl)
  adjR <- sm$adj.r.squared
  cf <- coef(sm)
  coefficient <- cf["condition", 1]
  pvalue <- cf["condition", 4]
  sub <- cbind(rownames(datExpr)[i], adjR, coefficient,  pvalue)
  result <- rbind(result, sub)
}

filt <- result[!is.na(result[, 4]), ]

if(length(filt)>0) {
   colnames(filt) <- c("Gene", "Adj.R2", "Coefficient", "Pvalue")
   FDR <- p.adjust(filt[, 4], "fdr")
   final <- data.frame(filt, FDR)
}

sig <- final[as.numeric(as.character(final$FDR)) <0.05, ]

write.table(sig, "../output/lm.DEG.FDR0.05.txt", row.names = FALSE, sep = "\t", quote = FALSE)


##identify DEGs by DESeq2
library(DESeq2)

data = read.table("../data/HTSeqExon.txt", sep="\t", header=T, stringsAsFactors=F, row.names=1)
head(data)
data=data[-2, ]

trait = data.frame(Tail=c(rep("N",9), rep("Y",12)))

ds <- DESeqDataSetFromMatrix(countData=data, colData=trait, design=~Tail)

ds<- DESeq(ds)
res <- results(ds, c("Tail","N","Y"))

res <- res[!is.na(res$pvalue), ]

sig<- res[which(res$padj < 0.05), ]
sig <-  as.data.frame(sig[order(sig$padj), ])
dim(sig)

write.csv(sig, file="../output/DESeq2v1_20_0.0.05FDR.csv", quote=F)
