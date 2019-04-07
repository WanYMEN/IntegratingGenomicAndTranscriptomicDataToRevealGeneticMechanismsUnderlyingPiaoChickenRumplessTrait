## adjust unwanted biological and technical covariates for the paper "Integrating Genomic and Transcriptomic Data to Reveal Genetic Mechanisms Underlying Piao Chicken Rumpless Trait"

## the path that stores the code folder 
setwd("IntegratingGenomicAndTranscriptomicDataToRevealGeneticMechanismsUnderlyingPiaoChickenRumplessTrait/code")

rm(list=ls())
options(stringsAsFactors=FALSE)

## load filtered and normalized data
load("../data/NormalizedExpressionData.Rdata")

## load sample information file
covar <- read.table("../data/Covariants.txt", sep="\t", header=TRUE)

condition <- as.numeric(as.factor(covar$condition))-1
age <- as.numeric(as.factor(covar$age))
breed <- as.numeric(as.factor(covar$breed))
lane <- as.numeric(as.factor(covar$lane))
regvars <- as.data.frame(cbind(condition, age, breed, lane))

## Run the regression and make the adjusted FPKM matrix
datExpr.reg <- matrix(NA, nrow=nrow(datExpr), ncol = ncol(datExpr))
rownames(datExpr.reg) <- rownames(datExpr)
colnames(datExpr.reg) <- colnames(datExpr)
coefmat <- matrix(NA, nrow = nrow(datExpr), ncol = ncol(regvars)+1)

for (i in 1:nrow(datExpr)) {
  lmmod1 <- lm(as.numeric(datExpr[i,])~condition+age+breed+lane, data = regvars)
  coef <- coef(lmmod1)
  coefmat[i,] <- coef
  datExpr.reg[i,] <- datExpr[i,] - coef["age"]*regvars[,"age"] - coef["breed"]*regvars[,"breed"] - coef["lane"]*regvars[,"lane"]
}

## This datExpr.reg is now a technical variable corrected matrix.
rownames(datExpr.reg) <- rownames(datExpr)
colnames(datExpr.reg) <- colnames(datExpr)

save(datExpr.reg, covar, file="../output/AdjustedExpressionData.Rdata")
