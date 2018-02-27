#!/usr/bin/env Rscript
# exomepeak Script 2
# R script
# Define parameters and load library
## IP vs input
library("DEGseq")
directory="/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/Hela/DEGseq"
setwd(directory)
featureCount="Hela_RNA-seq_featureCount.txt"
colCount <- ncol(read.table(file=featureCount, nrows=1))
countMatrix <- readGeneExp(file=featureCount, geneCol=1, valCol=c(7:colCount))

DEGexp(geneExpMatrix1=countMatrix, geneCol1=1, expCol1=c(4), groupLabel1="ip",
    geneExpMatrix2=countMatrix, geneCol2=1, expCol2=c(7), groupLabel2="input",
    method="MARS", outputDir=paste(directory, "shCont", sep="/"))

DEGexp(geneExpMatrix1=countMatrix, geneCol1=1, expCol1=c(10), groupLabel1="ip",
    geneExpMatrix2=countMatrix, geneCol2=1, expCol2=c(13), groupLabel2="input",
    method="MARS", outputDir=paste(directory, "shSetD2", sep="/"))

DEGexp(geneExpMatrix1=countMatrix, geneCol1=1, expCol1=c(16), groupLabel1="ip",
    geneExpMatrix2=countMatrix, geneCol2=1, expCol2=c(19), groupLabel2="input",
    method="MARS", outputDir=paste(directory, "shM14", sep="/"))

DEGexp(geneExpMatrix1=countMatrix, geneCol1=1, expCol1=c(22), groupLabel1="ip",
    geneExpMatrix2=countMatrix, geneCol2=1, expCol2=c(25), groupLabel2="input",
    method="MARS", outputDir=paste(directory, "shM3", sep="/"))

DEGexp(geneExpMatrix1=countMatrix, geneCol1=1, expCol1=c(28), groupLabel1="ip",
    geneExpMatrix2=countMatrix, geneCol2=1, expCol2=c(31), groupLabel2="input",
    method="MARS", outputDir=paste(directory, "shWTAP", sep="/"))


## input vs input
directory="/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/Hela/DEGseq/DeGenes"
setwd(directory)
featureCount="../Hela_RNA-seq_featureCount.txt"
colCount <- ncol(read.table(file=featureCount, nrows=1))
countMatrix <- readGeneExp(file=featureCount, geneCol=1, valCol=c(7:colCount))

DEGexp(geneExpMatrix1=countMatrix, geneCol1=1, expCol1=c(13), groupLabel1="shSetD2",
    geneExpMatrix2=countMatrix, geneCol2=1, expCol2=c(7), groupLabel2="shCont",
    method="MARS", outputDir=paste(directory, "shSetD2", sep="/"))

DEGexp(geneExpMatrix1=countMatrix, geneCol1=1, expCol1=c(19), groupLabel1="shM14",
    geneExpMatrix2=countMatrix, geneCol2=1, expCol2=c(7), groupLabel2="shCont",
    method="MARS", outputDir=paste(directory, "shM14", sep="/"))

DEGexp(geneExpMatrix1=countMatrix, geneCol1=1, expCol1=c(25), groupLabel1="shM3",
    geneExpMatrix2=countMatrix, geneCol2=1, expCol2=c(7), groupLabel2="shCont",
    method="MARS", outputDir=paste(directory, "shM3", sep="/"))

DEGexp(geneExpMatrix1=countMatrix, geneCol1=1, expCol1=c(31), groupLabel1="shWTAP",
    geneExpMatrix2=countMatrix, geneCol2=1, expCol2=c(7), groupLabel2="shCont",
    method="MARS", outputDir=paste(directory, "shWTAP", sep="/"))

