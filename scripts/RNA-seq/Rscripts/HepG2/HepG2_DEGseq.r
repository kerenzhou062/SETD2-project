#!/usr/bin/env Rscript
# exomepeak Script 2
# R script
# Define parameters and load library
library("DEGseq")
## IP vs input
directory="/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/HepG2/DEGseq"
setwd(directory)
featureCount="HepG2_RNA-seq_featureCount.txt"
colCount <- ncol(read.table(file=featureCount, nrows=1))
countMatrix <- readGeneExp(file=featureCount, geneCol=1, valCol=c(7:colCount))

DEGexp(geneExpMatrix1=countMatrix, geneCol1=1, expCol1=c(2,3,4), groupLabel1="ip",
    geneExpMatrix2=countMatrix, geneCol2=1, expCol2=c(5,6,7), groupLabel2="input",
    method="MARS", outputDir=paste(directory, "shCont", sep="/"))

DEGexp(geneExpMatrix1=countMatrix, geneCol1=1, expCol1=c(8,9,10), groupLabel1="ip",
    geneExpMatrix2=countMatrix, geneCol2=1, expCol2=c(11,12,13), groupLabel2="input",
    method="MARS", outputDir=paste(directory, "shSetD2", sep="/"))

DEGexp(geneExpMatrix1=countMatrix, geneCol1=1, expCol1=c(14,15,16), groupLabel1="ip",
    geneExpMatrix2=countMatrix, geneCol2=1, expCol2=c(17,18,19), groupLabel2="input",
    method="MARS", outputDir=paste(directory, "shM14", sep="/"))

DEGexp(geneExpMatrix1=countMatrix, geneCol1=1, expCol1=c(20,21,22), groupLabel1="ip",
    geneExpMatrix2=countMatrix, geneCol2=1, expCol2=c(23,24,25), groupLabel2="input",
    method="MARS", outputDir=paste(directory, "shM3", sep="/"))

DEGexp(geneExpMatrix1=countMatrix, geneCol1=1, expCol1=c(26,27,28), groupLabel1="ip",
    geneExpMatrix2=countMatrix, geneCol2=1, expCol2=c(29,30,31), groupLabel2="input",
    method="MARS", outputDir=paste(directory, "shWTAP", sep="/"))


## input vs input
directory="/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/HepG2/DEGseq/DeGenes"
setwd(directory)
featureCount="../HepG2_RNA-seq_featureCount.txt"
colCount <- ncol(read.table(file=featureCount, nrows=1))
countMatrix <- readGeneExp(file=featureCount, geneCol=1, valCol=c(7:colCount))

DEGexp(geneExpMatrix1=countMatrix, geneCol1=1, expCol1=c(11,12,13), groupLabel1="shSetD2",
    geneExpMatrix2=countMatrix, geneCol2=1, expCol2=c(5,6,7), groupLabel2="shCont",
    method="MARS", outputDir=paste(directory, "shSetD2", sep="/"))

DEGexp(geneExpMatrix1=countMatrix, geneCol1=1, expCol1=c(17,18,19), groupLabel1="shM14",
    geneExpMatrix2=countMatrix, geneCol2=1, expCol2=c(5,6,7), groupLabel2="shCont",
    method="MARS", outputDir=paste(directory, "shM14", sep="/"))

DEGexp(geneExpMatrix1=countMatrix, geneCol1=1, expCol1=c(23,24,25), groupLabel1="shM3",
    geneExpMatrix2=countMatrix, geneCol2=1, expCol2=c(5,6,7), groupLabel2="shCont",
    method="MARS", outputDir=paste(directory, "shM3", sep="/"))

DEGexp(geneExpMatrix1=countMatrix, geneCol1=1, expCol1=c(29,30,31), groupLabel1="shWTAP",
    geneExpMatrix2=countMatrix, geneCol2=1, expCol2=c(5,6,7), groupLabel2="shCont",
    method="MARS", outputDir=paste(directory, "shWTAP", sep="/"))


