#!/usr/bin/env Rscript
# exomepeak Script 2
# R script
# Define parameters and load library
## IP vs input
library("DEGseq")
directory="/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/mouse/DEGseq"
setwd(directory)
featureCount="mouse_RNA-seq_featureCount.txt"
colCount <- ncol(read.table(file=featureCount, nrows=1))
countMatrix <- readGeneExp(file=featureCount, geneCol=1, valCol=c(7:colCount))

DEGexp(geneExpMatrix1=countMatrix, geneCol1=1, expCol1=c(2), groupLabel1="ip",
    geneExpMatrix2=countMatrix, geneCol2=1, expCol2=c(3), groupLabel2="input",
    method="MARS", outputDir=paste(directory, "mESCsShContD0", sep="/"))

DEGexp(geneExpMatrix1=countMatrix, geneCol1=1, expCol1=c(4), groupLabel1="ip",
    geneExpMatrix2=countMatrix, geneCol2=1, expCol2=c(5), groupLabel2="input",
    method="MARS", outputDir=paste(directory, "mESCsShContD6", sep="/"))

DEGexp(geneExpMatrix1=countMatrix, geneCol1=1, expCol1=c(6), groupLabel1="ip",
    geneExpMatrix2=countMatrix, geneCol2=1, expCol2=c(7), groupLabel2="input",
    method="MARS", outputDir=paste(directory, "mESCsShSetD2D0", sep="/"))

DEGexp(geneExpMatrix1=countMatrix, geneCol1=1, expCol1=c(8), groupLabel1="ip",
    geneExpMatrix2=countMatrix, geneCol2=1, expCol2=c(9), groupLabel2="input",
    method="MARS", outputDir=paste(directory, "mESCsShSetD2D6", sep="/"))

DEGexp(geneExpMatrix1=countMatrix, geneCol1=1, expCol1=c(10), groupLabel1="ip",
    geneExpMatrix2=countMatrix, geneCol2=1, expCol2=c(11), groupLabel2="input",
    method="MARS", outputDir=paste(directory, "mEFsShCont", sep="/"))

DEGexp(geneExpMatrix1=countMatrix, geneCol1=1, expCol1=c(12), groupLabel1="ip",
    geneExpMatrix2=countMatrix, geneCol2=1, expCol2=c(13), groupLabel2="input",
    method="MARS", outputDir=paste(directory, "mEFsShSetD2", sep="/"))


## input vs input
directory="/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/mouse/DEGseq/DeGenes"
setwd(directory)
featureCount="../mouse_RNA-seq_featureCount.txt"
colCount <- ncol(read.table(file=featureCount, nrows=1))
countMatrix <- readGeneExp(file=featureCount, geneCol=1, valCol=c(7:colCount))

DEGexp(geneExpMatrix1=countMatrix, geneCol1=1, expCol1=c(4), groupLabel1="mESCsShContD6",
    geneExpMatrix2=countMatrix, geneCol2=1, expCol2=c(3), groupLabel2="mESCsShContD0",
    method="MARS", outputDir=paste(directory, "mESCsShContD6", sep="/"))

DEGexp(geneExpMatrix1=countMatrix, geneCol1=1, expCol1=c(7), groupLabel1="mESCsShSetD2D0",
    geneExpMatrix2=countMatrix, geneCol2=1, expCol2=c(3), groupLabel2="mESCsShContD0",
    method="MARS", outputDir=paste(directory, "mESCsShSetD2D0", sep="/"))

DEGexp(geneExpMatrix1=countMatrix, geneCol1=1, expCol1=c(9), groupLabel1="mESCsShSetD2D6",
    geneExpMatrix2=countMatrix, geneCol2=1, expCol2=c(3), groupLabel2="mESCsShContD0",
    method="MARS", outputDir=paste(directory, "mESCsShSetD2D6", sep="/"))

DEGexp(geneExpMatrix1=countMatrix, geneCol1=1, expCol1=c(13), groupLabel1="mEFsShSetD2",
    geneExpMatrix2=countMatrix, geneCol2=1, expCol2=c(11), groupLabel2="mEFsShCont",
    method="MARS", outputDir=paste(directory, "mEFsShSetD2", sep="/"))

