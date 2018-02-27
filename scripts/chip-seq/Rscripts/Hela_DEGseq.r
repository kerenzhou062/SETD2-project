#!/usr/bin/env Rscript
# exomepeak Script 2
# R script
# Define parameters and load library
library("DEGseq")
directory="/data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/Hela/DEGseq"
setwd(directory)
featureCount="Hela_ChIP-seq_featureCount.txt"
colCount <- ncol(read.table(file=featureCount, nrows=1))
countMatrix <- readGeneExp(file=featureCount, geneCol=1, valCol=c(7:colCount))

DEGexp(geneExpMatrix1=countMatrix, geneCol1=1, expCol1=c(3), groupLabel1="ip",
    geneExpMatrix2=countMatrix, geneCol2=1, expCol2=c(5), groupLabel2="input",
    method="MARS", outputDir=paste(directory, "shCont", sep="/"))

DEGexp(geneExpMatrix1=countMatrix, geneCol1=1, expCol1=c(7), groupLabel1="ip",
    geneExpMatrix2=countMatrix, geneCol2=1, expCol2=c(9), groupLabel2="input",
    method="MARS", outputDir=paste(directory, "shSetD2", sep="/"))

DEGexp(geneExpMatrix1=countMatrix, geneCol1=1, expCol1=c(10), groupLabel1="ip",
    geneExpMatrix2=countMatrix, geneCol2=1, expCol2=c(11), groupLabel2="input",
    method="MARS", outputDir=paste(directory, "shM14", sep="/"))

DEGexp(geneExpMatrix1=countMatrix, geneCol1=1, expCol1=c(12), groupLabel1="ip",
    geneExpMatrix2=countMatrix, geneCol2=1, expCol2=c(13), groupLabel2="input",
    method="MARS", outputDir=paste(directory, "shM3", sep="/"))

DEGexp(geneExpMatrix1=countMatrix, geneCol1=1, expCol1=c(14), groupLabel1="ip",
    geneExpMatrix2=countMatrix, geneCol2=1, expCol2=c(15), groupLabel2="input",
    method="MARS", outputDir=paste(directory, "shWTAP", sep="/"))

#DEGexp(geneExpMatrix1=geneExpMatrix1, expCol1=2, groupLabel1="ip1",
# geneExpMatrix2=geneExpMatrix2, expCol2=2, groupLabel2="ip2",
# replicateExpMatrix1=geneExpMatrix1, expColR1=3, replicateLabel1="input1",
# replicateExpMatrix2=geneExpMatrix2, expColR2=3, replicateLabel2="input2",
# method="MATR", outputDir=directory)

