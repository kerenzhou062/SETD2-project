#!/usr/bin/env Rscript
# exomepeak Script 2
# R script
# Define parameters and load library
library("DEGseq")
directory="/data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/mouse/DEGseq"
setwd(directory)
featureCount="mouse_ChIP-seq_featureCount.txt"
colCount <- ncol(read.table(file=featureCount, nrows=1))
countMatrix <- readGeneExp(file=featureCount, geneCol=1, valCol=c(7:colCount))

DEGexp(geneExpMatrix1=countMatrix, geneCol1=1, expCol1=c(2,3), groupLabel1="ip",
    geneExpMatrix2=countMatrix, geneCol2=1, expCol2=c(4,5), groupLabel2="input",
    method="MARS", outputDir=paste(directory, "mESCsShContD0", sep="/"))

DEGexp(geneExpMatrix1=countMatrix, geneCol1=1, expCol1=c(6,7), groupLabel1="ip",
    geneExpMatrix2=countMatrix, geneCol2=1, expCol2=c(8,9), groupLabel2="input",
    method="MARS", outputDir=paste(directory, "mESCsShSetD2D0", sep="/"))

DEGexp(geneExpMatrix1=countMatrix, geneCol1=1, expCol1=c(10,11), groupLabel1="ip",
    geneExpMatrix2=countMatrix, geneCol2=1, expCol2=c(12,13), groupLabel2="input",
    method="MARS", outputDir=paste(directory, "mEFsShCont", sep="/"))

DEGexp(geneExpMatrix1=countMatrix, geneCol1=1, expCol1=c(14,15), groupLabel1="ip",
    geneExpMatrix2=countMatrix, geneCol2=1, expCol2=c(16,17), groupLabel2="input",
    method="MARS", outputDir=paste(directory, "mEFsShSetD2", sep="/"))

#DEGexp(geneExpMatrix1=geneExpMatrix1, expCol1=2, groupLabel1="ip1",
# geneExpMatrix2=geneExpMatrix2, expCol2=2, groupLabel2="ip2",
# replicateExpMatrix1=geneExpMatrix1, expColR1=3, replicateLabel1="input1",
# replicateExpMatrix2=geneExpMatrix2, expColR2=3, replicateLabel2="input2",
# method="MATR", outputDir=directory)

