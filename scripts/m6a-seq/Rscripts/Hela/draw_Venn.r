#!/usr/bin/env Rscript
# exomepeak Script 2
# R script
# Define parameters and load library
library("VennDiagram")
library("grid")

setwd("/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/Hela/venn")
shSetD2_gene <- read.table("Hela_shSetD2_FC_gene.txt",sep="\t",header = FALSE, stringsAsFactors=TRUE)
shM14_gene <- read.table("Hela_shM14_FC_gene.txt",sep="\t",header = FALSE, stringsAsFactors=TRUE)
shM3_gene <- read.table("Hela_shM3_FC_gene.txt",sep="\t",header = FALSE, stringsAsFactors=TRUE)
shWTAP_gene <- read.table("Hela_shWTAP_FC_gene.txt",sep="\t",header = FALSE, stringsAsFactors=TRUE)

shSetD2_gene <- shSetD2_gene[,1]
shM14_gene <- shM14_gene[,1]
shM3_gene <- shM3_gene[,1]
shWTAP_gene <- shWTAP_gene[,1]

futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
venn.diagram(x= list(shSETD2 = shSetD2_gene,shMETTL14 = shM14_gene,shMETTL3 = shM3_gene,shWTAP = shWTAP_gene), filename = "Hela_venn.tiff",resolution = 300,
  imagetype="tiff", col="transparent", fill=c("cornflowerblue","seagreen3","darkorange1","darkorchid1"), cat.col=c("cornflowerblue","seagreen3","darkorange1","darkorchid1"),
  fontface = "bold", alpha = 0.5, cex=1.2, cat.cex=1.2)

