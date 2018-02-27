#!/usr/bin/env Rscript
# exomepeak Script 2
# R script
# Define parameters and load library
require(GenomicRanges)
library(ChIPpeakAnno)

setwd("/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/Hela/statistics")

## shCont replicates
bedFile <- read.table('Hela_shCont_m6a_rep1.bed', header=FALSE)
rep1_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

bedFile <- read.table('Hela_shCont_m6a_rep2.bed', header=FALSE)
rep2_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

bedFile <- read.table('Hela_shCont_m6a_rep3.bed', header=FALSE)
rep3_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

pdf("Hela_shCont_replicates_venn.pdf")
makeVennDiagram(Peaks=list(rep1_granges, rep2_granges, rep3_granges),
  NameOfPeaks=c("Rep1", "Rep2", "Rep3"), fill=c("orange","cyan","magenta"), margin=0.05, fontface = "bold", alpha = 0.5, cex=1.2, cat.cex=1.2)
dev.off()

## shSetD2 replicates
bedFile <- read.table('Hela_shSetD2_m6a_rep1.bed', header=FALSE)
rep1_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

bedFile <- read.table('Hela_shSetD2_m6a_rep2.bed', header=FALSE)
rep2_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

bedFile <- read.table('Hela_shSetD2_m6a_rep3.bed', header=FALSE)
rep3_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

pdf("Hela_shSetD2_replicates_venn.pdf")
makeVennDiagram(Peaks=list(rep1_granges, rep2_granges, rep3_granges),
  NameOfPeaks=c("Rep1", "Rep2", "Rep3"), fill=c("orange","cyan","magenta"), margin=0.05, fontface = "bold", alpha = 0.5, cex=1.2, cat.cex=1.2)
dev.off()

## shM14 replicates
bedFile <- read.table('Hela_shM14_m6a_rep1.bed', header=FALSE)
rep1_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

bedFile <- read.table('Hela_shM14_m6a_rep2.bed', header=FALSE)
rep2_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

bedFile <- read.table('Hela_shM14_m6a_rep3.bed', header=FALSE)
rep3_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

pdf("Hela_shM14_replicates_venn.pdf")
makeVennDiagram(Peaks=list(rep1_granges, rep2_granges, rep3_granges),
  NameOfPeaks=c("Rep1", "Rep2", "Rep3"), fill=c("orange","cyan","magenta"), margin=0.05, fontface = "bold", alpha = 0.5, cex=1.2, cat.cex=1.2)
dev.off()

## shM3 replicates
bedFile <- read.table('Hela_shM3_m6a_rep1.bed', header=FALSE)
rep1_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

bedFile <- read.table('Hela_shM3_m6a_rep2.bed', header=FALSE)
rep2_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

bedFile <- read.table('Hela_shM3_m6a_rep3.bed', header=FALSE)
rep3_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

pdf("Hela_shM3_replicates_venn.pdf")
makeVennDiagram(Peaks=list(rep1_granges, rep2_granges, rep3_granges),
  NameOfPeaks=c("Rep1", "Rep2", "Rep3"), fill=c("orange","cyan","magenta"), margin=0.05, fontface = "bold", alpha = 0.5, cex=1.2, cat.cex=1.2)
dev.off()

## shWTAP replicates
bedFile <- read.table('Hela_shWTAP_m6a_rep1.bed', header=FALSE)
rep1_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

bedFile <- read.table('Hela_shWTAP_m6a_rep2.bed', header=FALSE)
rep2_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

bedFile <- read.table('Hela_shWTAP_m6a_rep3.bed', header=FALSE)
rep3_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

pdf("Hela_shWTAP_replicates_venn.pdf")
makeVennDiagram(Peaks=list(rep1_granges, rep2_granges, rep3_granges),
  NameOfPeaks=c("Rep1", "Rep2", "Rep3"), fill=c("orange","cyan","magenta"), margin=0.05, fontface = "bold", alpha = 0.5, cex=1.2, cat.cex=1.2)
dev.off()

## shTargets, rep1-2
bedFile <- read.table('Hela_shCont_m6a_rep1-2.bed', header=FALSE)
shCont_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

bedFile <- read.table('Hela_shSetD2_m6a_rep1-2.bed', header=FALSE)
shSetD2_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

bedFile <- read.table('Hela_shM14_m6a_rep1-2.bed', header=FALSE)
shM14_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

bedFile <- read.table('Hela_shM3_m6a_rep1-2.bed', header=FALSE)
shM3_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

bedFile <- read.table('Hela_shWTAP_m6a_rep1-2.bed', header=FALSE)
shWTAP_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

pdf("Hela_m6a_shTargets_rep1-2_venn.pdf", width=8)
makeVennDiagram(Peaks=list(shCont_granges, shSetD2_granges, shM14_granges, shM3_granges, shWTAP_granges),
  NameOfPeaks=c("shNS", "shSETD2", "shMETTL14", "shMETTL3", "shWTAP"), margin=0.05,
  fill=c("orange","brown","cyan","green","magenta"), fontface = "bold", alpha = 0.5, cex=1, cat.cex=1)
dev.off()

## shTargets, rep1-3
bedFile <- read.table('Hela_shCont_m6a_rep1-3.bed', header=FALSE)
shCont_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

bedFile <- read.table('Hela_shSetD2_m6a_rep1-3.bed', header=FALSE)
shSetD2_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

bedFile <- read.table('Hela_shM14_m6a_rep1-3.bed', header=FALSE)
shM14_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

bedFile <- read.table('Hela_shM3_m6a_rep1-3.bed', header=FALSE)
shM3_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

bedFile <- read.table('Hela_shWTAP_m6a_rep1-3.bed', header=FALSE)
shWTAP_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

pdf("Hela_m6a_shTargets_rep1-3_venn.pdf", width=8)
makeVennDiagram(Peaks=list(shCont_granges, shSetD2_granges, shM14_granges, shM3_granges, shWTAP_granges),
  NameOfPeaks=c("shNS", "shSETD2", "shMETTL14", "shMETTL3", "shWTAP"), margin=0.05,
  fill=c("orange","brown","cyan","green","magenta"), fontface = "bold", alpha = 0.5, cex=1, cat.cex=1)
dev.off()

## shTargets, rep2-3
bedFile <- read.table('Hela_shCont_m6a_rep2-3.bed', header=FALSE)
shCont_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

bedFile <- read.table('Hela_shSetD2_m6a_rep2-3.bed', header=FALSE)
shSetD2_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

bedFile <- read.table('Hela_shM14_m6a_rep2-3.bed', header=FALSE)
shM14_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

bedFile <- read.table('Hela_shM3_m6a_rep2-3.bed', header=FALSE)
shM3_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

bedFile <- read.table('Hela_shWTAP_m6a_rep2-3.bed', header=FALSE)
shWTAP_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

pdf("Hela_m6a_shTargets_rep2-3_venn.pdf", width=8)
makeVennDiagram(Peaks=list(shCont_granges, shSetD2_granges, shM14_granges, shM3_granges, shWTAP_granges),
  NameOfPeaks=c("shNS", "shSETD2", "shMETTL14", "shMETTL3", "shWTAP"), margin=0.05,
  fill=c("orange","brown","cyan","green","magenta"), fontface = "bold", alpha = 0.5, cex=1, cat.cex=1)
dev.off()


## shTargets, rep1-2-3
bedFile <- read.table('Hela_shCont_m6a_rep1-2-3.bed', header=FALSE)
shCont_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

bedFile <- read.table('Hela_shSetD2_m6a_rep1-2-3.bed', header=FALSE)
shSetD2_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

bedFile <- read.table('Hela_shM14_m6a_rep1-2-3.bed', header=FALSE)
shM14_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

bedFile <- read.table('Hela_shM3_m6a_rep1-2-3.bed', header=FALSE)
shM3_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

bedFile <- read.table('Hela_shWTAP_m6a_rep1-2-3.bed', header=FALSE)
shWTAP_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

pdf("Hela_m6a_shTargets_rep1-2-3_venn.pdf", width=8)
makeVennDiagram(Peaks=list(shCont_granges, shSetD2_granges, shM14_granges, shM3_granges, shWTAP_granges),
  NameOfPeaks=c("shNS", "shSETD2", "shMETTL14", "shMETTL3", "shWTAP"), margin=0.05,
  fill=c("orange","brown","cyan","green","magenta"), fontface = "bold", alpha = 0.5, cex=1, cat.cex=1)
dev.off()

## shTargets, rep1
bedFile <- read.table('Hela_shCont_m6a_rep1.bed', header=FALSE)
shCont_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

bedFile <- read.table('Hela_shSetD2_m6a_rep1.bed', header=FALSE)
shSetD2_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

bedFile <- read.table('Hela_shM14_m6a_rep1.bed', header=FALSE)
shM14_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

bedFile <- read.table('Hela_shM3_m6a_rep1.bed', header=FALSE)
shM3_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

bedFile <- read.table('Hela_shWTAP_m6a_rep1.bed', header=FALSE)
shWTAP_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

pdf("Hela_m6a_shTargets_rep1_venn.pdf", width=8)
makeVennDiagram(Peaks=list(shCont_granges, shSetD2_granges, shM14_granges, shM3_granges, shWTAP_granges),
  NameOfPeaks=c("shNS", "shSETD2", "shMETTL14", "shMETTL3", "shWTAP"), margin=0.05,
  fill=c("orange","brown","cyan","green","magenta"), fontface = "bold", alpha = 0.5, cex=1, cat.cex=1)
dev.off()

## shTargets, rep2
bedFile <- read.table('Hela_shCont_m6a_rep2.bed', header=FALSE)
shCont_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

bedFile <- read.table('Hela_shSetD2_m6a_rep2.bed', header=FALSE)
shSetD2_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

bedFile <- read.table('Hela_shM14_m6a_rep2.bed', header=FALSE)
shM14_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

bedFile <- read.table('Hela_shM3_m6a_rep2.bed', header=FALSE)
shM3_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

bedFile <- read.table('Hela_shWTAP_m6a_rep2.bed', header=FALSE)
shWTAP_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

pdf("Hela_m6a_shTargets_rep2_venn.pdf", width=8)
makeVennDiagram(Peaks=list(shCont_granges, shSetD2_granges, shM14_granges, shM3_granges, shWTAP_granges),
  NameOfPeaks=c("shNS", "shSETD2", "shMETTL14", "shMETTL3", "shWTAP"), margin=0.05,
  fill=c("orange","brown","cyan","green","magenta"), fontface = "bold", alpha = 0.5, cex=1, cat.cex=1)
dev.off()

## shTargets, rep3
bedFile <- read.table('Hela_shCont_m6a_rep3.bed', header=FALSE)
shCont_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

bedFile <- read.table('Hela_shSetD2_m6a_rep3.bed', header=FALSE)
shSetD2_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

bedFile <- read.table('Hela_shM14_m6a_rep3.bed', header=FALSE)
shM14_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

bedFile <- read.table('Hela_shM3_m6a_rep3.bed', header=FALSE)
shM3_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

bedFile <- read.table('Hela_shWTAP_m6a_rep3.bed', header=FALSE)
shWTAP_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

pdf("Hela_m6a_shTargets_rep3_venn.pdf", width=8)
makeVennDiagram(Peaks=list(shCont_granges, shSetD2_granges, shM14_granges, shM3_granges, shWTAP_granges),
  NameOfPeaks=c("shNS", "shSETD2", "shMETTL14", "shMETTL3", "shWTAP"), margin=0.05,
  fill=c("orange","brown","cyan","green","magenta"), fontface = "bold", alpha = 0.5, cex=1, cat.cex=1)
dev.off()

