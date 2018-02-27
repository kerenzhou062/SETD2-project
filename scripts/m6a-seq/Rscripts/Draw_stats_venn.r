#!/usr/bin/env Rscript
# exomepeak Script 2
# R script
# Define parameters and load library
require(GenomicRanges)
library(ChIPpeakAnno)

setwd("/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/statistics")
bedFile <- read.table('HepG2_shCont_m6a.bed', header=FALSE)
shCont_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

bedFile <- read.table('HepG2_shSetD2_m6a.bed', header=FALSE)
shSetD2_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

bedFile <- read.table('HepG2_shM14_m6a.bed', header=FALSE)
shM14_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

bedFile <- read.table('HepG2_shM3_m6a.bed', header=FALSE)
shM3_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

bedFile <- read.table('HepG2_shWTAP_m6a.bed', header=FALSE)
shWTAP_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

VennForBeds <- makeVennDiagram(Peaks=list(shCont_granges, shSetD2_granges, shM14_granges, shM3_granges, shWTAP_granges),
  NameOfPeaks=c("", "", "", "", ""),
  fill=c("orange","brown","cyan","green","magenta"))
#NameOfPeaks=c("shNS", "shSETD2", "shMETTL14", "shMETTL3", "shWTAP"),
plot(VennForBeds)
dev.off()
