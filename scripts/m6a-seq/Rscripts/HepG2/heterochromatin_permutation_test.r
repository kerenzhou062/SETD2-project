#!/usr/bin/env Rscript
# exomepeak Script 2
# R script
# Define parameters and load library
setwd("/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/heterochomatin")
library("rtracklayer")
library("GenomicRanges")
library("regioneR")
genome <- filterChromosomes(getGenome("hg19"))
m6a <- import("HepG2_shCont_m6a.bed")
H3K36me3 <- import("HepG2_macs_shCont_peaks.bed")
H3K9me3 <- import("narrowPeak_H3K9me3.bed")
m6a_H3K36me3 <- overlapPermTest(A=m6a, B=H3K36me3, ntimes=1000, genome=genome, verbose=FALSE, alternative='greater')

sink("m6a_H3K36me3_permutation_test.txt")
print(summary(m6a_H3K36me3))
sink()
pdf("m6a_H3K36me3_permutation_test.pdf")
plot(m6a_H3K36me3)
dev.off()

m6a_H3K9me3 <- overlapPermTest(A=m6a, B=H3K9me3, ntimes=1000, genome=genome, verbose=FALSE, alternative='greater')
sink("m6a_H3K9me3_permutation_test.txt")
print(summary(m6a_H3K9me3))
sink()
pdf("m6a_H3K9me3_permutation_test.pdf")
plot(m6a_H3K9me3)
dev.off()
