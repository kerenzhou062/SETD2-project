#!/usr/bin/env Rscript
# exomepeak Script 2
# R script
# Define parameters and load library
library("ggplot2")
library("pheatmap")
library("data.table")
setwd("/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/Hela/xls/overlappedPeak")

myData <- read.table("Hela_m6a_peak_FC_log2.txt",sep="\t",header = TRUE)
pheatmapData <- myData[14:17]
setnames(pheatmapData, old=c("shSetD2", "shM14","shM3"), new=c("shSETD2", "shMETTL14", "shMETTL3"))
condition <- data.frame("KD" = c("shSETD2", "shMETTL14", "shMETTL3", "shWTAP"))
rownames(condition) <- colnames(pheatmapData)

pheatmap(pheatmapData, annotation = condition, cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = F,
	scale="none", show_rownames = F, legend_breaks = c(-5, -2, 2, 5, 6),
	main = "", legend_labels = c("-5", "-2", "2", "5", "log2FC\n"),
	legend = TRUE, color = colorRampPalette(c("blue", "white", "red"))(1000),
	filename = "Hela_m6a_peak_FC_log2_heatmap.pdf")
