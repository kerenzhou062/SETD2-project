#!/usr/bin/env Rscript
# exomepeak Script 2
# R script
# Define parameters and load library
library("ggplot2")
library("reshape2")
#library("devtools")
#install_github("const-ae/ggsignif")
library("ggsignif")
setwd("/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/mouse/xls/overlappedPeak")
myData <- read.table(file="mouse_all_samples_m6a_peak_FE_log2.txt",header=TRUE)
myLog2Data <- myData[15:ncol(myData)]
reshapeData <- melt(myLog2Data,id=NULL)
colnames(reshapeData) <- c("sample","value")
levels(reshapeData$sample)[levels(reshapeData$sample)=="mESCs.shCont.D0"] <- "mESCs shCont D0"
levels(reshapeData$sample)[levels(reshapeData$sample)=="mESCs.shCont.D6"] <- "mESCs shCont D6"
levels(reshapeData$sample)[levels(reshapeData$sample)=="mESCs.shSETD2.D0"] <- "mESCs shSETD2 D0"
levels(reshapeData$sample)[levels(reshapeData$sample)=="mESCs.shSETD2.D6"] <- "mESCs shSETD2 D6"
levels(reshapeData$sample)[levels(reshapeData$sample)=="mMEFs.shCont"] <- "mEFs shCont"
levels(reshapeData$sample)[levels(reshapeData$sample)=="mMEFs.shSETD2"] <- "mEFs shSETD2"

ylim1 <- boxplot.stats(reshapeData$value)$stats[c(1, 5)]
tipLen = c(0.01, 0.01)
vjust = 0.65

boxplot <- ggplot(reshapeData, aes(x = sample, y = value, fill=sample)) +
  geom_boxplot(outlier.shape=NA) +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  labs(y="Fold Enrichment of m6A Peaks (log2)", size = 13, face = "bold") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        panel.background = element_rect(fill = "#F2F2F2",
                                colour = "#F2F2F2",
                                size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                colour = "white"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                colour = "white")) +
  geom_signif(comparisons = list(c("mESCs shCont D0", "mESCs shCont D6")),
               map_signif_level = TRUE, test = "wilcox.test", textsize=5,
               y_position=ylim1[2]*1.03, tip_length = tipLen, vjust=vjust) +
  geom_signif(comparisons = list(c("mESCs shCont D0", "mESCs shSETD2 D0")),
               map_signif_level = TRUE, test = "wilcox.test", textsize=3,
               y_position=ylim1[2]*1.09, tip_length = tipLen, vjust=vjust) +
  geom_signif(comparisons = list(c("mESCs shCont D6", "mESCs shSETD2 D6")),
               map_signif_level = TRUE, test = "wilcox.test", textsize=5,
               y_position=ylim1[2]*1.13, tip_length = tipLen, vjust=vjust) +
  geom_signif(comparisons = list(c("mESCs shCont D0", "mESCs shSETD2 D6")),
               map_signif_level = TRUE, test = "wilcox.test", textsize=5,
               y_position=ylim1[2]*1.16, tip_length = tipLen, vjust=vjust) +
  geom_signif(comparisons = list(c("mESCs shSETD2 D6", "mEFs shCont")),
               map_signif_level = TRUE, test = "wilcox.test", textsize=3,
               y_position=ylim1[2]*1.08, tip_length = tipLen, vjust=vjust) +
  geom_signif(comparisons = list(c("mEFs shCont", "mEFs shSETD2")),
               map_signif_level = TRUE, test = "wilcox.test", textsize=3,
               y_position=ylim1[2]*1.03, tip_length = tipLen, vjust=vjust) +

  coord_cartesian(ylim = c(ylim1[1]*0.9, ylim1[2]*1.2))

pdf("mouse_all_sample_m6a_foldEnrichment_log2_boxplot.pdf", height=8)
plot(boxplot)
dev.off()
