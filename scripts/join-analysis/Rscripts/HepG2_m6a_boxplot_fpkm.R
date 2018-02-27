#!/usr/bin/env Rscript
# exomepeak Script 2
# R script
# Define parameters and load library
library("ggplot2")
library("reshape2")
#library("devtools")
#install_github("const-ae/ggsignif")
library("ggsignif")
setwd("/data/zhoukr/hhl_setd2_m6a/analysis/joint-analysis/HepG2/fpkm")

#peak fold enrichment
myData <- read.table(file="HepG2_join_exp-histone-m6A_foldEnrichment_log2.txt",header=TRUE)
myData$xorder <- factor(myData$expClass, levels = c("low","medium","high"))
ylim1 <- boxplot.stats(myData$H3K36me3)$stats[c(1, 5)]
tipLen = c(0.01, 0.01)
vjust = 0.3
boxplot <- ggplot(myData, aes(x = xorder, y = H3K36me3, fill=xorder)) +
  geom_boxplot(outlier.shape=NA) +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  labs(y="H3K36me3 Modification Levels", x="Gene Expression Level",
    size = 13, face = "bold") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        panel.background = element_rect(fill = "#F2F2F2",
                                colour = "#F2F2F2",
                                size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                colour = "white"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                colour = "white")) +
  geom_signif(comparisons = list(c("low", "medium")),
               map_signif_level = TRUE, test = "wilcox.test", textsize=5,
               y_position=ylim1[2]*1.1, tip_length = tipLen, vjust=vjust) +
  geom_signif(comparisons = list(c("medium", "high")),
               map_signif_level = TRUE, test = "wilcox.test", textsize=5,
               y_position=ylim1[2]*1.14, tip_length = tipLen, vjust=vjust) +
  geom_signif(comparisons = list(c("low", "high")),
               map_signif_level = TRUE, test = "wilcox.test", textsize=5,
               y_position=ylim1[2]*1.18, tip_length = tipLen, vjust=vjust) +
  coord_cartesian(ylim = c(ylim1[1]*0.9, ylim1[2]*1.2))

pdf("HepG2_join_exp-histone-m6A_foldEnrichment_log2_H3K36me3.pdf", height=8)
plot(boxplot)
dev.off()

myData$xorder <- factor(myData$expClass, levels = c("low","medium","high"))
ylim1 <- boxplot.stats(myData$m6A)$stats[c(1, 5)]
tipLen = c(0.01, 0.01)
vjust = 0.3
boxplot <- ggplot(myData, aes(x = xorder, y = m6A, fill=xorder)) +
  geom_boxplot(outlier.shape=NA) +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  labs(y="m6A Modification Levels", x="Gene Expression Level",
    size = 13, face = "bold") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        panel.background = element_rect(fill = "#F2F2F2",
                                colour = "#F2F2F2",
                                size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                colour = "white"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                colour = "white")) +
  geom_signif(comparisons = list(c("low", "medium")),
               map_signif_level = TRUE, test = "wilcox.test", textsize=5,
               y_position=ylim1[2]*1.1, tip_length = tipLen, vjust=vjust) +
  geom_signif(comparisons = list(c("medium", "high")),
               map_signif_level = TRUE, test = "wilcox.test", textsize=5,
               y_position=ylim1[2]*1.14, tip_length = tipLen, vjust=vjust) +
  geom_signif(comparisons = list(c("low", "high")),
               map_signif_level = TRUE, test = "wilcox.test", textsize=5,
               y_position=ylim1[2]*1.18, tip_length = tipLen, vjust=vjust) +
  coord_cartesian(ylim = c(ylim1[1]*0.9, ylim1[2]*1.2))

pdf("HepG2_join_exp-histone-m6A_foldEnrichment_log2_m6A.pdf", height=8)
plot(boxplot)
dev.off()

#gene fpkm
myData <- read.table(file="HepG2_join_exp-histone-m6A_fpkm_log2.txt",header=TRUE)
myData$xorder <- factor(myData$expClass, levels = c("low","medium","high"))
ylim1 <- boxplot.stats(myData$H3K36me3)$stats[c(1, 5)]
tipLen = c(0.01, 0.01)
vjust = 0.3
boxplot <- ggplot(myData, aes(x = xorder, y = H3K36me3, fill=xorder)) +
  geom_boxplot(outlier.shape=NA) +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  labs(y="H3K36me3 Modification Levels", x="Gene Expression Level",
    size = 13, face = "bold") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        panel.background = element_rect(fill = "#F2F2F2",
                                colour = "#F2F2F2",
                                size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                colour = "white"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                colour = "white")) +
  geom_signif(comparisons = list(c("low", "medium")),
               map_signif_level = TRUE, test = "wilcox.test", textsize=5,
               y_position=ylim1[2]*1.1, tip_length = tipLen, vjust=vjust) +
  geom_signif(comparisons = list(c("medium", "high")),
               map_signif_level = TRUE, test = "wilcox.test", textsize=5,
               y_position=ylim1[2]*1.14, tip_length = tipLen, vjust=vjust) +
  geom_signif(comparisons = list(c("low", "high")),
               map_signif_level = TRUE, test = "wilcox.test", textsize=5,
               y_position=ylim1[2]*1.18, tip_length = tipLen, vjust=vjust) +
  coord_cartesian(ylim = c(ylim1[1]*0.9, ylim1[2]*1.2))

pdf("HepG2_join_exp-histone-m6A_fpkm_log2_H3K36me3.pdf", height=8)
plot(boxplot)
dev.off()

myData$xorder <- factor(myData$expClass, levels = c("low","medium","high"))
ylim1 <- boxplot.stats(myData$m6A)$stats[c(1, 5)]
tipLen = c(0.01, 0.01)
vjust = 0.3
boxplot <- ggplot(myData, aes(x = xorder, y = m6A, fill=xorder)) +
  geom_boxplot(outlier.shape=NA) +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  labs(y="m6A Modification Levels", x="Gene Expression Level",
    size = 13, face = "bold") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        panel.background = element_rect(fill = "#F2F2F2",
                                colour = "#F2F2F2",
                                size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                colour = "white"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                colour = "white")) +
  geom_signif(comparisons = list(c("low", "medium")),
               map_signif_level = TRUE, test = "wilcox.test", textsize=5,
               y_position=ylim1[2]*1.2, tip_length = tipLen, vjust=vjust) +
  geom_signif(comparisons = list(c("medium", "high")),
               map_signif_level = TRUE, test = "wilcox.test", textsize=5,
               y_position=ylim1[2]*1.24, tip_length = tipLen, vjust=vjust) +
  geom_signif(comparisons = list(c("low", "high")),
               map_signif_level = TRUE, test = "wilcox.test", textsize=5,
               y_position=ylim1[2]*1.28, tip_length = tipLen, vjust=vjust) +
  coord_cartesian(ylim = c(ylim1[1]*0.9, ylim1[2]*1.3))

pdf("HepG2_join_exp-histone-m6A_fpkm_log2_m6A.pdf", height=8)
plot(boxplot)
dev.off()
