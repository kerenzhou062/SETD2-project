#!/usr/bin/env Rscript
# exomepeak Script 2
# R script
# Define parameters and load library
library("ggplot2")
library("reshape2")
#library("devtools")
#install_github("const-ae/ggsignif")
library("ggsignif")
setwd("/data/zhoukr/hhl_setd2_m6a/analysis/joint-analysis/Hela/DEGseq")

#gene fpkm
myData <- read.table(file="Hela_join_exp-m6A_DEGseq_log2.txt",header=TRUE)
myData$Type <- factor(myData$expClass, levels = c("low","medium","high"))
reshapeData <- melt(myData[5:ncol(myData)])
ylim1 <- boxplot.stats(reshapeData$value)$stats[c(1, 5)]
tipLen = c(0.01, 0.01)
vjust = 0.3

boxplot <- ggplot(reshapeData, aes(factor(variable), value)) +
  geom_boxplot(aes(fill = Type), outlier.shape=NA) +
  labs(y="m6A Modification Levels", x="",
    size = 13, face = "bold") +
  theme(axis.title.x=element_blank(),
        panel.background = element_rect(fill = "#F2F2F2",
                                colour = "#F2F2F2",
                                size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                colour = "white"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                colour = "white")) +
  coord_cartesian(ylim = c(ylim1[1]*0.8, ylim1[2]*1.4))

pdf("Hela_join_exp-m6a-DEGseq_log2.pdf", height=8)
plot(boxplot)
dev.off()

myData <- read.table(file="Hela_join_exp-histone_DEGseq_log2.txt",header=TRUE)
myData$Type <- factor(myData$expClass, levels = c("low","medium","high"))
reshapeData <- melt(myData[4:ncol(myData)])
ylim1 <- boxplot.stats(reshapeData$value)$stats[c(1, 5)]
tipLen = c(0.01, 0.01)
vjust = 0.3

boxplot <- ggplot(reshapeData, aes(factor(variable), value)) +
  geom_boxplot(aes(fill = Type), outlier.shape=NA) +
  labs(y="H3K36me3 Modification Levels", x="",
    size = 13, face = "bold") +
  theme(axis.title.x=element_blank(),
        panel.background = element_rect(fill = "#F2F2F2",
                                colour = "#F2F2F2",
                                size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                colour = "white"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                colour = "white")) +
  coord_cartesian(ylim = c(ylim1[1]*0.8, ylim1[2]*1.4))

pdf("Hela_join_exp-histone-DEGseq_log2.pdf", height=8)
plot(boxplot)
dev.off()
