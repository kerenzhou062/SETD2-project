#!/usr/bin/env Rscript
# exomepeak Script 2
# R script
# Define parameters and load library
library("ggplot2")
library("reshape2")

#### m6A and H3K36me3
setwd("/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/mMEFs/xls/m6a-H3K36me3")
# all m6A sites in shCont and shSetD2
shCont<- read.table("./mMEFs_shCont_m6a.txt",sep="\t",header = FALSE)
shSetD2 <- read.table("./mMEFs_shSetD2_m6a.txt",sep="\t",header = FALSE)
wilcoxTest <- wilcox.test(shCont$V1,shSetD2$V1, alternative="greater", paired = FALSE, exact = TRUE, correct = TRUE)
wilcoxPVal = pValFunction(wilcoxTest$p.value)

myData <- read.table("mMEFs_m6A_all_cumulativePlot.txt",sep="\t",header = TRUE)
reshapeData <- melt(myData,id="FE")
levels(reshapeData$variable)[levels(reshapeData$variable)=="shCont"] <- "Control"
levels(reshapeData$variable)[levels(reshapeData$variable)=="shSetD2"] <- paste("shSETD2", " (p ", wilcoxPVal,")", sep = "")
cumulativePlot <- ggplot(data=reshapeData, aes(x=FE, y=value, colour=variable)) + geom_line(size=1.5)+ labs(x="Fold Enrichment(log2)", y = "Cumulative frequency") +
  theme(panel.background = element_rect(fill = "white", colour = "black"), axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 11, face = "bold"), legend.text = element_text(size = 8), legend.title=element_blank(),
    legend.key.width=unit(1,"line"),legend.key.height=unit(1,"line"), legend.background=element_rect(fill = "transparent", colour = "transparent"),
    legend.justification=c(0,1), legend.position=c(0.01,0.99), panel.grid.major.x = element_line(color = "grey80"),
        panel.grid.major.y = element_line(color = "grey80"))
pdf("mMEFs_m6A_all_cumulativePlot.pdf", height=5)
plot(cumulativePlot)
dev.off()

## overlapped m6A sites in shCont and shSetD2 in H3K36me3
shCont<- read.table("./mMEFs_intersect_shCont_m6a+H3K36me3.txt",sep="\t",header = FALSE)
shSetD2 <- read.table("./mMEFs_intersect_shSetD2_m6a+H3K36me3.txt",sep="\t",header = FALSE)
wilcoxTest <- wilcox.test(shCont$V1,shSetD2$V1, alternative="greater", paired = FALSE, exact = TRUE, correct = TRUE)
wilcoxPVal = pValFunction(wilcoxTest$p.value)

myData <- read.table("mMEFs_intersect_m6A_H3K36me3_cumulativePlot.txt",sep="\t",header = TRUE)
reshapeData <- melt(myData,id="FE")
levels(reshapeData$variable)[levels(reshapeData$variable)=="shCont"] <- "Control"
levels(reshapeData$variable)[levels(reshapeData$variable)=="shSetD2"] <- paste("shSETD2", " (p ", wilcoxPVal,")", sep = "")
cumulativePlot <- ggplot(data=reshapeData, aes(x=FE, y=value, colour=variable)) + geom_line(size=1.5)+ labs(x="Fold Enrichment (log2)", y = "Cumulative frequency") +
  theme(panel.background = element_rect(fill = "white", colour = "black"), axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 11, face = "bold"), legend.text = element_text(size = 8), legend.title=element_blank(),
    legend.key.width=unit(1,"line"),legend.key.height=unit(1,"line"), legend.background=element_rect(fill = "transparent", colour = "transparent"),
    legend.justification=c(0,1), legend.position=c(0.01,0.99), panel.grid.major.x = element_line(color = "grey80"),
        panel.grid.major.y = element_line(color = "grey80"))
pdf("mMEFs_intersect_m6A_H3K36me3_cumulativePlot.pdf", height=5)
plot(cumulativePlot)
dev.off()

myData <- read.table("mMEFs_intersect_m6A_H3K36me3.txt",sep="\t",header = TRUE)
wilcoxTest <- wilcox.test(myData$shCont,myData$shSETD2, alternative="greater", paired = TRUE, exact = TRUE, correct = TRUE)
wilcoxPVal = pValFunction(wilcoxTest$p.value)
p1 <- ggplot(myData, aes(x=shCont, y=shSETD2)) +
  geom_point(colour="blue",size=0.5) + xlim(0, 10) + ylim(0,10) + geom_abline(slope=1, intercept=0, color='red', size=1) +
  annotate("text", x = 1, y = 10 , label = paste("p ", wilcoxPVal, sep = ""), fontface="bold", size=4) +
  theme_bw() + labs(x="Fold Enrichment (H3K36me3+) of Control", y = "Fold Enrichment (H3K36me3+) of shSETD2") +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 13, face = "bold"))
pdf("mMEFs_intersect_m6A_H3K36me3_scatter.pdf")
plot(p1)
dev.off()
