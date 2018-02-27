#!/usr/bin/env Rscript
# exomepeak Script 2
# R script
# Define parameters and load library
library("ggplot2")
library("reshape2")

setwd("/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/Hela/xls/overlappedPeak")

M14TargetFC <- read.table("./target/Hela_shM14_FC_log2.txt",sep="\t",header = FALSE)
M3TargetFC <- read.table("./target/Hela_shM3_FC_log2.txt",sep="\t",header = FALSE)
WTAPTargetFC <- read.table("./target/Hela_shWTAP_FC_log2.txt",sep="\t",header = FALSE)
shareTargetFC <- read.table("./target/Hela_shareTargets_FC_log2.txt",sep="\t",header = FALSE)
nonTargetFC <- read.table("./target/Hela_nonTargets_FC_log2.txt",sep="\t",header = FALSE)

pValFunction <- function(pval) {
  if (pval == 0) {
    return ("< 2.22e-308")
  }else if (pval == 1) {
    return ("= 1")
  }else{
    exactVal = sprintf("%.2e", pval)
    exactVal = paste("= ", exactVal, sep = "")
    return (exactVal)
  }
}

myData <- read.table("Hela_cumulativePlot.txt",sep="\t",header = TRUE)

wilcoxTest <- wilcox.test(M14TargetFC$V1,nonTargetFC$V1, alternative="less", paired = FALSE, exact = TRUE, correct = TRUE)
wilcoxShM14PVal = pValFunction(wilcoxTest$p.value)
wilcoxTest <- wilcox.test(M3TargetFC$V1,nonTargetFC$V1, alternative="less", paired = FALSE, exact = TRUE, correct = TRUE)
wilcoxshM3PVal = pValFunction(wilcoxTest$p.value)
wilcoxTest <- wilcox.test(WTAPTargetFC$V1,nonTargetFC$V1, alternative="less", paired = FALSE, exact = TRUE, correct = TRUE)
wilcoxshWTAPPVal = pValFunction(wilcoxTest$p.value)
wilcoxTest <- wilcox.test(shareTargetFC$V1,nonTargetFC$V1, alternative="less", paired = FALSE, exact = TRUE, correct = TRUE)
wilcoxsharePVal = pValFunction(wilcoxTest$p.value)

#wilcoxTest <- wilcox.test(myData$shM14,myData$nonTarget, alternative="greater", paired = FALSE, exact = TRUE, correct = TRUE)
#wilcoxShM14PVal = pValFunction(wilcoxTest$p.value)
#wilcoxTest <- wilcox.test(myData$shM3,myData$nonTarget, alternative="greater", paired = FALSE, exact = TRUE, correct = TRUE)
#wilcoxshM3PVal = pValFunction(wilcoxTest$p.value)
#wilcoxTest <- wilcox.test(myData$shWTAP,myData$nonTarget, alternative="greater", paired = FALSE, exact = TRUE, correct = TRUE)
#wilcoxshWTAPPVal = pValFunction(wilcoxTest$p.value)
#wilcoxTest <- wilcox.test(myData$share,myData$nonTarget, alternative="greater", paired = FALSE, exact = TRUE, correct = TRUE)
#wilcoxsharePVal = pValFunction(wilcoxTest$p.value)

reshapeData <- melt(myData,id="FC")

levels(reshapeData$variable)[levels(reshapeData$variable)=="shM14"] <- paste("METTL14 target", " (p ", wilcoxShM14PVal,")", sep = "")
levels(reshapeData$variable)[levels(reshapeData$variable)=="shM3"] <- paste("METTL3 target", " (p ", wilcoxshM3PVal,")", sep = "")
levels(reshapeData$variable)[levels(reshapeData$variable)=="shWTAP"] <- paste("WTAP target", " (p ", wilcoxshWTAPPVal,")", sep = "")
levels(reshapeData$variable)[levels(reshapeData$variable)=="share"] <- paste("share target", " (p ", wilcoxsharePVal,")", sep = "")
levels(reshapeData$variable)[levels(reshapeData$variable)=="nonTarget"] <- "non target"

overlappedPeakCum <- ggplot(data=reshapeData, aes(x=FC, y=value, colour=variable)) + geom_line(size=1.5)+ labs(x="Fold Enrichment(log2)", y = "Cumulative frequency") +
  theme(panel.background = element_rect(fill = "white", colour = "black"), axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 11, face = "bold"), legend.text = element_text(size = 8), legend.title=element_blank(),
    legend.key.width=unit(1,"line"),legend.key.height=unit(1,"line"), legend.background=element_rect(fill = "transparent", colour = "transparent"),
    legend.justification=c(0,1), legend.position=c(0.01,0.99), panel.grid.major.x = element_line(color = "grey80"),
        panel.grid.major.y = element_line(color = "grey80"))

pdf("Hela_overlappedCumulativePlot.pdf", height=5)
plot(overlappedPeakCum)
dev.off()

#### m6A and H3K36me3
setwd("/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/Hela/xls/m6a-H3K36me3")
# all m6A sites in shCont and shSetD2
shCont<- read.table("./Hela_shCont_m6a.txt",sep="\t",header = FALSE)
shSetD2 <- read.table("./Hela_shSetD2_m6a.txt",sep="\t",header = FALSE)
wilcoxTest <- wilcox.test(shCont$V1,shSetD2$V1, alternative="greater", paired = FALSE, exact = TRUE, correct = TRUE)
wilcoxPVal = pValFunction(wilcoxTest$p.value)

myData <- read.table("Hela_m6A_all_cumulativePlot.txt",sep="\t",header = TRUE)
reshapeData <- melt(myData,id="FE")
levels(reshapeData$variable)[levels(reshapeData$variable)=="shCont"] <- "Control"
levels(reshapeData$variable)[levels(reshapeData$variable)=="shSetD2"] <- paste("shSETD2", " (p ", wilcoxPVal,")", sep = "")
cumulativePlot <- ggplot(data=reshapeData, aes(x=FE, y=value, colour=variable)) + geom_line(size=1.5)+ labs(x="Fold Enrichment(log2)", y = "Cumulative frequency") +
  theme(panel.background = element_rect(fill = "white", colour = "black"), axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 11, face = "bold"), legend.text = element_text(size = 8), legend.title=element_blank(),
    legend.key.width=unit(1,"line"),legend.key.height=unit(1,"line"), legend.background=element_rect(fill = "transparent", colour = "transparent"),
    legend.justification=c(0,1), legend.position=c(0.01,0.99), panel.grid.major.x = element_line(color = "grey80"),
        panel.grid.major.y = element_line(color = "grey80"))
pdf("Hela_m6A_all_cumulativePlot.pdf", height=5)
plot(cumulativePlot)
dev.off()

## overlapped m6A sites in shCont and shSetD2 in H3K36me3
shCont<- read.table("./Hela_intersect_shCont_m6a+H3K36me3.txt",sep="\t",header = FALSE)
shSetD2 <- read.table("./Hela_intersect_shSetD2_m6a+H3K36me3.txt",sep="\t",header = FALSE)
wilcoxTest <- wilcox.test(shCont$V1,shSetD2$V1, alternative="greater", paired = FALSE, exact = TRUE, correct = TRUE)
wilcoxPVal = pValFunction(wilcoxTest$p.value)

myData <- read.table("Hela_intersect_m6A_H3K36me3_cumulativePlot.txt",sep="\t",header = TRUE)
reshapeData <- melt(myData,id="FE")
levels(reshapeData$variable)[levels(reshapeData$variable)=="shCont"] <- "Control"
levels(reshapeData$variable)[levels(reshapeData$variable)=="shSetD2"] <- paste("shSETD2", " (p ", wilcoxPVal,")", sep = "")
cumulativePlot <- ggplot(data=reshapeData, aes(x=FE, y=value, colour=variable)) + geom_line(size=1.5)+ labs(x="Fold Enrichment (log2)", y = "Cumulative frequency") +
  theme(panel.background = element_rect(fill = "white", colour = "black"), axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 11, face = "bold"), legend.text = element_text(size = 8), legend.title=element_blank(),
    legend.key.width=unit(1,"line"),legend.key.height=unit(1,"line"), legend.background=element_rect(fill = "transparent", colour = "transparent"),
    legend.justification=c(0,1), legend.position=c(0.01,0.99), panel.grid.major.x = element_line(color = "grey80"),
        panel.grid.major.y = element_line(color = "grey80"))
pdf("Hela_intersect_m6A_H3K36me3_cumulativePlot.pdf", height=5)
plot(cumulativePlot)
dev.off()

myData <- read.table("Hela_intersect_m6A_H3K36me3.txt",sep="\t",header = TRUE)
wilcoxTest <- wilcox.test(myData$shCont,myData$shSETD2, alternative="greater", paired = TRUE, exact = TRUE, correct = TRUE)
wilcoxPVal = pValFunction(wilcoxTest$p.value)
p1 <- ggplot(myData, aes(x=shCont, y=shSETD2)) +
  geom_point(colour="blue",size=0.8) + xlim(1, 8) + ylim(1,8) + geom_abline(slope=1, intercept=0, color='red', size=1) +
  annotate("text", x = 1, y = max(myData$shSETD2) , label = paste("p ", wilcoxPVal, sep = ""), fontface="bold", size=4) +
  theme_bw() + labs(x="Fold Enrichment (H3K36me3+) of Control", y = "Fold Enrichment (H3K36me3+) of shSETD2") +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 13, face = "bold"))
pdf("Hela_intersect_m6A_H3K36me3_scatter.pdf")
plot(p1)
dev.off()

## overlapped m6A sites in shCont and shSetD2
shCont<- read.table("./Hela_intersect_shCont_m6a+shSetD2.txt",sep="\t",header = FALSE)
shSetD2 <- read.table("./Hela_intersect_shSetD2_m6a+shCont.txt",sep="\t",header = FALSE)
wilcoxTest <- wilcox.test(shCont$V1,shSetD2$V1, alternative="greater", paired = FALSE, exact = TRUE, correct = TRUE)
wilcoxPVal = pValFunction(wilcoxTest$p.value)

myData <- read.table("Hela_intersect_m6A_cumulativePlot.txt",sep="\t",header = TRUE)
reshapeData <- melt(myData,id="FE")
levels(reshapeData$variable)[levels(reshapeData$variable)=="shCont"] <- "Control"
levels(reshapeData$variable)[levels(reshapeData$variable)=="shSetD2"] <- paste("shSETD2", " (p ", wilcoxPVal,")", sep = "")
cumulativePlot <- ggplot(data=reshapeData, aes(x=FE, y=value, colour=variable)) + geom_line(size=1.5)+ labs(x="Fold Enrichment (log2)", y = "Cumulative frequency") +
  theme(panel.background = element_rect(fill = "white", colour = "black"), axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 11, face = "bold"), legend.text = element_text(size = 8), legend.title=element_blank(),
    legend.key.width=unit(1,"line"),legend.key.height=unit(1,"line"), legend.background=element_rect(fill = "transparent", colour = "transparent"),
    legend.justification=c(0,1), legend.position=c(0.01,0.99), panel.grid.major.x = element_line(color = "grey80"),
        panel.grid.major.y = element_line(color = "grey80"))
pdf("Hela_intersect_m6A_cumulativePlot.pdf", height=5)
plot(cumulativePlot)
dev.off()

myData <- read.table("Hela_intersect_m6A.txt",sep="\t",header = TRUE)
wilcoxTest <- wilcox.test(myData$shCont,myData$shSETD2, alternative="greater", paired = TRUE, exact = TRUE, correct = TRUE)
wilcoxPVal = pValFunction(wilcoxTest$p.value)
p1 <- ggplot(myData, aes(x=shCont, y=shSETD2)) +
  geom_point(colour="blue",size=0.8) + xlim(1, 8) + ylim(1,8) + geom_abline(slope=1, intercept=0, color='red', size=1) +
  annotate("text", x = 1.5, y = max(myData$shSETD2)-1 , label = paste("p ", wilcoxPVal, sep = ""), fontface="bold", size=4) +
  theme_bw() + labs(x="Fold Enrichment of Control", y = "Fold Enrichment of shSETD2") +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 13, face = "bold"))
pdf("Hela_intersect_m6A_scatter.pdf")
plot(p1)
dev.off()
