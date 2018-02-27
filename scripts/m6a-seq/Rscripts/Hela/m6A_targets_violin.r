#!/usr/bin/env Rscript
# exomepeak Script 2
# R script
# Define parameters and load library
library("ggplot2")
library("reshape2")

setwd("/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/Hela/xls/allPeak/")

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

wilcoxTest <- wilcox.test(M14TargetFC$V1,nonTargetFC$V1, alternative="less", paired = FALSE, exact = TRUE, correct = TRUE)
wilcoxShM14PVal = pValFunction(wilcoxTest$p.value)
wilcoxTest <- wilcox.test(M3TargetFC$V1,nonTargetFC$V1, alternative="less", paired = FALSE, exact = TRUE, correct = TRUE)
wilcoxshM3PVal = pValFunction(wilcoxTest$p.value)
wilcoxTest <- wilcox.test(WTAPTargetFC$V1,nonTargetFC$V1, alternative="less", paired = FALSE, exact = TRUE, correct = TRUE)
wilcoxshWTAPPVal = pValFunction(wilcoxTest$p.value)
wilcoxTest <- wilcox.test(shareTargetFC$V1,nonTargetFC$V1, alternative="less", paired = FALSE, exact = TRUE, correct = TRUE)
wilcoxsharePVal = pValFunction(wilcoxTest$p.value)

targetFC <- read.table("./target/Hela_m6A_target_violin_Plot.txt",sep="\t",header = TRUE)
targetFC$Type <- factor(targetFC$Type, levels = c("Share","METTL14","METTL3", "WTAP", "Non-Target"))
ylim1 <- boxplot.stats(targetFC$FC)$stats[c(1, 5)]
levels(targetFC$Type)[levels(targetFC$Type)=="METTL14"] <- paste("METTL14", " (p ", wilcoxShM14PVal,")", sep = "")
levels(targetFC$Type)[levels(targetFC$Type)=="METTL3"] <- paste("METTL3", " (p ", wilcoxshM3PVal,")", sep = "")
levels(targetFC$Type)[levels(targetFC$Type)=="WTAP"] <- paste("WTAP", " (p ", wilcoxshWTAPPVal,")", sep = "")
levels(targetFC$Type)[levels(targetFC$Type)=="Share"] <- paste("Share", " (p ", wilcoxsharePVal,")", sep = "")
levels(targetFC$Type)[levels(targetFC$Type)=="Non-Target"] <- "Non-target"
violin <- ggplot(targetFC, aes(factor(Type), FC)) +
  geom_violin(aes(fill = Type),stat = "ydensity",
    position = "dodge", adjust=1/2) +
  labs(y="Fold Change of m6A peaks (log2(shSETD2/shCont))", x="",
    size = 13, face = "bold") +
  theme(axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    legend.title=element_blank(),
    legend.key.width=unit(1,"line"),legend.key.height=unit(1,"line"),
    legend.background=element_rect(fill = "transparent", colour = "transparent"),
    legend.justification=c(0,1), legend.position=c(0.01,1),
    panel.background = element_rect(fill = "#F2F2F2",
        colour = "#F2F2F2",
        size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
        colour = "white"),
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
        colour = "white")) +
  coord_cartesian(ylim = c(ylim1[1]*1.1, ylim1[2]*2.2))

pdf("Hela_m6A_target_violin_Plot.pdf", height=6)
plot(violin)
dev.off()
