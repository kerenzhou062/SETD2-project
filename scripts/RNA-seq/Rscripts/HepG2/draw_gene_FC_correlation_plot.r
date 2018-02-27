#!/usr/bin/env Rscript
# exomepeak Script 2
# R script
# Define parameters and load library
library("ggplot2")
library("corrplot")

pValFunction <- function(pval) {
  if (pval == 0) {
    return ("<2.22e-308")
  }else if (pval == 1) {
    return ("=1")
  }else{
    exactVal = sprintf("%.2e", pval)
    exactVal = paste("=", exactVal, sep = "")
    return (exactVal)
  }
}

setwd("/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/HepG2/input")

myData <- read.table("HepG2_RNA-seq_shSETD2vsshMETTL14_gene_FC.txt",sep="\t",header = TRUE)
log2FC <- myData[3:4]
shMETTL14 <- log2FC$shMETTL14
shSETD2 <- log2FC$shSETD2
cor_test_shMETTL14 <- cor.test(shMETTL14,shSETD2,method="pearson",exact=TRUE)
cor_test_shMETTL14_cor = sprintf("%.2f",cor_test_shMETTL14$estimate)
cor_test_shMETTL14_pval = pValFunction(cor_test_shMETTL14$p.value)

p <- ggplot(log2FC, aes(x=shSETD2, y=shMETTL14)) +
  geom_point(colour="blue", size = 2) + geom_smooth(method=lm,se=FALSE,colour="red",size = 2) +
  annotate("text", x = 5.5, y = -4.7, label = paste("r=", cor_test_shMETTL14_cor, sep = ""),fontface="bold") +
  annotate("text", x = 6.0, y = -5.0, label = paste("p", cor_test_shMETTL14_pval, sep = ""),fontface="bold") +
  theme_bw() + labs(x="shSETD2 (log2FC)", y = "shMETTL14 (log2FC)") +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 13, face = "bold"))

pdf("HepG2_RNA-seq_FC_correlation_shMETTL14vs.shSETD2.pdf")
plot(p)
dev.off()

myData <- read.table("HepG2_RNA-seq_shSETD2vsshMETTL3_gene_FC.txt",sep="\t",header = TRUE)
log2FC <- myData[3:4]
shMETTL3 <- log2FC$shMETTL3
shSETD2 <- log2FC$shSETD2
cor_test_shMETTL3 <- cor.test(shMETTL3,shSETD2,method="pearson",exact=TRUE)
cor_test_shMETTL3_cor = sprintf("%.2f",cor_test_shMETTL3$estimate)
cor_test_shMETTL3_pval = pValFunction(cor_test_shMETTL3$p.value)

p <- ggplot(log2FC, aes(x=shSETD2, y=shMETTL3)) +
  geom_point(colour="blue", size = 2) + geom_smooth(method=lm,se=FALSE,colour="red",size = 2) +
  annotate("text", x = 5.5, y = -4.7, label = paste("r=", cor_test_shMETTL3_cor, sep = ""),fontface="bold") +
  annotate("text", x = 6.0, y = -5.0, label = paste("p", cor_test_shMETTL3_pval, sep = ""),fontface="bold") +
  theme_bw() + labs(x="shSETD2 (log2FC)", y = "shMETTL3 (log2FC)") +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 13, face = "bold"))

pdf("HepG2_RNA-seq_FC_correlation_shMETTL3vs.shSETD2.pdf")
plot(p)
dev.off()

myData <- read.table("HepG2_RNA-seq_shSETD2vsshWTAP_gene_FC.txt",sep="\t",header = TRUE)
log2FC <- myData[3:4]
shWTAP <- log2FC$shWTAP
shSETD2 <- log2FC$shSETD2
cor_test_shWTAP <- cor.test(shWTAP,shSETD2,method="pearson",exact=TRUE)
cor_test_shWTAP_cor = sprintf("%.2f",cor_test_shWTAP$estimate)
cor_test_shWTAP_pval = pValFunction(cor_test_shWTAP$p.value)

p <- ggplot(log2FC, aes(x=shSETD2, y=shWTAP)) +
  geom_point(colour="blue", size = 2) + geom_smooth(method=lm,se=FALSE,colour="red",size = 2) +
  annotate("text", x = 5.5, y = -4.7, label = paste("r=", cor_test_shWTAP_cor, sep = ""),fontface="bold") +
  annotate("text", x = 5.9, y = -5.0, label = paste("p", cor_test_shWTAP_pval, sep = ""),fontface="bold") +
  theme_bw() + labs(x="shSETD2 (log2FC)", y = "shWTAP (log2FC)") +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 13, face = "bold"))

pdf("HepG2_RNA-seq_FC_correlation_shWTAPvs.shSETD2.pdf")
plot(p)
dev.off()

