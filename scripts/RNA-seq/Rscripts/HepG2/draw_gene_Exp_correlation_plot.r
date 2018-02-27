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

setwd("/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/HepG2/RSEM-correlation")

myData <- read.table("HepG2_RNA-seq_RSEM_shSetD2vsshCont_log2_input.txt",sep="\t",header = TRUE)
log2FPKM <- myData[2:3]
shCont <- log2FPKM$shCont
shSETD2 <- log2FPKM$shSETD2
cor_test_shCont <- cor.test(shCont,shSETD2,method="pearson",exact=TRUE)
cor_test_shCont_cor = sprintf("%.2f",cor_test_shCont$estimate)
cor_test_shCont_pval = pValFunction(cor_test_shCont$p.value)

p <- ggplot(log2FPKM, aes(x=shSETD2, y=shCont)) +
  geom_point(colour="blue", size = 1) + geom_smooth(method=lm,se=FALSE,colour="red",size = 2) +
  annotate("text", x = 8.5, y = 1, label = paste("r=", cor_test_shCont_cor, sep = ""),fontface="bold") +
  annotate("text", x = 9.0, y = 0.5, label = paste("p", cor_test_shCont_pval, sep = ""),fontface="bold") +
  theme_bw() + labs(x="shSETD2 (log2FPKM)", y = "shCont (log2FPKM)") +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 13, face = "bold"))

pdf("HepG2_RNA-seq_RSEM_shSetD2vsshCont_log2_input_correlation.pdf")
plot(p)
dev.off()


myData <- read.table("HepG2_RNA-seq_RSEM_shSetD2vsshM14_log2_input.txt",sep="\t",header = TRUE)
log2FPKM <- myData[2:3]
shMETTL14 <- log2FPKM$shMETTL14
shSETD2 <- log2FPKM$shSETD2
cor_test_shMETTL14 <- cor.test(shMETTL14,shSETD2,method="pearson",exact=TRUE)
cor_test_shMETTL14_cor = sprintf("%.2f",cor_test_shMETTL14$estimate)
cor_test_shMETTL14_pval = pValFunction(cor_test_shMETTL14$p.value)

p <- ggplot(log2FPKM, aes(x=shSETD2, y=shMETTL14)) +
  geom_point(colour="blue", size = 1) + geom_smooth(method=lm,se=FALSE,colour="red",size = 2) +
  annotate("text", x = 8.5, y = 1, label = paste("r=", cor_test_shMETTL14_cor, sep = ""),fontface="bold") +
  annotate("text", x = 9.0, y = 0.5, label = paste("p", cor_test_shMETTL14_pval, sep = ""),fontface="bold") +
  theme_bw() + labs(x="shSETD2 (log2FPKM)", y = "shMETTL14 (log2FPKM)") +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 13, face = "bold"))

pdf("HepG2_RNA-seq_RSEM_shSetD2vsshM14_log2_input_correlation.pdf")
plot(p)
dev.off()


myData <- read.table("HepG2_RNA-seq_RSEM_shSetD2vsshM3_log2_input.txt",sep="\t",header = TRUE)
log2FPKM <- myData[2:3]
shMETTL3 <- log2FPKM$shMETTL3
shSETD2 <- log2FPKM$shSETD2
cor_test_shMETTL3 <- cor.test(shMETTL3,shSETD2,method="pearson",exact=TRUE)
cor_test_shMETTL3_cor = sprintf("%.2f",cor_test_shMETTL3$estimate)
cor_test_shMETTL3_pval = pValFunction(cor_test_shMETTL3$p.value)

p <- ggplot(log2FPKM, aes(x=shSETD2, y=shMETTL3)) +
  geom_point(colour="blue", size = 1) + geom_smooth(method=lm,se=FALSE,colour="red",size = 2) +
  annotate("text", x = 8.5, y = 1, label = paste("r=", cor_test_shMETTL3_cor, sep = ""),fontface="bold") +
  annotate("text", x = 9.0, y = 0.5, label = paste("p", cor_test_shMETTL3_pval, sep = ""),fontface="bold") +
  theme_bw() + labs(x="shSETD2 (log2FPKM)", y = "shMETTL3 (log2FPKM)") +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 13, face = "bold"))

pdf("HepG2_RNA-seq_RSEM_shSetD2vsshM3_log2_input_correlation.pdf")
plot(p)
dev.off()


myData <- read.table("HepG2_RNA-seq_RSEM_shSetD2vsshWTAP_log2_input.txt",sep="\t",header = TRUE)
log2FPKM <- myData[2:3]
shWTAP <- log2FPKM$shWTAP
shSETD2 <- log2FPKM$shSETD2
cor_test_shWTAP <- cor.test(shWTAP,shSETD2,method="pearson",exact=TRUE)
cor_test_shWTAP_cor = sprintf("%.2f",cor_test_shWTAP$estimate)
cor_test_shWTAP_pval = pValFunction(cor_test_shWTAP$p.value)

p <- ggplot(log2FPKM, aes(x=shSETD2, y=shWTAP)) +
  geom_point(colour="blue", size = 1) + geom_smooth(method=lm,se=FALSE,colour="red",size = 2) +
  annotate("text", x = 8.5, y = 1, label = paste("r=", cor_test_shWTAP_cor, sep = ""),fontface="bold") +
  annotate("text", x = 9.0, y = 0.5, label = paste("p", cor_test_shWTAP_pval, sep = ""),fontface="bold") +
  theme_bw() + labs(x="shSETD2 (log2FPKM)", y = "shWTAP (log2FPKM)") +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 13, face = "bold"))

pdf("HepG2_RNA-seq_RSEM_shSetD2vsshWTAP_log2_input_correlation.pdf")
plot(p)
dev.off()

