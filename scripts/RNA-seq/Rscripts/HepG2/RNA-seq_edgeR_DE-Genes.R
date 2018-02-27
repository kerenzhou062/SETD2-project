#! /usr/bin/env Rscript
# exomepeak Script 2
# R script
# Define parameters and load library
source("https://bioconductor.org/biocLite.R")
library("edgeR")

setwd("/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/HepG2/input")
rawdata <- read.delim("gene/HepG2_RNA-seq_shCont-shSetD2_htseq-count.txt",check.names=FALSE, stringsAsFactors=FALSE);
y <- DGEList(counts=rawdata[,2:ncol(rawdata)], genes=rawdata[,1])
Tissue <- factor(c("N","N","N","T","T","T"))
data.frame(Sample=colnames(y),Tissue)
design <- model.matrix(~Tissue)
rownames(design) <- colnames(y)
y <- estimateDisp(y, design, robust=TRUE)
#y$common.dispersion
#plotBCV(y)
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
#topTags(lrt)
out <- topTags(lrt, n=Inf, adjust.method="BH")
write.table(out, file="HepG2_RNA-seq_shCont-shSetD2_htseq-count-gene-edgeR.FC.txt", quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")

rawdata <- read.delim("gene/HepG2_RNA-seq_shCont-shM14_htseq-count.txt",check.names=FALSE, stringsAsFactors=FALSE);
y <- DGEList(counts=rawdata[,2:ncol(rawdata)], genes=rawdata[,1])
Tissue <- factor(c("N","N","N","T","T","T"))
data.frame(Sample=colnames(y),Tissue)
design <- model.matrix(~Tissue)
rownames(design) <- colnames(y)
y <- estimateDisp(y, design, robust=TRUE)
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
out <- topTags(lrt, n=Inf, adjust.method="BH")
write.table(out, file="HepG2_RNA-seq_shCont-shM14_htseq-count-gene-edgeR.FC.txt", quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")

rawdata <- read.delim("gene/HepG2_RNA-seq_shCont-shM3_htseq-count.txt",check.names=FALSE, stringsAsFactors=FALSE);
y <- DGEList(counts=rawdata[,2:ncol(rawdata)], genes=rawdata[,1])
Tissue <- factor(c("N","N","N","T","T","T"))
data.frame(Sample=colnames(y),Tissue)
design <- model.matrix(~Tissue)
rownames(design) <- colnames(y)
y <- estimateDisp(y, design, robust=TRUE)
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
out <- topTags(lrt, n=Inf, adjust.method="BH")
write.table(out, file="HepG2_RNA-seq_shCont-shM3_htseq-count-gene-edgeR.FC.txt", quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")

rawdata <- read.delim("gene/HepG2_RNA-seq_shCont-shWTAP_htseq-count.txt",check.names=FALSE, stringsAsFactors=FALSE);
y <- DGEList(counts=rawdata[,2:ncol(rawdata)], genes=rawdata[,1])
Tissue <- factor(c("N","N","N","T","T","T"))
data.frame(Sample=colnames(y),Tissue)
design <- model.matrix(~Tissue)
rownames(design) <- colnames(y)
y <- estimateDisp(y, design, robust=TRUE)
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
out <- topTags(lrt, n=Inf, adjust.method="BH")
write.table(out, file="HepG2_RNA-seq_shCont-shWTAP_htseq-count-gene-edgeR.FC.txt", quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")

rawdata <- read.delim("gene/HepG2_RNA-seq_shCont-shSetD2_htseq-count.txt",check.names=FALSE, stringsAsFactors=FALSE);
y <- DGEList(counts=rawdata[,2:ncol(rawdata)], genes=rawdata[,1])
Tissue <- factor(c("N","N","N","T","T","T"))
data.frame(Sample=colnames(y),Tissue)
design <- model.matrix(~Tissue)
rownames(design) <- colnames(y)
y <- estimateDisp(y, design, robust=TRUE)
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
out <- topTags(lrt, n=Inf, adjust.method="BH")
write.table(out, file="HepG2_RNA-seq_shCont-shSetD2_htseq-count-gene-edgeR.FC.txt", quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")
