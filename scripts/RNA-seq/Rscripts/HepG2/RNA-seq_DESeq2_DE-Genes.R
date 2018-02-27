#! /usr/bin/env Rscript
# exomepeak Script 2
# R script
# Define parameters and load library
source("https://bioconductor.org/biocLite.R")
library("DESeq2")

setwd("/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/HepG2/input")
condition <- c("treated","treated","treated","untreated","untreated","untreated","untreated")

sampleTable <- read.delim("gene/shCont.vs.shSetD2.sampleTable.txt", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE);
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,directory = "gene",design= ~ condition)
dds <- DESeq(ddsHTSeq)
res <- results(dds)
resOrdered <- res[order(res$pvalue),]
write.table(as.data.frame(resOrdered), file="HepG2_RNA-seq_shCont-shSetD2_htseq-count-gene-DESeq2.FC.txt",quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")

sampleTable <- read.delim("gene/shCont.vs.shM14.sampleTable.txt", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE);
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,directory = "gene",design= ~ condition)
dds <- DESeq(ddsHTSeq)
res <- results(dds)
resOrdered <- res[order(res$pvalue),]
write.table(as.data.frame(resOrdered), file="HepG2_RNA-seq_shCont-shM14_htseq-count-gene-DESeq2.FC.txt",quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")

sampleTable <- read.delim("gene/shCont.vs.shM3.sampleTable.txt", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE);
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,directory = "gene",design= ~ condition)
dds <- DESeq(ddsHTSeq)
res <- results(dds)
resOrdered <- res[order(res$pvalue),]
write.table(as.data.frame(resOrdered), file="HepG2_RNA-seq_shCont-shM3_htseq-count-gene-DESeq2.FC.txt",quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")

sampleTable <- read.delim("gene/shCont.vs.shWTAP.sampleTable.txt", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE);
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,directory = "gene",design= ~ condition)
dds <- DESeq(ddsHTSeq)
res <- results(dds)
resOrdered <- res[order(res$pvalue),]
write.table(as.data.frame(resOrdered), file="HepG2_RNA-seq_shCont-shWTAP_htseq-count-gene-DESeq2.FC.txt",quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")
