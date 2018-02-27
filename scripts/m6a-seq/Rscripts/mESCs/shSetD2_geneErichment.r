#!/usr/bin/env Rscript
# exomepeak Script 2
# R script
# Define parameters and load library
library("clusterProfiler")
library("org.Mm.eg.db")

setwd("/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/mESCs/xls/GOSET")
myData <- read.table("mESCs_D0_-F0.5_shSetD2_hypomethylation_geneID.txt",sep="\t",header = TRUE)
geneTransform <- bitr(myData$ENSEMBL, fromType="ENSEMBL", toType=c("ENTREZID"), OrgDb="org.Mm.eg.db")
geneENTREZID=as.character(geneTransform[,2])
geneTransform <- bitr(geneENTREZID, fromType="ENTREZID", toType=c("ENSEMBL"), OrgDb="org.Mm.eg.db")
gene=as.character(geneTransform[,2])

GO_result <- enrichGO(gene         = gene,
                 OrgDb         = org.Mm.eg.db,
                 keytype       = 'ENSEMBL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.001,
                 qvalueCutoff  = 0.005,
                 readable      = TRUE)

write.table(GO_result, file="mESCs_D0_-F0.5_shSetD2_hypomethylation_gene_GO_result.xls", quote=FALSE, sep="\t")

myData <- read.table("mESCs_D0_exomePeak_shSetD2_hypomethylation_geneID.txt",sep="\t",header = TRUE)
geneTransform <- bitr(myData$ENSEMBL, fromType="ENSEMBL", toType=c("ENTREZID"), OrgDb="org.Mm.eg.db")
geneENTREZID=as.character(geneTransform[,2])
geneTransform <- bitr(geneENTREZID, fromType="ENTREZID", toType=c("ENSEMBL"), OrgDb="org.Mm.eg.db")
gene=as.character(geneTransform[,2])

GO_result <- enrichGO(gene         = gene,
                 OrgDb         = org.Mm.eg.db,
                 keytype       = 'ENSEMBL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.001,
                 qvalueCutoff  = 0.005,
                 readable      = TRUE)

write.table(GO_result, file="mESCs_D0_exomePeak_shSetD2_hypomethylation_gene_GO_result.xls", quote=FALSE, sep="\t")


myData <- read.table("mESCs_D6_-F0.5_shSetD2_hypomethylation_geneID.txt",sep="\t",header = TRUE)
geneTransform <- bitr(myData$ENSEMBL, fromType="ENSEMBL", toType=c("ENTREZID"), OrgDb="org.Mm.eg.db")
geneENTREZID=as.character(geneTransform[,2])
geneTransform <- bitr(geneENTREZID, fromType="ENTREZID", toType=c("ENSEMBL"), OrgDb="org.Mm.eg.db")
gene=as.character(geneTransform[,2])

GO_result <- enrichGO(gene         = gene,
                 OrgDb         = org.Mm.eg.db,
                 keytype       = 'ENSEMBL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)

write.table(GO_result, file="mESCs_D6_-F0.5_shSetD2_hypomethylation_gene_GO_result.xls", quote=FALSE, sep="\t")

myData <- read.table("mESCs_D6_exomePeak_shSetD2_hypomethylation_geneID.txt",sep="\t",header = TRUE)
geneTransform <- bitr(myData$ENSEMBL, fromType="ENSEMBL", toType=c("ENTREZID"), OrgDb="org.Mm.eg.db")
geneENTREZID=as.character(geneTransform[,2])
geneTransform <- bitr(geneENTREZID, fromType="ENTREZID", toType=c("ENSEMBL"), OrgDb="org.Mm.eg.db")
gene=as.character(geneTransform[,2])

GO_result <- enrichGO(gene         = gene,
                 OrgDb         = org.Mm.eg.db,
                 keytype       = 'ENSEMBL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)

write.table(GO_result, file="mESCs_D6_exomePeak_shSetD2_hypomethylation_gene_GO_result.xls", quote=FALSE, sep="\t")

myData <- read.table("mESCs_shSetD2vsshCont_D6vsD0_-F0.5_shSetD2_hypomethylation_geneID.txt",sep="\t",header = TRUE)
geneTransform <- bitr(myData$ENSEMBL, fromType="ENSEMBL", toType=c("ENTREZID"), OrgDb="org.Mm.eg.db")
geneENTREZID=as.character(geneTransform[,2])
geneTransform <- bitr(geneENTREZID, fromType="ENTREZID", toType=c("ENSEMBL"), OrgDb="org.Mm.eg.db")
gene=as.character(geneTransform[,2])

GO_result <- enrichGO(gene         = gene,
                 OrgDb         = org.Mm.eg.db,
                 keytype       = 'ENSEMBL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)

write.table(GO_result, file="mESCs_shSetD2vsshCont_D6vsD0_-F0.5_shSetD2_hypomethylation_gene_GO_result.xls", quote=FALSE, sep="\t")
