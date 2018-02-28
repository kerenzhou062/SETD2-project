#! /usr/bin/env Rscript
# exomepeak Script 2
# R script
# Define parameters and load library
library("exomePeak")
setwd("./SetD2-KO_NA/")
IP_BAM = "mEF_m6A-seq_SetD2-KO_NA_IP.sorted.bam"
INPUT_BAM = "mEF_m6A-seq_SetD2-KO_NA_input.sorted.bam"
GTF_ANNO = "/data/sunwj/genomeInfo/mouse/mm10/annotation/gencode.vM9.annotation.gtf"
# comparison
exomepeak(GENE_ANNO_GTF=GTF_ANNO, IP_BAM=IP_BAM, INPUT_BAM=INPUT_BAM, EXPERIMENT_NAME="SetD2-KO_NA_m6A")
