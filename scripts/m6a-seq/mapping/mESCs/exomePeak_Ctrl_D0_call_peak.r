#! /usr/bin/env Rscript
# exomepeak Script 2
# R script
# Define parameters and load library
library("exomePeak")
setwd("./Ctrl_D0/")
IP_BAM = "mES_m6A-seq_Ctrl_D0_IP.sorted.bam"
INPUT_BAM = "mES_m6A-seq_Ctrl_D0_input.sorted.bam"
GTF_ANNO = "/data/sunwj/genomeInfo/mouse/mm10/annotation/gencode.vM9.annotation.gtf"
# comparison
exomepeak(GENE_ANNO_GTF=GTF_ANNO, IP_BAM=IP_BAM, INPUT_BAM=INPUT_BAM, EXPERIMENT_NAME="Ctrl_D0_m6A")
