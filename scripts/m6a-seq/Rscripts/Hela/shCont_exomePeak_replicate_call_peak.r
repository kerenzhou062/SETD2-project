#!/usr/bin/env Rscript
# exomepeak Script 2
# R script
# Define parameters and load library
library("exomePeak")
setwd("./shCont/")

INPUT1_BAM = "Hela_m6A-seq_shCont_input_rep1.fastq.sorted.bam"
INPUT2_BAM = "Hela_m6A-seq_shCont_input_rep2.fastq.sorted.bam"
INPUT3_BAM = "Hela_m6A-seq_shCont_input_rep3.fastq.sorted.bam"
IP1_BAM = "Hela_m6A-seq_shCont_IP_rep1.fastq.sorted.bam"
IP2_BAM = "Hela_m6A-seq_shCont_IP_rep2.fastq.sorted.bam"
IP3_BAM = "Hela_m6A-seq_shCont_IP_rep3.fastq.sorted.bam"

GTF_ANNO = "/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.gtf"
# comparison
#exomepeak(GENE_ANNO_GTF=GTF_ANNO, IP_BAM=c(IP1_BAM, IP2_BAM), INPUT_BAM=c(INPUT1_BAM, INPUT2_BAM), EXPERIMENT_NAME="shCont_m6A_rep1-2")
#exomepeak(GENE_ANNO_GTF=GTF_ANNO, IP_BAM=c(IP1_BAM, IP3_BAM), INPUT_BAM=c(INPUT1_BAM, INPUT3_BAM), EXPERIMENT_NAME="shCont_m6A_rep1-3")
exomepeak(GENE_ANNO_GTF=GTF_ANNO, IP_BAM=c(IP2_BAM, IP3_BAM), INPUT_BAM=c(INPUT2_BAM, INPUT3_BAM), EXPERIMENT_NAME="shCont_m6A_rep2-3")
#exomepeak(GENE_ANNO_GTF=GTF_ANNO, IP_BAM=c(IP1_BAM, IP2_BAM, IP3_BAM), INPUT_BAM=c(INPUT1_BAM, INPUT2_BAM, INPUT3_BAM), EXPERIMENT_NAME="shCont_m6A_rep1-2-3")
#exomepeak(GENE_ANNO_GTF=GTF_ANNO, IP_BAM=c(IP1_BAM), INPUT_BAM=c(INPUT1_BAM), EXPERIMENT_NAME="shCont_m6A_rep1")
#exomepeak(GENE_ANNO_GTF=GTF_ANNO, IP_BAM=c(IP2_BAM), INPUT_BAM=c(INPUT2_BAM), EXPERIMENT_NAME="shCont_m6A_rep2")
#exomepeak(GENE_ANNO_GTF=GTF_ANNO, IP_BAM=c(IP3_BAM), INPUT_BAM=c(INPUT3_BAM), EXPERIMENT_NAME="shCont_m6A_rep3")
