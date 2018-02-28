#!/bin/sh
## ref. Identifying ChIP-seq enrichment using MACS, Nat Protoc., 2012

samtools merge /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shCont/input/Hela_ChIP-seq_shCont_input.fastq.merge.bam /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shCont/input_rep1/Hela_ChIP-seq_shCont_input_rep1.fastq.sorted.bam /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shCont/input_rep2/Hela_ChIP-seq_shCont_input_rep2.fastq.sorted.bam
samtools merge /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shCont/IP/Hela_ChIP-seq_shCont_IP.fastq.merge.bam /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shCont/IP_rep1/Hela_ChIP-seq_shCont_IP_rep1.fastq.sorted.bam /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shCont/IP_rep2/Hela_ChIP-seq_shCont_IP_rep2.fastq.sorted.bam
samtools merge /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shSetD2/input/Hela_ChIP-seq_shSetD2_input.fastq.merge.bam /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shSetD2/input_rep1/Hela_ChIP-seq_shSetD2_input_rep1.fastq.sorted.bam /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shSetD2/input_rep2/Hela_ChIP-seq_shSetD2_input_rep2.fastq.sorted.bam
samtools merge /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shSetD2/IP/Hela_ChIP-seq_shSetD2_IP.fastq.merge.bam /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shSetD2/IP_rep1/Hela_ChIP-seq_shSetD2_IP_rep1.fastq.sorted.bam /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shSetD2/IP_rep2/Hela_ChIP-seq_shSetD2_IP_rep2.fastq.sorted.bam

