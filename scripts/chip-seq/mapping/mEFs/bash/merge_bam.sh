#!/bin/sh
## ref. Identifying ChIP-seq enrichment using MACS, Nat Protoc., 2012
cd /data/zhoukr/hhl_setd2_m6a/mouse_MEFs_ChIP-seq/SetD2-KO_NA
samtools merge -f mMEFs_ChIP-seq_SetD2-KO_NA_ChIP.merge.bam mMEFs_ChIP-seq_SetD2-KO_NA_ChIP-rep1.sorted.bam mMEFs_ChIP-seq_SetD2-KO_NA_ChIP-rep2.sorted.bam
samtools merge -f mMEFs_ChIP-seq_SetD2-KO_NA_input.merge.bam mMEFs_ChIP-seq_SetD2-KO_NA_input-rep1.sorted.bam mMEFs_ChIP-seq_SetD2-KO_NA_input-rep2.sorted.bam

cd /data/zhoukr/hhl_setd2_m6a/mouse_MEFs_ChIP-seq/SetD2-WT_NA
samtools merge -f mMEFs_ChIP-seq_SetD2-WT_NA_ChIP.merge.bam mMEFs_ChIP-seq_SetD2-WT_NA_ChIP-rep1.sorted.bam mMEFs_ChIP-seq_SetD2-WT_NA_ChIP-rep2.sorted.bam
samtools merge -f mMEFs_ChIP-seq_SetD2-WT_NA_input.merge.bam mMEFs_ChIP-seq_SetD2-WT_NA_input-rep1.sorted.bam mMEFs_ChIP-seq_SetD2-WT_NA_input-rep2.sorted.bam

