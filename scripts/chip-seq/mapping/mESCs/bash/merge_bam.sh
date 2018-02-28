#!/bin/sh
## ref. Identifying ChIP-seq enrichment using MACS, Nat Protoc., 2012
cd /data/zhoukr/hhl_setd2_m6a/mouse_ESCs_ChIP-seq/Ctrl_D0
samtools merge -f mESCs_ChIP-seq_Ctrl_D0_input.merge.bam mESCs_ChIP-seq_Ctrl_D0_input-rep1.sorted.bam mESCs_ChIP-seq_Ctrl_D0_input-rep2.sorted.bam
samtools merge -f mESCs_ChIP-seq_Ctrl_D0_ChIP.merge.bam mESCs_ChIP-seq_Ctrl_D0_ChIP-rep1.sorted.bam mESCs_ChIP-seq_Ctrl_D0_ChIP-rep2.sorted.bam

cd /data/zhoukr/hhl_setd2_m6a/mouse_ESCs_ChIP-seq/SetD2-KD_D0
samtools merge -f mESCs_ChIP-seq_SetD2-KD_D0_input.merge.bam mESCs_ChIP-seq_SetD2-KD_D0_input-rep1.sorted.bam mESCs_ChIP-seq_SetD2-KD_D0_input-rep2.sorted.bam
samtools merge -f mESCs_ChIP-seq_SetD2-KD_D0_ChIP.merge.bam mESCs_ChIP-seq_SetD2-KD_D0_ChIP-rep1.sorted.bam mESCs_ChIP-seq_SetD2-KD_D0_ChIP-rep2.sorted.bam

