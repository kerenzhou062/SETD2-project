#!/bin/sh

fastqc HepG2_ChIP-seq_shM14_IP_rep1.fastq HepG2_ChIP-seq_shSetD2_input_rep2.fastq HepG2_ChIP-seq_shCont_input_rep2.fastq HepG2_ChIP-seq_shCont_IP_rep2.fastq HepG2_ChIP-seq_shM3_IP_rep1.fastq HepG2_ChIP-seq_shSetD2_IP_rep2.fastq HepG2_ChIP-seq_shCont_IP_rep1.fastq HepG2_ChIP-seq_shWTAP_input_rep1.fastq HepG2_ChIP-seq_shM14_input_rep1.fastq HepG2_ChIP-seq_shM3_input_rep1.fastq HepG2_ChIP-seq_shSetD2_input_rep1.fastq HepG2_ChIP-seq_shWTAP_IP_rep1.fastq HepG2_ChIP-seq_shSetD2_IP_rep1.fastq HepG2_ChIP-seq_shCont_input_rep1.fastq -o 0_fastqc_output > fastqc.run.log 2>&1
