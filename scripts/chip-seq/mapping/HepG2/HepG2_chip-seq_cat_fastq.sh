#!/bin/sh

cat HepG2_ChIP-seq_shCont_IP_rep1.fastq HepG2_ChIP-seq_shCont_IP_rep2.fastq > HepG2_ChIP-seq_shCont_IP.fastq
cat HepG2_ChIP-seq_shCont_input_rep1.fastq HepG2_ChIP-seq_shCont_input_rep2.fastq > HepG2_ChIP-seq_shCont_input.fastq

cat HepG2_ChIP-seq_shSetD2_IP_rep1.fastq HepG2_ChIP-seq_shSetD2_IP_rep2.fastq > HepG2_ChIP-seq_shSetD2_IP.fastq
cat HepG2_ChIP-seq_shSetD2_input_rep1.fastq HepG2_ChIP-seq_shSetD2_input_rep2.fastq > HepG2_ChIP-seq_shSetD2_input.fastq
