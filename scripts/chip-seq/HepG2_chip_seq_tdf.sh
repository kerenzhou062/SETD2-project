#!/bin/sh
bamFolder=/data/zhoukr/hhl_setd2_m6a/HepG2_ChIP-seq
prefix=HepG2_ChIP-seq
gsize=hg19.chrom.sizes

rm -rf /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/tdf
mkdir -p /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/tdf
cd /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/tdf
cp /data/zhoukr/reference/genome/genome_size/number/human_hg19_chrsize.txt $gsize


nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/shCont/input_rep1/${prefix}_shCont_input_rep1.fastq.sorted.bam  ${prefix}_shCont_input_rep1.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/shCont/input_rep2/${prefix}_shCont_input_rep2.fastq.sorted.bam  ${prefix}_shCont_input_rep2.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/shCont/input/${prefix}_shCont_input.fastq.merge.bam  ${prefix}_shCont_input_merge.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/shCont/IP_rep1/${prefix}_shCont_IP_rep1.fastq.sorted.bam  ${prefix}_shCont_IP_rep1.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/shCont/IP_rep2/${prefix}_shCont_IP_rep2.fastq.sorted.bam  ${prefix}_shCont_IP_rep2.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/shCont/IP/${prefix}_shCont_IP.fastq.merge.bam  ${prefix}_shCont_IP_merge.cov.tdf $gsize > /dev/null 2>&1 &


nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/shSetD2/input_rep1/${prefix}_shSetD2_input_rep1.fastq.sorted.bam  ${prefix}_shSetD2_input_rep1.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/shSetD2/input_rep2/${prefix}_shSetD2_input_rep2.fastq.sorted.bam  ${prefix}_shSetD2_input_rep2.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/shSetD2/input/${prefix}_shSetD2_input.fastq.merge.bam  ${prefix}_shSetD2_input_merge.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/shSetD2/IP_rep1/${prefix}_shSetD2_IP_rep1.fastq.sorted.bam  ${prefix}_shSetD2_IP_rep1.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/shSetD2/IP_rep2/${prefix}_shSetD2_IP_rep2.fastq.sorted.bam  ${prefix}_shSetD2_IP_rep2.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/shSetD2/IP/${prefix}_shSetD2_IP.fastq.merge.bam  ${prefix}_shSetD2_IP_merge.cov.tdf $gsize > /dev/null 2>&1 &

