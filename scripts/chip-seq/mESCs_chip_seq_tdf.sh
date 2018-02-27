#!/bin/sh
bamFolder=/data/zhoukr/hhl_setd2_m6a/mouse_ESCs_ChIP-seq
prefix=mESCs_ChIP-seq
gsize=mm10.chrom.sizes

rm -rf /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/mESCs/tdf
mkdir -p /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/mESCs/tdf
cd /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/mESCs/tdf
cp /data/zhoukr/reference/genome/genome_size/number/mouse_mm10_chrsize.txt $gsize

nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/Ctrl_D0/${prefix}_Ctrl_D0_ChIP-rep1.sorted.bam  ${prefix}_Ctrl_D0_ChIP-rep1.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/Ctrl_D0/${prefix}_Ctrl_D0_ChIP-rep2.sorted.bam  ${prefix}_Ctrl_D0_ChIP-rep2.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/Ctrl_D0/${prefix}_Ctrl_D0_input-rep1.sorted.bam  ${prefix}_Ctrl_D0_input-rep1.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/Ctrl_D0/${prefix}_Ctrl_D0_input-rep2.sorted.bam  ${prefix}_Ctrl_D0_input-rep2.cov.tdf $gsize > /dev/null 2>&1 &

nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/SetD2-KD_D0/${prefix}_SetD2-KD_D0_ChIP-rep1.sorted.bam  ${prefix}_SetD2-KD_D0_ChIP-rep1.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/SetD2-KD_D0/${prefix}_SetD2-KD_D0_ChIP-rep2.sorted.bam  ${prefix}_SetD2-KD_D0_ChIP-rep2.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/SetD2-KD_D0/${prefix}_SetD2-KD_D0_input-rep1.sorted.bam  ${prefix}_SetD2-KD_D0_input-rep1.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/SetD2-KD_D0/${prefix}_SetD2-KD_D0_input-rep2.sorted.bam  ${prefix}_SetD2-KD_D0_input-rep2.cov.tdf $gsize > /dev/null 2>&1 &
