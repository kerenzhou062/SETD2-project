#!/bin/sh
bamFolder=/data/zhoukr/hhl_setd2_m6a/mouse_ESCs_ChIP-seq
prefix=mouse_ESCs_ChIP-seq
gsize=mm10.chrom.sizes

rm -rf /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/mESCs/tdf
mkdir -p /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/mESCs/tdf
cd /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/mESCs/tdf
cp /data/zhoukr/reference/genome/genome_size/number/mouse_mm10_chrsize.txt $gsize


nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/Ctrl_D0/${prefix}_Ctrl_D0_input.sorted.bam  ${prefix}_Ctrl_D0_input.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/Ctrl_D0/${prefix}_Ctrl_D0_IP.sorted.bam  ${prefix}_Ctrl_D0_IP.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/Ctrl_D6/${prefix}_Ctrl_D6_input.sorted.bam  ${prefix}_Ctrl_D6_input.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/Ctrl_D6/${prefix}_Ctrl_D6_IP.sorted.bam  ${prefix}_Ctrl_D6_IP.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/SetD2-KD_D0/${prefix}_SetD2-KD_D0_input.sorted.bam  ${prefix}_SetD2-KD_D0_input.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/SetD2-KD_D0/${prefix}_SetD2-KD_D0_IP.sorted.bam  ${prefix}_SetD2-KD_D0_IP.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/SetD2-KD_D6/${prefix}_SetD2-KD_D6_input.sorted.bam  ${prefix}_SetD2-KD_D6_input.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/SetD2-KD_D6/${prefix}_SetD2-KD_D6_IP.sorted.bam  ${prefix}_SetD2-KD_D6_IP.cov.tdf $gsize > /dev/null 2>&1 &


