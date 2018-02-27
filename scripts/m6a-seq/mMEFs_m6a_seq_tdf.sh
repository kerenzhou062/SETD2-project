#!/bin/sh
bamFolder=/data/zhoukr/hhl_setd2_m6a/mouse_MEFs_m6A-seq
prefix=mMEFs_m6A-seq
gsize=mm10.chrom.sizes

rm -rf /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/mMEFs/tdf
mkdir -p /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/mMEFs/tdf
cd /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/mMEFs/tdf
cp /data/zhoukr/reference/genome/genome_size/number/mouse_mm10_chrsize.txt $gsize


nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/SetD2-WT_NA/${prefix}_SetD2-WT_NA_input.sorted.bam  ${prefix}_SetD2-WT_NA_input.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/SetD2-WT_NA/${prefix}_SetD2-WT_NA_IP.sorted.bam  ${prefix}_SetD2-WT_NA_IP.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/SetD2-KO_NA/${prefix}_SetD2-KO_NA_input.sorted.bam  ${prefix}_SetD2-KO_NA_input.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/SetD2-KO_NA/${prefix}_SetD2-KO_NA_IP.sorted.bam  ${prefix}_SetD2-KO_NA_IP.cov.tdf $gsize > /dev/null 2>&1 &
