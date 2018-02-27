#!/bin/sh
bamFolder=/data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq
prefix=HepG2_m6A-seq
gsize=hg19.chrom.sizes

rm -rf /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/tdf
mkdir -p /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/tdf
cd /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/tdf
cp /data/zhoukr/reference/genome/genome_size/number/human_hg19_chrsize.txt $gsize


nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/shCont/${prefix}_shCont_input_rep1.fastq.sorted.bam  ${prefix}_shCont_input_rep1.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/shCont/${prefix}_shCont_input_rep2.fastq.sorted.bam  ${prefix}_shCont_input_rep2.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/shCont/${prefix}_shCont_input_rep3.fastq.sorted.bam  ${prefix}_shCont_input_rep3.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/shCont/${prefix}_shCont_IP_rep1.fastq.sorted.bam  ${prefix}_shCont_IP_rep1.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/shCont/${prefix}_shCont_IP_rep2.fastq.sorted.bam  ${prefix}_shCont_IP_rep2.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/shCont/${prefix}_shCont_IP_rep3.fastq.sorted.bam  ${prefix}_shCont_IP_rep3.cov.tdf $gsize > /dev/null 2>&1 &

nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/shSetD2/${prefix}_shSetD2_input_rep1.fastq.sorted.bam  ${prefix}_shSetD2_input_rep1.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/shSetD2/${prefix}_shSetD2_input_rep2.fastq.sorted.bam  ${prefix}_shSetD2_input_rep2.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/shSetD2/${prefix}_shSetD2_input_rep3.fastq.sorted.bam  ${prefix}_shSetD2_input_rep3.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/shSetD2/${prefix}_shSetD2_IP_rep1.fastq.sorted.bam  ${prefix}_shSetD2_IP_rep1.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/shSetD2/${prefix}_shSetD2_IP_rep2.fastq.sorted.bam  ${prefix}_shSetD2_IP_rep2.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/shSetD2/${prefix}_shSetD2_IP_rep3.fastq.sorted.bam  ${prefix}_shSetD2_IP_rep3.cov.tdf $gsize > /dev/null 2>&1 &

nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/shM14/${prefix}_shM14_input_rep1.fastq.sorted.bam  ${prefix}_shM14_input_rep1.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/shM14/${prefix}_shM14_input_rep2.fastq.sorted.bam  ${prefix}_shM14_input_rep2.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/shM14/${prefix}_shM14_input_rep3.fastq.sorted.bam  ${prefix}_shM14_input_rep3.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/shM14/${prefix}_shM14_IP_rep1.fastq.sorted.bam  ${prefix}_shM14_IP_rep1.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/shM14/${prefix}_shM14_IP_rep2.fastq.sorted.bam  ${prefix}_shM14_IP_rep2.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/shM14/${prefix}_shM14_IP_rep3.fastq.sorted.bam  ${prefix}_shM14_IP_rep3.cov.tdf $gsize > /dev/null 2>&1 &

nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/shM3/${prefix}_shM3_input_rep1.fastq.sorted.bam  ${prefix}_shM3_input_rep1.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/shM3/${prefix}_shM3_input_rep2.fastq.sorted.bam  ${prefix}_shM3_input_rep2.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/shM3/${prefix}_shM3_input_rep3.fastq.sorted.bam  ${prefix}_shM3_input_rep3.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/shM3/${prefix}_shM3_IP_rep1.fastq.sorted.bam  ${prefix}_shM3_IP_rep1.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/shM3/${prefix}_shM3_IP_rep2.fastq.sorted.bam  ${prefix}_shM3_IP_rep2.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/shM3/${prefix}_shM3_IP_rep3.fastq.sorted.bam  ${prefix}_shM3_IP_rep3.cov.tdf $gsize > /dev/null 2>&1 &

nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/shWTAP/${prefix}_shWTAP_input_rep1.fastq.sorted.bam  ${prefix}_shWTAP_input_rep1.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/shWTAP/${prefix}_shWTAP_input_rep2.fastq.sorted.bam  ${prefix}_shWTAP_input_rep2.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/shWTAP/${prefix}_shWTAP_input_rep3.fastq.sorted.bam  ${prefix}_shWTAP_input_rep3.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/shWTAP/${prefix}_shWTAP_IP_rep1.fastq.sorted.bam  ${prefix}_shWTAP_IP_rep1.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/shWTAP/${prefix}_shWTAP_IP_rep2.fastq.sorted.bam  ${prefix}_shWTAP_IP_rep2.cov.tdf $gsize > /dev/null 2>&1 &
nohup igvtools count -z 7 -w 25 -e 250 $bamFolder/shWTAP/${prefix}_shWTAP_IP_rep3.fastq.sorted.bam  ${prefix}_shWTAP_IP_rep3.cov.tdf $gsize > /dev/null 2>&1 &
