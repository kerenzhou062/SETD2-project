#!/bin/sh

rm -rf /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shCont/macs1.4/rep1
mkdir -p /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shCont/macs1.4/rep1
cd /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shCont/macs1.4/rep1
macs14 -t /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shCont/IP_rep1/Hela_ChIP-seq_shCont_IP_rep1.fastq.sorted.bam -c /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shCont/input_rep1/Hela_ChIP-seq_shCont_input_rep1.fastq.sorted.bam  -g hs -n shCont_rep1 --nomodel --shiftsize 147 -B -S --call-subpeaks > shCont_rep1.sorted.macs14.log 2>&1 &

rm -rf /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shCont/macs1.4/rep2
mkdir -p /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shCont/macs1.4/rep2
cd /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shCont/macs1.4/rep2
macs14 -t /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shCont/IP_rep2/Hela_ChIP-seq_shCont_IP_rep2.fastq.sorted.bam -c /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shCont/input_rep2/Hela_ChIP-seq_shCont_input_rep2.fastq.sorted.bam  -g hs -n shCont_rep2 --nomodel --shiftsize 147 -B -S --call-subpeaks > shCont_rep2.sorted.macs14.log 2>&1 &


rm -rf /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shSetD2/macs1.4/rep1
mkdir -p /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shSetD2/macs1.4/rep1
cd /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shSetD2/macs1.4/rep1
macs14 -t /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shSetD2/IP_rep1/Hela_ChIP-seq_shSetD2_IP_rep1.fastq.sorted.bam -c /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shSetD2/input_rep1/Hela_ChIP-seq_shSetD2_input_rep1.fastq.sorted.bam  -g hs -n shSetD2_rep1 --nomodel --shiftsize 147 -B -S --call-subpeaks > shSetD2_rep1.sorted.macs14.log 2>&1 &

rm -rf /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shSetD2/macs1.4/rep2
mkdir -p /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shSetD2/macs1.4/rep2
cd /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shSetD2/macs1.4/rep2
macs14 -t /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shSetD2/IP_rep2/Hela_ChIP-seq_shSetD2_IP_rep2.fastq.sorted.bam -c /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shSetD2/input_rep2/Hela_ChIP-seq_shSetD2_input_rep2.fastq.sorted.bam  -g hs -n shSetD2_rep2 --nomodel --shiftsize 147 -B -S --call-subpeaks > shSetD2_rep2.sorted.macs14.log 2>&1 &


rm -rf /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shWTAP/macs1.4
mkdir /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shWTAP/macs1.4
cd /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shWTAP/macs1.4
macs14 -t /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shWTAP/IP_rep1/Hela_ChIP-seq_shWTAP_IP_rep1.fastq.sorted.bam -c /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shWTAP/input_rep1/Hela_ChIP-seq_shWTAP_input_rep1.fastq.sorted.bam  -g hs -n shWTAP_rep1 --nomodel --shiftsize 147 -B -S --call-subpeaks > shWTAP_rep1.macs14.log 2>&1 &


rm -rf /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shM14/macs1.4
mkdir /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shM14/macs1.4
cd /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shM14/macs1.4
macs14 -t /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shM14/IP_rep1/Hela_ChIP-seq_shM14_IP_rep1.fastq.sorted.bam -c /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shM14/input_rep1/Hela_ChIP-seq_shM14_input_rep1.fastq.sorted.bam  -g hs -n shM14_rep1 --nomodel --shiftsize 147 -B -S --call-subpeaks > shM14_rep1.macs14.log 2>&1 &


rm -rf /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shM3/macs1.4
mkdir /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shM3/macs1.4
cd /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shM3/macs1.4
macs14 -t /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shM3/IP_rep1/Hela_ChIP-seq_shM3_IP_rep1.fastq.sorted.bam -c /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shM3/input_rep1/Hela_ChIP-seq_shM3_input_rep1.fastq.sorted.bam  -g hs -n shM3_rep1 --nomodel --shiftsize 147 -B -S --call-subpeaks > shM3_rep1.macs14.log 2>&1 &

