#!/bin/sh

exomePeakPath=/data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq
rm -rf /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/statistics
mkdir /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/statistics
cd /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/statistics

## rep1
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x $exomePeakPath/shCont/shCont_m6A_rep1/con_peak.xls -o HepG2_shCont_m6a.tmp
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x $exomePeakPath/shSetD2/shSetD2_m6A_rep1/con_peak.xls -o HepG2_shSetD2_m6a.tmp
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x $exomePeakPath/shM14/shM14_m6A_rep1/con_peak.xls -o HepG2_shM14_m6a.tmp
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x $exomePeakPath/shM3/shM3_m6A_rep1/con_peak.xls -o HepG2_shM3_m6a.tmp
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x $exomePeakPath/shWTAP/shWTAP_m6A_rep1/con_peak.xls -o HepG2_shWTAP_m6a.tmp

awk '{key=$1""$2""$3; hash[key]=$1"\t"$2"\t"$3"\t"$6;}END{for(i in hash){print hash[i];}}' OFS="\t" HepG2_shCont_m6a.tmp | sort | uniq > HepG2_shCont_m6a_rep1.bed
awk '{key=$1""$2""$3; hash[key]=$1"\t"$2"\t"$3"\t"$6;}END{for(i in hash){print hash[i];}}' OFS="\t" HepG2_shSetD2_m6a.tmp | sort | uniq > HepG2_shSetD2_m6a_rep1.bed
awk '{key=$1""$2""$3; hash[key]=$1"\t"$2"\t"$3"\t"$6;}END{for(i in hash){print hash[i];}}' OFS="\t" HepG2_shM14_m6a.tmp | sort | uniq > HepG2_shM14_m6a_rep1.bed
awk '{key=$1""$2""$3; hash[key]=$1"\t"$2"\t"$3"\t"$6;}END{for(i in hash){print hash[i];}}' OFS="\t" HepG2_shM3_m6a.tmp | sort | uniq > HepG2_shM3_m6a_rep1.bed
awk '{key=$1""$2""$3; hash[key]=$1"\t"$2"\t"$3"\t"$6;}END{for(i in hash){print hash[i];}}' OFS="\t" HepG2_shWTAP_m6a.tmp | sort | uniq > HepG2_shWTAP_m6a_rep1.bed

## rep2
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x $exomePeakPath/shCont/shCont_m6A_rep2/con_peak.xls -o HepG2_shCont_m6a.tmp
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x $exomePeakPath/shSetD2/shSetD2_m6A_rep2/con_peak.xls -o HepG2_shSetD2_m6a.tmp
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x $exomePeakPath/shM14/shM14_m6A_rep2/con_peak.xls -o HepG2_shM14_m6a.tmp
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x $exomePeakPath/shM3/shM3_m6A_rep2/con_peak.xls -o HepG2_shM3_m6a.tmp
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x $exomePeakPath/shWTAP/shWTAP_m6A_rep2/con_peak.xls -o HepG2_shWTAP_m6a.tmp

awk '{key=$1""$2""$3; hash[key]=$1"\t"$2"\t"$3"\t"$6;}END{for(i in hash){print hash[i];}}' OFS="\t" HepG2_shCont_m6a.tmp | sort | uniq > HepG2_shCont_m6a_rep2.bed
awk '{key=$1""$2""$3; hash[key]=$1"\t"$2"\t"$3"\t"$6;}END{for(i in hash){print hash[i];}}' OFS="\t" HepG2_shSetD2_m6a.tmp | sort | uniq > HepG2_shSetD2_m6a_rep2.bed
awk '{key=$1""$2""$3; hash[key]=$1"\t"$2"\t"$3"\t"$6;}END{for(i in hash){print hash[i];}}' OFS="\t" HepG2_shM14_m6a.tmp | sort | uniq > HepG2_shM14_m6a_rep2.bed
awk '{key=$1""$2""$3; hash[key]=$1"\t"$2"\t"$3"\t"$6;}END{for(i in hash){print hash[i];}}' OFS="\t" HepG2_shM3_m6a.tmp | sort | uniq > HepG2_shM3_m6a_rep2.bed
awk '{key=$1""$2""$3; hash[key]=$1"\t"$2"\t"$3"\t"$6;}END{for(i in hash){print hash[i];}}' OFS="\t" HepG2_shWTAP_m6a.tmp | sort | uniq > HepG2_shWTAP_m6a_rep2.bed

## rep3
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x $exomePeakPath/shCont/shCont_m6A_rep3/con_peak.xls -o HepG2_shCont_m6a.tmp
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x $exomePeakPath/shSetD2/shSetD2_m6A_rep3/con_peak.xls -o HepG2_shSetD2_m6a.tmp
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x $exomePeakPath/shM14/shM14_m6A_rep3/con_peak.xls -o HepG2_shM14_m6a.tmp
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x $exomePeakPath/shM3/shM3_m6A_rep3/con_peak.xls -o HepG2_shM3_m6a.tmp
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x $exomePeakPath/shWTAP/shWTAP_m6A_rep3/con_peak.xls -o HepG2_shWTAP_m6a.tmp

awk '{key=$1""$2""$3; hash[key]=$1"\t"$2"\t"$3"\t"$6;}END{for(i in hash){print hash[i];}}' OFS="\t" HepG2_shCont_m6a.tmp | sort | uniq > HepG2_shCont_m6a_rep3.bed
awk '{key=$1""$2""$3; hash[key]=$1"\t"$2"\t"$3"\t"$6;}END{for(i in hash){print hash[i];}}' OFS="\t" HepG2_shSetD2_m6a.tmp | sort | uniq > HepG2_shSetD2_m6a_rep3.bed
awk '{key=$1""$2""$3; hash[key]=$1"\t"$2"\t"$3"\t"$6;}END{for(i in hash){print hash[i];}}' OFS="\t" HepG2_shM14_m6a.tmp | sort | uniq > HepG2_shM14_m6a_rep3.bed
awk '{key=$1""$2""$3; hash[key]=$1"\t"$2"\t"$3"\t"$6;}END{for(i in hash){print hash[i];}}' OFS="\t" HepG2_shM3_m6a.tmp | sort | uniq > HepG2_shM3_m6a_rep3.bed
awk '{key=$1""$2""$3; hash[key]=$1"\t"$2"\t"$3"\t"$6;}END{for(i in hash){print hash[i];}}' OFS="\t" HepG2_shWTAP_m6a.tmp | sort | uniq > HepG2_shWTAP_m6a_rep3.bed

## rep1-2
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x $exomePeakPath/shCont/shCont_m6A_rep1-2/con_peak.xls -o HepG2_shCont_m6a.tmp
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x $exomePeakPath/shSetD2/shSetD2_m6A_rep1-2/con_peak.xls -o HepG2_shSetD2_m6a.tmp
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x $exomePeakPath/shM14/shM14_m6A_rep1-2/con_peak.xls -o HepG2_shM14_m6a.tmp
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x $exomePeakPath/shM3/shM3_m6A_rep1-2/con_peak.xls -o HepG2_shM3_m6a.tmp
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x $exomePeakPath/shWTAP/shWTAP_m6A_rep1-2/con_peak.xls -o HepG2_shWTAP_m6a.tmp

awk '{key=$1""$2""$3; hash[key]=$1"\t"$2"\t"$3"\t"$6;}END{for(i in hash){print hash[i];}}' OFS="\t" HepG2_shCont_m6a.tmp | sort | uniq > HepG2_shCont_m6a_rep1-2.bed
awk '{key=$1""$2""$3; hash[key]=$1"\t"$2"\t"$3"\t"$6;}END{for(i in hash){print hash[i];}}' OFS="\t" HepG2_shSetD2_m6a.tmp | sort | uniq > HepG2_shSetD2_m6a_rep1-2.bed
awk '{key=$1""$2""$3; hash[key]=$1"\t"$2"\t"$3"\t"$6;}END{for(i in hash){print hash[i];}}' OFS="\t" HepG2_shM14_m6a.tmp | sort | uniq > HepG2_shM14_m6a_rep1-2.bed
awk '{key=$1""$2""$3; hash[key]=$1"\t"$2"\t"$3"\t"$6;}END{for(i in hash){print hash[i];}}' OFS="\t" HepG2_shM3_m6a.tmp | sort | uniq > HepG2_shM3_m6a_rep1-2.bed
awk '{key=$1""$2""$3; hash[key]=$1"\t"$2"\t"$3"\t"$6;}END{for(i in hash){print hash[i];}}' OFS="\t" HepG2_shWTAP_m6a.tmp | sort | uniq > HepG2_shWTAP_m6a_rep1-2.bed

## rep1-3
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x $exomePeakPath/shCont/shCont_m6A_rep1-3/con_peak.xls -o HepG2_shCont_m6a.tmp
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x $exomePeakPath/shSetD2/shSetD2_m6A_rep1-3/con_peak.xls -o HepG2_shSetD2_m6a.tmp
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x $exomePeakPath/shM14/shM14_m6A_rep1-3/con_peak.xls -o HepG2_shM14_m6a.tmp
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x $exomePeakPath/shM3/shM3_m6A_rep1-3/con_peak.xls -o HepG2_shM3_m6a.tmp
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x $exomePeakPath/shWTAP/shWTAP_m6A_rep1-3/con_peak.xls -o HepG2_shWTAP_m6a.tmp

awk '{key=$1""$2""$3; hash[key]=$1"\t"$2"\t"$3"\t"$6;}END{for(i in hash){print hash[i];}}' OFS="\t" HepG2_shCont_m6a.tmp | sort | uniq > HepG2_shCont_m6a_rep1-3.bed
awk '{key=$1""$2""$3; hash[key]=$1"\t"$2"\t"$3"\t"$6;}END{for(i in hash){print hash[i];}}' OFS="\t" HepG2_shSetD2_m6a.tmp | sort | uniq > HepG2_shSetD2_m6a_rep1-3.bed
awk '{key=$1""$2""$3; hash[key]=$1"\t"$2"\t"$3"\t"$6;}END{for(i in hash){print hash[i];}}' OFS="\t" HepG2_shM14_m6a.tmp | sort | uniq > HepG2_shM14_m6a_rep1-3.bed
awk '{key=$1""$2""$3; hash[key]=$1"\t"$2"\t"$3"\t"$6;}END{for(i in hash){print hash[i];}}' OFS="\t" HepG2_shM3_m6a.tmp | sort | uniq > HepG2_shM3_m6a_rep1-3.bed
awk '{key=$1""$2""$3; hash[key]=$1"\t"$2"\t"$3"\t"$6;}END{for(i in hash){print hash[i];}}' OFS="\t" HepG2_shWTAP_m6a.tmp | sort | uniq > HepG2_shWTAP_m6a_rep1-3.bed

## rep2-3
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x $exomePeakPath/shCont/shCont_m6A_rep2-3/con_peak.xls -o HepG2_shCont_m6a.tmp
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x $exomePeakPath/shSetD2/shSetD2_m6A_rep2-3/con_peak.xls -o HepG2_shSetD2_m6a.tmp
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x $exomePeakPath/shM14/shM14_m6A_rep2-3/con_peak.xls -o HepG2_shM14_m6a.tmp
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x $exomePeakPath/shM3/shM3_m6A_rep2-3/con_peak.xls -o HepG2_shM3_m6a.tmp
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x $exomePeakPath/shWTAP/shWTAP_m6A_rep2-3/con_peak.xls -o HepG2_shWTAP_m6a.tmp

awk '{key=$1""$2""$3; hash[key]=$1"\t"$2"\t"$3"\t"$6;}END{for(i in hash){print hash[i];}}' OFS="\t" HepG2_shCont_m6a.tmp | sort | uniq > HepG2_shCont_m6a_rep2-3.bed
awk '{key=$1""$2""$3; hash[key]=$1"\t"$2"\t"$3"\t"$6;}END{for(i in hash){print hash[i];}}' OFS="\t" HepG2_shSetD2_m6a.tmp | sort | uniq > HepG2_shSetD2_m6a_rep2-3.bed
awk '{key=$1""$2""$3; hash[key]=$1"\t"$2"\t"$3"\t"$6;}END{for(i in hash){print hash[i];}}' OFS="\t" HepG2_shM14_m6a.tmp | sort | uniq > HepG2_shM14_m6a_rep2-3.bed
awk '{key=$1""$2""$3; hash[key]=$1"\t"$2"\t"$3"\t"$6;}END{for(i in hash){print hash[i];}}' OFS="\t" HepG2_shM3_m6a.tmp | sort | uniq > HepG2_shM3_m6a_rep2-3.bed
awk '{key=$1""$2""$3; hash[key]=$1"\t"$2"\t"$3"\t"$6;}END{for(i in hash){print hash[i];}}' OFS="\t" HepG2_shWTAP_m6a.tmp | sort | uniq > HepG2_shWTAP_m6a_rep2-3.bed

## rep1-2-3
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x $exomePeakPath/shCont/shCont_m6A_rep1-2-3/con_peak.xls -o HepG2_shCont_m6a.tmp
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x $exomePeakPath/shSetD2/shSetD2_m6A_rep1-2-3/con_peak.xls -o HepG2_shSetD2_m6a.tmp
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x $exomePeakPath/shM14/shM14_m6A_rep1-2-3/con_peak.xls -o HepG2_shM14_m6a.tmp
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x $exomePeakPath/shM3/shM3_m6A_rep1-2-3/con_peak.xls -o HepG2_shM3_m6a.tmp
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x $exomePeakPath/shWTAP/shWTAP_m6A_rep1-2-3/con_peak.xls -o HepG2_shWTAP_m6a.tmp

awk '{key=$1""$2""$3; hash[key]=$1"\t"$2"\t"$3"\t"$6;}END{for(i in hash){print hash[i];}}' OFS="\t" HepG2_shCont_m6a.tmp | sort | uniq > HepG2_shCont_m6a_rep1-2-3.bed
awk '{key=$1""$2""$3; hash[key]=$1"\t"$2"\t"$3"\t"$6;}END{for(i in hash){print hash[i];}}' OFS="\t" HepG2_shSetD2_m6a.tmp | sort | uniq > HepG2_shSetD2_m6a_rep1-2-3.bed
awk '{key=$1""$2""$3; hash[key]=$1"\t"$2"\t"$3"\t"$6;}END{for(i in hash){print hash[i];}}' OFS="\t" HepG2_shM14_m6a.tmp | sort | uniq > HepG2_shM14_m6a_rep1-2-3.bed
awk '{key=$1""$2""$3; hash[key]=$1"\t"$2"\t"$3"\t"$6;}END{for(i in hash){print hash[i];}}' OFS="\t" HepG2_shM3_m6a.tmp | sort | uniq > HepG2_shM3_m6a_rep1-2-3.bed
awk '{key=$1""$2""$3; hash[key]=$1"\t"$2"\t"$3"\t"$6;}END{for(i in hash){print hash[i];}}' OFS="\t" HepG2_shWTAP_m6a.tmp | sort | uniq > HepG2_shWTAP_m6a_rep1-2-3.bed

find ./ -type f -name "*.tmp" | xargs -I {} rm -f {}
