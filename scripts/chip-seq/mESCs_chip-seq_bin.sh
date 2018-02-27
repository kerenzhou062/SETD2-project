#!/bin/sh

cd /data/zhoukr/hhl_setd2_m6a/mouse_ESCs_ChIP-seq/
rm -rf /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/mESCs/bed6
mkdir /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/mESCs/bed6

macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x Ctrl_D0/macs1.4-merge/shCont_peaks.xls -name shCont_macs -o /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/mESCs/bed6/mESCs_macs_shCont_peaks.bed
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x SetD2-KD_D0/macs1.4-merge/shSetD2_peaks.xls -name shSetD2_macs -o /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/mESCs/bed6/mESCs_macs_shSetD2_peaks.bed

cd /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/mESCs/bed6

bed2exclude.pl -a mESCs_macs_shCont_peaks.bed -b mESCs_macs_shSetD2_peaks.bed -o ./mESCs_shCont_SetD2_responsive.bed


#####

###draw bin distribution based on count and bed12 for writer targets

bedBinDistribution.pl --feature 5PStart,geneBody,3PEnd -span 5 -smooth move -t percentage --input mESCs_macs_shCont_peaks.bed -bed6 /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.vM9.mRNA.geneFeature.bed6 -o ./mESCs_chip_shCont.bin
bedBinDistribution.pl --feature 5PStart,geneBody,3PEnd -span 5 -smooth move -t percentage --input mESCs_macs_shSetD2_peaks.bed -bed6 /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.vM9.mRNA.geneFeature.bed6 -o ./mESCs_chip_shSetD2.bin
paste mESCs_chip_shCont.bin mESCs_chip_shSetD2.bin | cut -f 1,2,3,6 > mESCs_chip_all.bin
sed  -i '1i region\tbin\tshCont\tshSetD2' mESCs_chip_all.bin
rm -f mESCs_chip_shCont.bin mESCs_chip_shSetD2.bin

#### gene type counts and region counts of SetD2 dependent H3K36me3 sites

regionDistribution.pl --feature 5utr,cds,stopCodon,3utr --input mESCs_shCont_SetD2_responsive.bed -bed6 /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.vM9.annotation.all.longest.exon.bed6 -o ./mESCs_macs_shCont_SetD2_responsive.region
sed  -i '1i region\tpeakNumber\tenrichment' ./mESCs_macs_shCont_SetD2_responsive.region
