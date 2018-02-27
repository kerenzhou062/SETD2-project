#!/bin/sh

cd /data/zhoukr/hhl_setd2_m6a/mouse_MEFs_ChIP-seq/
rm -rf /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/mMEFs/bed6
mkdir /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/mMEFs/bed6

macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x SetD2-WT_NA/macs1.4-merge/shCont_peaks.xls -name shCont_macs -o /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/mMEFs/bed6/mMEFs_macs_shCont_peaks.bed
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x SetD2-KO_NA/macs1.4-merge/shSetD2_peaks.xls -name shSetD2_macs -o /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/mMEFs/bed6/mMEFs_macs_shSetD2_peaks.bed

cd /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/mMEFs/bed6

bed2exclude.pl -a mMEFs_macs_shCont_peaks.bed -b mMEFs_macs_shSetD2_peaks.bed -o ./mMEFs_shCont_SetD2_responsive.bed


#####

###draw bin distribution based on count and bed12 for writer targets

bedBinDistribution.pl --feature 5PStart,geneBody,3PEnd -span 5 -smooth move -t count --input mMEFs_macs_shCont_peaks.bed -bed6 /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.vM9.mRNA.geneFeature.bed6 -o ./mMEFs_chip_shCont.bin
bedBinDistribution.pl --feature 5PStart,geneBody,3PEnd -span 5 -smooth move -t count --input mMEFs_macs_shSetD2_peaks.bed -bed6 /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.vM9.mRNA.geneFeature.bed6 -o ./mMEFs_chip_shSetD2.bin
paste mMEFs_chip_shCont.bin mMEFs_chip_shSetD2.bin | cut -f 1,2,3,6 > mMEFs_chip_all.bin
sed  -i '1i region\tbin\tshCont\tshSetD2' mMEFs_chip_all.bin
rm -f mMEFs_chip_shCont.bin mMEFs_chip_shSetD2.bin

#### gene type counts and region counts of SetD2 responsive H3K36me3 sites

regionDistribution.pl --feature 5utr,cds,stopCodon,3utr --input mMEFs_shCont_SetD2_responsive.bed -bed6 /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.vM9.annotation.all.longest.exon.bed6 -o ./mMEFs_macs_shCont_SetD2_responsive.region
sed  -i '1i region\tpeakNumber\tenrichment' ./mMEFs_macs_shCont_SetD2_responsive.region
