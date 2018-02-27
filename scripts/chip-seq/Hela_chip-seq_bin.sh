#!/bin/sh

cd /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/
rm -rf /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/Hela/bed6
mkdir /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/Hela/bed6

macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x shCont/macs1.4/rep2/shCont_rep2_peaks.xls -name shCont_macs -o /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/Hela/bed6/Hela_macs_shCont_peaks.bed
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x shSetD2/macs1.4/rep2/shSetD2_rep2_peaks.xls -name shSetD2_macs -o /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/Hela/bed6/Hela_macs_shSetD2_peaks.bed
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x shM14/macs1.4/shM14_rep1_peaks.xls -name shM14_macs -o /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/Hela/bed6/Hela_macs_shM14_peaks.bed
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x shM3/macs1.4/shM3_rep1_peaks.xls -name shM3_macs -o /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/Hela/bed6/Hela_macs_shM3_peaks.bed
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x shWTAP/macs1.4/shWTAP_rep1_peaks.xls -name shWTAP_macs -o /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/Hela/bed6/Hela_macs_shWTAP_peaks.bed

cd /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/Hela/bed6

bed2exclude.pl -a Hela_macs_shCont_peaks.bed -b Hela_macs_shSetD2_peaks.bed -o ./Hela_shCont_SetD2_responsive.bed

bed2exclude.pl -a Hela_macs_shCont_peaks.bed -b Hela_macs_shM14_peaks.bed -o ./Hela_shCont_M14_responsive.bed
bed2overlap.pl -aOver -a Hela_macs_shSetD2_peaks.bed -b Hela_shCont_M14_responsive.bed -o ./Hela_shSetD2_M14_responsive.bed

bed2exclude.pl -a Hela_macs_shCont_peaks.bed -b Hela_macs_shM3_peaks.bed -o ./Hela_shCont_M3_responsive.bed
bed2overlap.pl -aOver -a Hela_macs_shSetD2_peaks.bed -b Hela_shCont_M3_responsive.bed -o ./Hela_shSetD2_M3_responsive.bed

bed2exclude.pl -a Hela_macs_shCont_peaks.bed -b Hela_macs_shWTAP_peaks.bed -o ./Hela_shCont_WTAP_responsive.bed
bed2overlap.pl -aOver -a Hela_macs_shSetD2_peaks.bed -b Hela_shCont_WTAP_responsive.bed -o ./Hela_shSetD2_WTAP_responsive.bed

#####
cat Hela_shCont_M14_responsive.bed Hela_shCont_M3_responsive.bed Hela_shCont_WTAP_responsive.bed | sort -t $'\t' -k1,1V -k2,2n | uniq > Hela_shCont_allWriter_targets.bed
bed2exclude.pl -a Hela_macs_shCont_peaks.bed -b Hela_shCont_allWriter_targets.bed -o ./Hela_chip_shCont_nonWriter_responsive.bed
bed2overlap.pl -aOver -a Hela_macs_shSetD2_peaks.bed -b Hela_chip_shCont_nonWriter_responsive.bed -o ./Hela_chip_shSetD2_nonWriter_responsive.bed

###draw bin distribution based on count and bed12 for writer targets

bedBinDistribution.pl --feature 5PStart,geneBody,3PEnd -span 5 -smooth move -t count --input Hela_shCont_M14_responsive.bed -bed6 /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.mRNA.geneFeature.bed6 -o ./Hela_shCont_M14_responsive.bin
bedBinDistribution.pl --feature 5PStart,geneBody,3PEnd -span 5 -smooth move -t count --input Hela_shSetD2_M14_responsive.bed -bed6 /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.mRNA.geneFeature.bed6 -o ./Hela_shSetD2_M14_responsive.bin
paste Hela_shCont_M14_responsive.bin Hela_shSetD2_M14_responsive.bin | cut -f 1,2,3,6 > Hela_chip_M14_responsive.bin
sed  -i '1i region\tbin\tshCont\tshSetD2' Hela_chip_M14_responsive.bin
rm -f Hela_shCont_M14_responsive.bin Hela_shSetD2_M14_responsive.bin

bedBinDistribution.pl --feature 5PStart,geneBody,3PEnd -span 5 -smooth move -t count --input Hela_shCont_M3_responsive.bed -bed6 /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.mRNA.geneFeature.bed6 -o ./Hela_shCont_M3_responsive.bin
bedBinDistribution.pl --feature 5PStart,geneBody,3PEnd -span 5 -smooth move -t count --input Hela_shSetD2_M3_responsive.bed -bed6 /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.mRNA.geneFeature.bed6 -o ./Hela_shSetD2_M3_responsive.bin
paste Hela_shCont_M3_responsive.bin Hela_shSetD2_M3_responsive.bin | cut -f 1,2,3,6 > Hela_chip_M3_responsive.bin
sed  -i '1i region\tbin\tshCont\tshSetD2' Hela_chip_M3_responsive.bin
rm -f Hela_shCont_M3_responsive.bin Hela_shSetD2_M3_responsive.bin

bedBinDistribution.pl --feature 5PStart,geneBody,3PEnd -span 5 -smooth move -t count --input Hela_shCont_WTAP_responsive.bed -bed6 /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.mRNA.geneFeature.bed6 -o ./Hela_shCont_WTAP_responsive.bin
bedBinDistribution.pl --feature 5PStart,geneBody,3PEnd -span 5 -smooth move -t count --input Hela_shSetD2_WTAP_responsive.bed -bed6 /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.mRNA.geneFeature.bed6 -o ./Hela_shSetD2_WTAP_responsive.bin
paste Hela_shCont_WTAP_responsive.bin Hela_shSetD2_WTAP_responsive.bin | cut -f 1,2,3,6 > Hela_chip_WTAP_responsive.bin
sed  -i '1i region\tbin\tshCont\tshSetD2' Hela_chip_WTAP_responsive.bin
rm -f Hela_shCont_WTAP_responsive.bin Hela_shSetD2_WTAP_responsive.bin

bedBinDistribution.pl --feature 5PStart,geneBody,3PEnd -span 5 -smooth move -t count --input Hela_chip_shCont_nonWriter_responsive.bed -bed6 /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.mRNA.geneFeature.bed6 -o ./Hela_chip_shCont_nonWriter_responsive.bin
bedBinDistribution.pl --feature 5PStart,geneBody,3PEnd -span 5 -smooth move -t count --input Hela_chip_shSetD2_nonWriter_responsive.bed -bed6 /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.mRNA.geneFeature.bed6 -o ./Hela_chip_shSetD2_nonWriter_responsive.bin
paste Hela_chip_shCont_nonWriter_responsive.bin Hela_chip_shSetD2_nonWriter_responsive.bin | cut -f 1,2,3,6 > Hela_chip_nonWriter_responsive.bin
sed  -i '1i region\tbin\tshCont\tshSetD2' Hela_chip_nonWriter_responsive.bin
rm -f Hela_chip_shCont_nonWriter_responsive.bin Hela_chip_shSetD2_nonWriter_responsive.bin

bedBinDistribution.pl --feature 5PStart,geneBody,3PEnd -span 5 -smooth move -t percentage --input Hela_macs_shCont_peaks.bed -bed6 /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.mRNA.geneFeature.bed6 -o ./Hela_chip_shCont.bin
bedBinDistribution.pl --feature 5PStart,geneBody,3PEnd -span 5 -smooth move -t percentage --input Hela_macs_shSetD2_peaks.bed -bed6 /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.mRNA.geneFeature.bed6 -o ./Hela_chip_shSetD2.bin
bedBinDistribution.pl --feature 5PStart,geneBody,3PEnd -span 5 -smooth move -t percentage --input Hela_macs_shM14_peaks.bed -bed6 /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.mRNA.geneFeature.bed6 -o ./Hela_chip_shM14.bin
bedBinDistribution.pl --feature 5PStart,geneBody,3PEnd -span 5 -smooth move -t percentage --input Hela_macs_shM3_peaks.bed -bed6 /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.mRNA.geneFeature.bed6 -o ./Hela_chip_shM3.bin
bedBinDistribution.pl --feature 5PStart,geneBody,3PEnd -span 5 -smooth move -t percentage --input Hela_macs_shWTAP_peaks.bed -bed6 /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.mRNA.geneFeature.bed6 -o ./Hela_chip_shWTAP.bin
paste Hela_chip_shCont.bin Hela_chip_shSetD2.bin Hela_chip_shM14.bin Hela_chip_shM3.bin Hela_chip_shWTAP.bin | cut -f 1,2,3,6,9,12,15 > Hela_chip_all.bin
sed  -i '1i region\tbin\tshCont\tshSetD2\tshM14\tshM3\tshWTAP' Hela_chip_all.bin
rm -f Hela_chip_shCont.bin Hela_chip_shSetD2.bin Hela_chip_shM14.bin Hela_chip_shM3.bin Hela_chip_shWTAP.bin

#### gene type counts and region counts of SetD2 responsive H3K36me3 sites

regionDistribution.pl --feature 5utr,cds,stopCodon,3utr --input Hela_shCont_SetD2_responsive.bed -bed6 /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.mRNA.longest.exon+promoter.bed6 -o ./Hela_macs_shCont_SetD2_responsive.region
sed  -i '1i region\tpeakNumber\tenrichment' ./Hela_macs_shCont_SetD2_responsive.region
