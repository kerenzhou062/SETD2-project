#!/bin/sh

#binAnno=/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.mRNA.geneFeature.bed6
#binFeature='5PStart,geneBody,3PEnd'
binAnno=/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.mRNA.longest.exon+promoter.bed6
binFeature='promoter,5utr,cds,3utr'
regionAnno=/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.mRNA.longest.exon+promoter.bed6
regionFeature='promoter,5utr,cds,stopCodon,3utr'

cd /data/zhoukr/hhl_setd2_m6a/HepG2_ChIP-seq/
rm -rf /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/bed6
mkdir /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/bed6

macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x shCont/macs1.4-merge/shCont_peaks.xls -name shCont_macs -o /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/bed6/HepG2_macs_shCont_peaks.bed
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x shSetD2/macs1.4-merge/shSetD2_peaks.xls -name shSetD2_macs -o /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/bed6/HepG2_macs_shSetD2_peaks.bed
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x shM14/macs1.4/shM14_rep1_peaks.xls -name shM14_macs -o /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/bed6/HepG2_macs_shM14_peaks.bed
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x shM3/macs1.4/shM3_rep1_peaks.xls -name shM3_macs -o /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/bed6/HepG2_macs_shM3_peaks.bed
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x shWTAP/macs1.4/shWTAP_rep1_peaks.xls -name shWTAP_macs -o /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/bed6/HepG2_macs_shWTAP_peaks.bed

cd /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/bed6

bed2exclude.pl -a HepG2_macs_shCont_peaks.bed -b HepG2_macs_shSetD2_peaks.bed -o ./HepG2_shCont_SetD2_responsive.bed

bed2exclude.pl -a HepG2_macs_shCont_peaks.bed -b HepG2_macs_shM14_peaks.bed -o ./HepG2_shCont_M14_responsive.bed
bed2overlap.pl -aOver -a HepG2_macs_shSetD2_peaks.bed -b HepG2_shCont_M14_responsive.bed -o ./HepG2_shSetD2_M14_responsive.bed

bed2exclude.pl -a HepG2_macs_shCont_peaks.bed -b HepG2_macs_shM3_peaks.bed -o ./HepG2_shCont_M3_responsive.bed
bed2overlap.pl -aOver -a HepG2_macs_shSetD2_peaks.bed -b HepG2_shCont_M3_responsive.bed -o ./HepG2_shSetD2_M3_responsive.bed

bed2exclude.pl -a HepG2_macs_shCont_peaks.bed -b HepG2_macs_shWTAP_peaks.bed -o ./HepG2_shCont_WTAP_responsive.bed
bed2overlap.pl -aOver -a HepG2_macs_shSetD2_peaks.bed -b HepG2_shCont_WTAP_responsive.bed -o ./HepG2_shSetD2_WTAP_responsive.bed

#####
cat HepG2_shCont_M14_responsive.bed HepG2_shCont_M3_responsive.bed HepG2_shCont_WTAP_responsive.bed | sort -t $'\t' -k1,1V -k2,2n | uniq > HepG2_shCont_allWriter_targets.bed
bed2exclude.pl -a HepG2_macs_shCont_peaks.bed -b HepG2_shCont_allWriter_targets.bed -o ./HepG2_chip_shCont_nonWriter_responsive.bed
bed2overlap.pl -aOver -a HepG2_macs_shSetD2_peaks.bed -b HepG2_chip_shCont_nonWriter_responsive.bed -o ./HepG2_chip_shSetD2_nonWriter_responsive.bed

###draw bin distribution based on count and bed12 for writer targets

bedBinDistribution.pl --feature $binFeature -span 5 -smooth move -t count --input HepG2_shCont_M14_responsive.bed -bed6 $binAnno -o ./HepG2_shCont_M14_responsive.bin
bedBinDistribution.pl --feature $binFeature -span 5 -smooth move -t count --input HepG2_shSetD2_M14_responsive.bed -bed6 $binAnno -o ./HepG2_shSetD2_M14_responsive.bin
paste HepG2_shCont_M14_responsive.bin HepG2_shSetD2_M14_responsive.bin | cut -f 1,2,3,6 > HepG2_chip_M14_responsive.bin
sed  -i '1i region\tbin\tshCont\tshSetD2' HepG2_chip_M14_responsive.bin
rm -f HepG2_shCont_M14_responsive.bin HepG2_shSetD2_M14_responsive.bin

bedBinDistribution.pl --feature $binFeature -span 5 -smooth move -t count --input HepG2_shCont_M3_responsive.bed -bed6 $binAnno -o ./HepG2_shCont_M3_responsive.bin
bedBinDistribution.pl --feature $binFeature -span 5 -smooth move -t count --input HepG2_shSetD2_M3_responsive.bed -bed6 $binAnno -o ./HepG2_shSetD2_M3_responsive.bin
paste HepG2_shCont_M3_responsive.bin HepG2_shSetD2_M3_responsive.bin | cut -f 1,2,3,6 > HepG2_chip_M3_responsive.bin
sed  -i '1i region\tbin\tshCont\tshSetD2' HepG2_chip_M3_responsive.bin
rm -f HepG2_shCont_M3_responsive.bin HepG2_shSetD2_M3_responsive.bin

bedBinDistribution.pl --feature $binFeature -span 5 -smooth move -t count --input HepG2_shCont_WTAP_responsive.bed -bed6 $binAnno -o ./HepG2_shCont_WTAP_responsive.bin
bedBinDistribution.pl --feature $binFeature -span 5 -smooth move -t count --input HepG2_shSetD2_WTAP_responsive.bed -bed6 $binAnno -o ./HepG2_shSetD2_WTAP_responsive.bin
paste HepG2_shCont_WTAP_responsive.bin HepG2_shSetD2_WTAP_responsive.bin | cut -f 1,2,3,6 > HepG2_chip_WTAP_responsive.bin
sed  -i '1i region\tbin\tshCont\tshSetD2' HepG2_chip_WTAP_responsive.bin
rm -f HepG2_shCont_WTAP_responsive.bin HepG2_shSetD2_WTAP_responsive.bin

bedBinDistribution.pl --feature $binFeature -span 5 -smooth move -t count --input HepG2_chip_shCont_nonWriter_responsive.bed -bed6 $binAnno -o ./HepG2_chip_shCont_nonWriter_responsive.bin
bedBinDistribution.pl --feature $binFeature -span 5 -smooth move -t count --input HepG2_chip_shSetD2_nonWriter_responsive.bed -bed6 $binAnno -o ./HepG2_chip_shSetD2_nonWriter_responsive.bin
paste HepG2_chip_shCont_nonWriter_responsive.bin HepG2_chip_shSetD2_nonWriter_responsive.bin | cut -f 1,2,3,6 > HepG2_chip_nonWriter_responsive.bin
sed  -i '1i region\tbin\tshCont\tshSetD2' HepG2_chip_nonWriter_responsive.bin
rm -f HepG2_chip_shCont_nonWriter_responsive.bin HepG2_chip_shSetD2_nonWriter_responsive.bin

bedBinDistribution.pl --feature $binFeature -span 5 -smooth move -t count --input HepG2_macs_shCont_peaks.bed -bed6 $binAnno -o ./HepG2_chip_shCont.bin
bedBinDistribution.pl --feature $binFeature -span 5 -smooth move -t count --input HepG2_macs_shSetD2_peaks.bed -bed6 $binAnno -o ./HepG2_chip_shSetD2.bin
bedBinDistribution.pl --feature $binFeature -span 5 -smooth move -t count --input HepG2_macs_shM14_peaks.bed -bed6 $binAnno -o ./HepG2_chip_shM14.bin
bedBinDistribution.pl --feature $binFeature -span 5 -smooth move -t count --input HepG2_macs_shM3_peaks.bed -bed6 $binAnno -o ./HepG2_chip_shM3.bin
bedBinDistribution.pl --feature $binFeature -span 5 -smooth move -t count --input HepG2_macs_shWTAP_peaks.bed -bed6 $binAnno -o ./HepG2_chip_shWTAP.bin
paste HepG2_chip_shCont.bin HepG2_chip_shSetD2.bin HepG2_chip_shM14.bin HepG2_chip_shM3.bin HepG2_chip_shWTAP.bin | cut -f 1,2,3,6,9,12,15 > HepG2_chip_all.bin
sed  -i '1i region\tbin\tshCont\tshSetD2\tshM14\tshM3\tshWTAP' HepG2_chip_all.bin
rm -f HepG2_chip_shCont.bin HepG2_chip_shSetD2.bin HepG2_chip_shM14.bin HepG2_chip_shM3.bin HepG2_chip_shWTAP.bin

#### gene type counts and region counts of SetD2 responsive H3K36me3 sites

regionDistribution.pl --feature $regionFeature --input HepG2_shCont_SetD2_responsive.bed -bed6 $regionAnno -o ./HepG2_macs_shCont_SetD2_responsive.region
sed  -i '1i region\tpeakNumber\tenrichment' ./HepG2_macs_shCont_SetD2_responsive.region
