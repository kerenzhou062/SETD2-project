#!/bin/sh

annotation=/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.mRNA.longest.exon.bed6
histonePeak=/data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/bed6/HepG2_macs_shCont_peaks.bed

### m6a-H3K36me3
rm -rf /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/m6a-H3K36me3/
mkdir /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/m6a-H3K36me3/
cd /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/m6a-H3K36me3/
### control and hypomethylated m6A
cp ./../bed6/HepG2_m6a_shCont_SetD2_dependent.bed6 HepG2_m6a_shCont_SetD2_dependent.bed6
bed2exclude.pl -s -a ./../bed6/HepG2_shCont_m6a.bed6 -b HepG2_m6a_shCont_SetD2_dependent.bed6 -o HepG2_m6a_shCont_Non_SetD2_dependent.bed6

bed2overlap.pl -aOver -a HepG2_m6a_shCont_SetD2_dependent.bed6 -b $histonePeak -o ./HepG2_m6a_shCont_SetD2_dependent+H3K36me3.bed6
bed2overlap.pl -aOver -a ./HepG2_m6a_shCont_Non_SetD2_dependent.bed6 -b $histonePeak -o ./HepG2_m6a_shCont_Non_SetD2_dependent+H3K36me3.bed6

bedBinDistribution.pl -strand -smooth move -t count --input HepG2_m6a_shCont_SetD2_dependent+H3K36me3.bed6 -bed6 $annotation -o ./HepG2_m6a_shCont_SetD2_dependent+H3K36me3.bed6bin
bedBinDistribution.pl -strand -smooth move -t count --input HepG2_m6a_shCont_Non_SetD2_dependent+H3K36me3.bed6 -bed6 $annotation -o ./HepG2_m6a_shCont_Non_SetD2_dependent+H3K36me3.bed6bin

paste HepG2_m6a_shCont_SetD2_dependent+H3K36me3.bed6bin HepG2_m6a_shCont_Non_SetD2_dependent+H3K36me3.bed6bin | cut -f 1,2,3,6 > HepG2_m6a+H3K36me3.bed6bin
sed  -i '1i region\tbin\tSETD2 Dependent m6A(H3K36me3+)\tNon-SETD2 Dependent m6A(H3K36me3+)' HepG2_m6a+H3K36me3.bed6bin
rm -f HepG2_m6a_shCont_SetD2_dependent+H3K36me3.bed6bin HepG2_m6a_shCont_Non_SetD2_dependent+H3K36me3.bed6bin

bed2exclude.pl -a HepG2_m6a_shCont_SetD2_dependent.bed6 -b $histonePeak -o ./HepG2_m6a_shCont_SetD2_dependent-H3K36me3.bed6
bed2exclude.pl -a HepG2_m6a_shCont_Non_SetD2_dependent.bed6 -b $histonePeak -o ./HepG2_m6a_shCont_Non_SetD2_dependent-H3K36me3.bed6

bedBinDistribution.pl -strand -smooth move -t count --input HepG2_m6a_shCont_SetD2_dependent-H3K36me3.bed6 -bed6 $annotation -o ./HepG2_m6a_shCont_SetD2_dependent-H3K36me3.bed6bin
bedBinDistribution.pl -strand -smooth move -t count --input HepG2_m6a_shCont_Non_SetD2_dependent-H3K36me3.bed6 -bed6 $annotation -o ./HepG2_m6a_shCont_Non_SetD2_dependent-H3K36me3.bed6bin

paste HepG2_m6a_shCont_SetD2_dependent-H3K36me3.bed6bin HepG2_m6a_shCont_Non_SetD2_dependent-H3K36me3.bed6bin | cut -f 1,2,3,6 > HepG2_m6a-H3K36me3.bed6bin
sed  -i '1i region\tbin\tSETD2 Dependent m6A(H3K36me3-)\tNon-SETD2 Dependent m6A(H3K36me3-)' HepG2_m6a-H3K36me3.bed6bin
rm -f HepG2_m6a_shCont_SetD2_dependent+H3K36me3.bed6bin HepG2_m6a_shCont_Non_SetD2_dependen-H3K36me3.bed6bin

paste HepG2_m6a+H3K36me3.bed6bin HepG2_m6a-H3K36me3.bed6bin | cut -f 1,2,3,4,7,8 > HepG2_m6a_H3K36me3.bed6bin
rm -f  HepG2_m6a+H3K36me3.bed6bin HepG2_m6a-H3K36me3.bed6bin

setD2_dependent_m6a_plus_H3K36me3_peakNum=`awk 'END{print NR}' HepG2_m6a_shCont_SetD2_dependent+H3K36me3.bed6`
setD2_dependent_m6a_minus_H3K36me3_peakNum=`awk 'END{print NR}' HepG2_m6a_shCont_SetD2_dependent-H3K36me3.bed6`
non_setD2_dependent_m6a_plus_H3K36me3_peakNum=`awk 'END{print NR}' HepG2_m6a_shCont_Non_SetD2_dependent+H3K36me3.bed6`
non_setD2_dependent_m6a_minus_H3K36me3_peakNum=`awk 'END{print NR}' HepG2_m6a_shCont_Non_SetD2_dependent-H3K36me3.bed6`

touch HepG2_shCont_m6a_H3K36me3_peakNum.txt
echo -e "type\tH3K36me3+\tH3K36me3-" > HepG2_shCont_m6a_H3K36me3_peakNum.txt
echo -e "SetD2 Dependent\t${setD2_dependent_m6a_plus_H3K36me3_peakNum}\t${setD2_dependent_m6a_minus_H3K36me3_peakNum}" >> HepG2_shCont_m6a_H3K36me3_peakNum.txt
echo -e "Non-SetD2 Dependent\t${non_setD2_dependent_m6a_plus_H3K36me3_peakNum}\t${non_setD2_dependent_m6a_minus_H3K36me3_peakNum}" >> HepG2_shCont_m6a_H3K36me3_peakNum.txt
