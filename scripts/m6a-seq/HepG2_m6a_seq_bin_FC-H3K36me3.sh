#!/bin/sh

annotation=/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.mRNA.longest.exon.bed6
histonePeak=/data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/bed6/HepG2_macs_shCont_peaks.bed

### m6a-H3K36me3
rm -rf /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/m6a-H3K36me3/
mkdir /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/m6a-H3K36me3/
cd /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/m6a-H3K36me3/
### control and hypomethylated m6A
awk '{if(FNR>1){if($14<0.5) {for(i=1;i<12;i++) {printf("%s\t",$i);} printf("%s", $12); printf "\n";} } }' ./../xls/allPeak/HepG2_m6a_peak_-F0.5_FC.txt > ./HepG2_m6a_shCont_SetD2_dependent.bed6


bed2overlap.pl -aOver -a ./../bed6/HepG2_shCont_m6a.bed6 -b $histonePeak -o ./HepG2_shCont_m6a+H3K36me3.bed6
bed2overlap.pl -aOver -a HepG2_m6a_shCont_SetD2_dependent.bed6 -b $histonePeak -o ./HepG2_SetD2_dependent_m6a+H3K36me3.bed6
#exomePeakToSummit.pl -p 0.05 --name shCont_m6a+H3K36me3 -b HepG2_shCont_m6a+H3K36me3.bed6 -o ./HepG2_shCont_m6a+H3K36me3.bed6
#exomePeakToSummit.pl -p 0.05 --name SetD2_dependent_m6a+H3K36me3 -b HepG2_SetD2_dependent_m6a+H3K36me3.bed6 -o ./HepG2_SetD2_dependent_m6a+H3K36me3.bed6

bedBinDistribution.pl -strand -smooth move -t percentage --input HepG2_shCont_m6a+H3K36me3.bed6 -bed6 $annotation -o ./HepG2_shCont_m6a+H3K36me3.bed6bin
bedBinDistribution.pl -strand -smooth move -t percentage --input HepG2_SetD2_dependent_m6a+H3K36me3.bed6 -bed6 $annotation -o ./HepG2_SetD2_dependent_m6a+H3K36me3.bed6bin

paste HepG2_shCont_m6a+H3K36me3.bed6bin HepG2_SetD2_dependent_m6a+H3K36me3.bed6bin | cut -f 1,2,3,6 > HepG2_m6a+H3K36me3.bed6bin
sed  -i '1i region\tbin\tcontrol m6A(H3K36me3+)\thypomethylated m6A(H3K36me3+)' HepG2_m6a+H3K36me3.bed6bin
rm -f HepG2_shCont_m6a+H3K36me3.bed6bin HepG2_SetD2_dependent_m6a+H3K36me3.bed6bin

bed2exclude.pl -a ./../bed6/HepG2_shCont_m6a.bed6 -b $histonePeak -o ./HepG2_shCont_m6a-H3K36me3.bed6
bed2exclude.pl -a HepG2_m6a_shCont_SetD2_dependent.bed6 -b $histonePeak -o ./HepG2_SetD2_dependent_m6a-H3K36me3.bed6
#exomePeakToSummit.pl -p 0.05 --name shCont_m6a-H3K36me3 -b HepG2_shCont_m6a-H3K36me3.bed6 -o ./HepG2_shCont_m6a-H3K36me3.bed6
#exomePeakToSummit.pl -p 0.05 --name SetD2_dependent_m6a-H3K36me3 -b HepG2_SetD2_dependent_m6a-H3K36me3.bed6 -o ./HepG2_SetD2_dependent_m6a-H3K36me3.bed6

bedBinDistribution.pl -strand -smooth move -t percentage --input HepG2_shCont_m6a-H3K36me3.bed6 -bed6 $annotation -o ./HepG2_shCont_m6a-H3K36me3.bed6bin
bedBinDistribution.pl -strand -smooth move -t percentage --input HepG2_SetD2_dependent_m6a-H3K36me3.bed6 -bed6 $annotation -o ./HepG2_SetD2_dependent_m6a-H3K36me3.bed6bin

paste HepG2_shCont_m6a-H3K36me3.bed6bin HepG2_SetD2_dependent_m6a-H3K36me3.bed6bin | cut -f 1,2,3,6 > HepG2_m6a-H3K36me3.bed6bin
sed -i '1i region\tbin\tcontrol m6A(H3K36me3-)\thypomethylated m6A(H3K36me3-)' HepG2_m6a-H3K36me3.bed6bin
rm -f  HepG2_shCont_m6a-H3K36me3.bed6bin HepG2_SetD2_dependent_m6a-H3K36me3.bed6bin

paste HepG2_m6a+H3K36me3.bed6bin HepG2_m6a-H3K36me3.bed6bin | cut -f 1,2,3,4,7,8 > HepG2_SetD2_m6a_H3K36me3.bed6bin
rm -f  HepG2_m6a+H3K36me3.bed6bin HepG2_m6a-H3K36me3.bed6bin

shCont_m6a_plus_H3K36me3_peakNum=`awk 'END{print NR}' HepG2_shCont_m6a+H3K36me3.bed6`
shCont_m6a_minus_H3K36me3_peakNum=`awk 'END{print NR}' HepG2_shCont_m6a-H3K36me3.bed6`
hypomethylated_m6a_plus_H3K36me3_peakNum=`awk 'END{print NR}' HepG2_SetD2_dependent_m6a+H3K36me3.bed6`
hypomethylated_m6a_minus_H3K36me3_peakNum=`awk 'END{print NR}' HepG2_SetD2_dependent_m6a-H3K36me3.bed6`

touch HepG2_shCont_SetD2_m6a_H3K36me3_peakNum.txt
echo -e "type\tH3K36me3+\tH3K36me3-" > HepG2_shCont_SetD2_m6a_H3K36me3_peakNum.txt
echo -e "control\t${shCont_m6a_plus_H3K36me3_peakNum}\t${shCont_m6a_minus_H3K36me3_peakNum}" >> HepG2_shCont_SetD2_m6a_H3K36me3_peakNum.txt
echo -e "hypomethylated\t${hypomethylated_m6a_plus_H3K36me3_peakNum}\t${hypomethylated_m6a_minus_H3K36me3_peakNum}" >> HepG2_shCont_SetD2_m6a_H3K36me3_peakNum.txt


### Control and shSetD2

bed2overlap.pl -aOver -a ./../bed6/HepG2_shCont_m6a.bed6 -b $histonePeak -o ./HepG2_shCont_m6a+H3K36me3.bed6
bed2overlap.pl -aOver -a ./../bed6/HepG2_shSetD2_m6a.bed6 -b $histonePeak -o ./HepG2_shSetD2_m6a+H3K36me3.bed6
#exomePeakToSummit.pl -p 0.05 --name shCont_m6a+H3K36me3 -b HepG2_shCont_m6a+H3K36me3.bed6 -o ./HepG2_shCont_m6a+H3K36me3.bed6
#exomePeakToSummit.pl -p 0.05 --name shSetD2_m6a+H3K36me3 -b HepG2_shSetD2_m6a+H3K36me3.bed6 -o ./HepG2_shSetD2_m6a+H3K36me3.bed6

bedBinDistribution.pl -strand -smooth move -t percentage --input HepG2_shCont_m6a+H3K36me3.bed6 -bed6 $annotation -o ./HepG2_shCont_m6a+H3K36me3.bed6bin
bedBinDistribution.pl -strand -smooth move -t percentage --input HepG2_shSetD2_m6a+H3K36me3.bed6 -bed6 $annotation -o ./HepG2_shSetD2_m6a+H3K36me3.bed6bin

paste HepG2_shCont_m6a+H3K36me3.bed6bin HepG2_shSetD2_m6a+H3K36me3.bed6bin | cut -f 1,2,3,6 > HepG2_m6a+H3K36me3.bed6bin
sed  -i '1i region\tbin\tshCont(H3K36me3+)\tshSetD2(H3K36me3+)' HepG2_m6a+H3K36me3.bed6bin
rm -f HepG2_shCont_m6a+H3K36me3.bed6bin HepG2_shSetD2_m6a+H3K36me3.bed6bin

bed2exclude.pl -a ./../bed6/HepG2_shCont_m6a.bed6 -b $histonePeak -o ./HepG2_shCont_m6a-H3K36me3.bed6
bed2exclude.pl -a ./../bed6/HepG2_shSetD2_m6a.bed6 -b $histonePeak -o ./HepG2_shSetD2_m6a-H3K36me3.bed6
#exomePeakToSummit.pl -p 0.05 --name shCont_m6a-H3K36me3 -b HepG2_shCont_m6a-H3K36me3.bed6 -o ./HepG2_shCont_m6a-H3K36me3.bed6
#exomePeakToSummit.pl -p 0.05 --name shSetD2_m6a-H3K36me3 -b HepG2_shSetD2_m6a-H3K36me3.bed6 -o ./HepG2_shSetD2_m6a-H3K36me3.bed6

bedBinDistribution.pl -strand -smooth move -t percentage --input HepG2_shCont_m6a-H3K36me3.bed6 -bed6 $annotation -o ./HepG2_shCont_m6a-H3K36me3.bed6bin
bedBinDistribution.pl -strand -smooth move -t percentage --input HepG2_shSetD2_m6a-H3K36me3.bed6 -bed6 $annotation -o ./HepG2_shSetD2_m6a-H3K36me3.bed6bin

paste HepG2_shCont_m6a-H3K36me3.bed6bin HepG2_shSetD2_m6a-H3K36me3.bed6bin | cut -f 1,2,3,6 > HepG2_m6a-H3K36me3.bed6bin
sed  -i '1i region\tbin\tshCont(H3K36me3-)\tshSetD2(H3K36me3-)' HepG2_m6a-H3K36me3.bed6bin
rm -f  HepG2_shCont_m6a-H3K36me3.bed6bin HepG2_shSetD2_m6a-H3K36me3.bed6bin

paste HepG2_m6a+H3K36me3.bed6bin HepG2_m6a-H3K36me3.bed6bin | cut -f 1,2,3,4,7,8 > HepG2_m6a_H3K36me3.bed6bin
rm -f HepG2_m6a+H3K36me3.bed6bin HepG2_m6a-H3K36me3.bed6bin

shCont_m6a_plus_H3K36me3_peakNum=`awk 'END{print NR}' HepG2_shCont_m6a+H3K36me3.bed6`
shCont_m6a_minus_H3K36me3_peakNum=`awk 'END{print NR}' HepG2_shCont_m6a-H3K36me3.bed6`
shSetD2_m6a_plus_H3K36me3_peakNum=`awk 'END{print NR}' HepG2_shSetD2_m6a+H3K36me3.bed6`
shSetD2_m6a_minus_H3K36me3_peakNum=`awk 'END{print NR}' HepG2_shSetD2_m6a-H3K36me3.bed6`

touch HepG2_m6a_H3K36me3_peakNum.txt
echo -e "type\tH3K36me3+\tH3K36me3-" > HepG2_m6a_H3K36me3_peakNum.txt
echo -e "shCont\t${shCont_m6a_plus_H3K36me3_peakNum}\t${shCont_m6a_minus_H3K36me3_peakNum}" >> HepG2_m6a_H3K36me3_peakNum.txt
echo -e "shSetD2\t${shSetD2_m6a_plus_H3K36me3_peakNum}\t${shSetD2_m6a_minus_H3K36me3_peakNum}" >> HepG2_m6a_H3K36me3_peakNum.txt
