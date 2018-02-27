#!/bin/sh
mRNAAnnotation=/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.mRNA.longest.exon.bed6
allAnnotation=/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.all.longest.exon.bed6
histonePeak=/data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/bed6/HepG2_macs_shCont_peaks.bed
bed12Path=/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/bed12
bed6Path=/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/bed6

cd /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq
### bed12
rm -rf $bed12Path
mkdir $bed12Path
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shCont/shCont_m6A_rep1-2-3/peak.xls -o $bed12Path/HepG2_shCont_m6a.bed12
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shSetD2/shSetD2_m6A_rep1-2-3/peak.xls -o $bed12Path/HepG2_shSetD2_m6a.bed12
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shM14/shM14_m6A_rep1-2-3/peak.xls -o $bed12Path/HepG2_shM14_m6a.bed12
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shM3/shM3_m6A_rep1-2-3/peak.xls -o $bed12Path/HepG2_shM3_m6a.bed12
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shWTAP/shWTAP_m6A_rep1-2-3/peak.xls -o $bed12Path/HepG2_shWTAP_m6a.bed12
bedtools intersect -a $bed12Path/HepG2_shCont_m6a.bed12 -b $histonePeak -wa | sort -t $'\t' -k 1,1 -k 2,2n | uniq > $bed12Path/HepG2_shCont_m6a_H3K36me3.bed12
bedtools intersect -a $bed12Path/HepG2_shSetD2_m6a.bed12 -b $histonePeak -wa | sort -t $'\t' -k 1,1 -k 2,2n | uniq > $bed12Path/HepG2_shSetD2_m6a_H3K36me3.bed12
bedtools intersect -a $bed12Path/HepG2_shM14_m6a.bed12 -b $histonePeak -wa | sort -t $'\t' -k 1,1 -k 2,2n | uniq > $bed12Path/HepG2_shM14_m6a_H3K36me3.bed12
bedtools intersect -a $bed12Path/HepG2_shM3_m6a.bed12 -b $histonePeak -wa | sort -t $'\t' -k 1,1 -k 2,2n | uniq > $bed12Path/HepG2_shM3_m6a_H3K36me3.bed12
bedtools intersect -a $bed12Path/HepG2_shWTAP_m6a.bed12 -b $histonePeak -wa | sort -t $'\t' -k 1,1 -k 2,2n | uniq > $bed12Path/HepG2_shWTAP_m6a_H3K36me3.bed12

### bed6
rm -rf $bed6Path
mkdir $bed6Path
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -x shCont/shCont_m6A_rep1-2-3/peak.xls -o $bed6Path/HepG2_shCont_m6a.bed6
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -x shSetD2/shSetD2_m6A_rep1-2-3/peak.xls -o $bed6Path/HepG2_shSetD2_m6a.bed6
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -x shM14/shM14_m6A_rep1-2-3/peak.xls -o $bed6Path/HepG2_shM14_m6a.bed6
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -x shM3/shM3_m6A_rep1-2-3/peak.xls -o $bed6Path/HepG2_shM3_m6a.bed6
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -x shWTAP/shWTAP_m6A_rep1-2-3/peak.xls -o $bed6Path/HepG2_shWTAP_m6a.bed6
bedtools intersect -a $bed6Path/HepG2_shCont_m6a.bed6 -b $histonePeak -wa | sort -t $'\t' -k 1,1 -k 2,2n | uniq > $bed6Path/HepG2_shCont_m6a_H3K36me3.bed6
bedtools intersect -a $bed6Path/HepG2_shSetD2_m6a.bed6 -b $histonePeak -wa | sort -t $'\t' -k 1,1 -k 2,2n | uniq > $bed6Path/HepG2_shSetD2_m6a_H3K36me3.bed6
bedtools intersect -a $bed6Path/HepG2_shM14_m6a.bed6 -b $histonePeak -wa | sort -t $'\t' -k 1,1 -k 2,2n | uniq > $bed6Path/HepG2_shM14_m6a_H3K36me3.bed6
bedtools intersect -a $bed6Path/HepG2_shM3_m6a.bed6 -b $histonePeak -wa | sort -t $'\t' -k 1,1 -k 2,2n | uniq > $bed6Path/HepG2_shM3_m6a_H3K36me3.bed6
bedtools intersect -a $bed6Path/HepG2_shWTAP_m6a.bed6 -b $histonePeak -wa | sort -t $'\t' -k 1,1 -k 2,2n | uniq > $bed6Path/HepG2_shWTAP_m6a_H3K36me3.bed6


### identify writer targets
cd $bed12Path

#### SetD2 responsive m6a sites
awk '{if(FNR>1){if($14<(1/2)) {for(i=1;i<12;i++) {printf("%s\t",$i);} printf("%s", $12); printf "\n";}}}' ./../xls/allPeak/HepG2_m6a_peak_FC.txt > ./HepG2_m6a_shCont_SetD2_responsive.bed12
exomePeakToSummit.pl --name SetD2_responsive_m6a -b HepG2_m6a_shCont_SetD2_responsive.bed12 -o $bed6Path/HepG2_m6a_shCont_SetD2_responsive.bed6
contrlPeak=`cat HepG2_shCont_m6a.bed12 | wc -l`
setd2Peak=`cat HepG2_shSetD2_m6a.bed12 | wc -l`

awk '{if(FNR>1){if($15<(1/2)) {for(i=1;i<12;i++) {printf("%s\t",$i);} printf("%s", $12); printf "\n";}}}' ./../xls/allPeak/HepG2_m6a_peak_FC.txt > HepG2_m6a_shCont_M14_responsive.bed12
bedtools intersect -a HepG2_shSetD2_m6a.bed12 -b HepG2_m6a_shCont_M14_responsive.bed12 -s -split -u | sort  > ./temp.bed12
bedtools intersect -a temp.bed12 -b HepG2_m6a_shCont_SetD2_responsive.bed12 -s -split -v | sort | uniq > HepG2_m6a_shSetD2_M14_responsive.bed12

bedtools intersect -a HepG2_m6a_shCont_M14_responsive.bed12 -b $histonePeak -wa | sort -t $'\t' -k 1,1 -k 2,2n | uniq > HepG2_m6a_shCont_M14_responsive_H3K36me3.bed12
bedtools intersect -a HepG2_m6a_shSetD2_M14_responsive.bed12 -b $histonePeak -wa | sort -t $'\t' -k 1,1 -k 2,2n | uniq > HepG2_m6a_shSetD2_M14_responsive_H3K36me3.bed12
exomePeakToSummit.pl --name shCont_M14_m6a -b HepG2_m6a_shCont_M14_responsive.bed12 -o $bed6Path/HepG2_m6a_shCont_M14_responsive.bed6
exomePeakToSummit.pl --name shSetD2_M14_m6a -b HepG2_m6a_shSetD2_M14_responsive.bed12 -o $bed6Path/HepG2_m6a_shSetD2_M14_responsive.bed6
exomePeakToSummit.pl --name shCont_M14_m6a -b HepG2_m6a_shCont_M14_responsive_H3K36me3.bed12 -o $bed6Path/HepG2_m6a_shCont_M14_responsive_H3K36me3.bed6
exomePeakToSummit.pl --name shSetD2_M14_m6a -b HepG2_m6a_shSetD2_M14_responsive_H3K36me3.bed12 -o $bed6Path/HepG2_m6a_shSetD2_M14_responsive_H3K36me3.bed6


awk '{if(FNR>1){if($16<(1/2)) {for(i=1;i<12;i++) {printf("%s\t",$i);} printf("%s", $12); printf "\n";}}}' ./../xls/allPeak/HepG2_m6a_peak_FC.txt > HepG2_m6a_shCont_M3_responsive.bed12
bedtools intersect -a HepG2_shSetD2_m6a.bed12 -b HepG2_m6a_shCont_M3_responsive.bed12 -s -split -u | sort  > ./temp.bed12
bedtools intersect -a temp.bed12 -b HepG2_m6a_shCont_SetD2_responsive.bed12 -s -split -v | sort | uniq > HepG2_m6a_shSetD2_M3_responsive.bed12

bedtools intersect -a HepG2_m6a_shCont_M3_responsive.bed12 -b $histonePeak -wa | sort -t $'\t' -k 1,1 -k 2,2n | uniq > HepG2_m6a_shCont_M3_responsive_H3K36me3.bed12
bedtools intersect -a HepG2_m6a_shSetD2_M3_responsive.bed12 -b $histonePeak -wa | sort -t $'\t' -k 1,1 -k 2,2n | uniq > HepG2_m6a_shSetD2_M3_responsive_H3K36me3.bed12
exomePeakToSummit.pl --name shCont_M3_m6a -b HepG2_m6a_shCont_M3_responsive.bed12 -o $bed6Path/HepG2_m6a_shCont_M3_responsive.bed6
exomePeakToSummit.pl --name shSetD2_M3_m6a -b HepG2_m6a_shSetD2_M3_responsive.bed12 -o $bed6Path/HepG2_m6a_shSetD2_M3_responsive.bed6
exomePeakToSummit.pl --name shCont_M3_m6a -b HepG2_m6a_shCont_M3_responsive_H3K36me3.bed12 -o $bed6Path/HepG2_m6a_shCont_M3_responsive_H3K36me3.bed6
exomePeakToSummit.pl --name shSetD2_M3_m6a -b HepG2_m6a_shSetD2_M3_responsive_H3K36me3.bed12 -o $bed6Path/HepG2_m6a_shSetD2_M3_responsive_H3K36me3.bed6

awk '{if(FNR>1){if($17<(1/2)) {for(i=1;i<12;i++) {printf("%s\t",$i);} printf("%s", $12); printf "\n";}}}' ./../xls/allPeak/HepG2_m6a_peak_FC.txt > HepG2_m6a_shCont_WTAP_responsive.bed12
bedtools intersect -a HepG2_shSetD2_m6a.bed12 -b HepG2_m6a_shCont_WTAP_responsive.bed12 -s -split -u | sort  > ./temp.bed12
bedtools intersect -a temp.bed12 -b HepG2_m6a_shCont_SetD2_responsive.bed12 -s -split -v | sort | uniq > HepG2_m6a_shSetD2_WTAP_responsive.bed12

bedtools intersect -a HepG2_m6a_shCont_WTAP_responsive.bed12 -b $histonePeak -wa | sort -t $'\t' -k 1,1 -k 2,2n | uniq > HepG2_m6a_shCont_WTAP_responsive_H3K36me3.bed12
bedtools intersect -a HepG2_m6a_shSetD2_WTAP_responsive.bed12 -b $histonePeak -wa | sort -t $'\t' -k 1,1 -k 2,2n | uniq > HepG2_m6a_shSetD2_WTAP_responsive_H3K36me3.bed12
exomePeakToSummit.pl --name shCont_WTAP_m6a -b HepG2_m6a_shCont_WTAP_responsive.bed12 -o $bed6Path/HepG2_m6a_shCont_WTAP_responsive.bed6
exomePeakToSummit.pl --name shSetD2_WTAP_m6a -b HepG2_m6a_shSetD2_WTAP_responsive.bed12 -o $bed6Path/HepG2_m6a_shSetD2_WTAP_responsive.bed6
exomePeakToSummit.pl --name shCont_WTAP_m6a -b HepG2_m6a_shCont_WTAP_responsive_H3K36me3.bed12 -o $bed6Path/HepG2_m6a_shCont_WTAP_responsive_H3K36me3.bed6
exomePeakToSummit.pl --name shSetD2_WTAP_m6a -b HepG2_m6a_shSetD2_WTAP_responsive_H3K36me3.bed12 -o $bed6Path/HepG2_m6a_shSetD2_WTAP_responsive_H3K36me3.bed6

awk '{if(FNR>1){if($15<(1/2) && $16<(1/2) && $17<(1/2)) {for(i=1;i<12;i++) {printf("%s\t",$i);} printf("%s", $12); printf "\n";}}}' ./../xls/allPeak/HepG2_m6a_peak_FC.txt > HepG2_m6a_shCont_writerShare_responsive.bed12
bedtools intersect -a HepG2_shSetD2_m6a.bed12 -b HepG2_m6a_shCont_writerShare_responsive.bed12 -s -split -u | sort  > ./temp.bed12
bedtools intersect -a temp.bed12 -b HepG2_m6a_shCont_SetD2_responsive.bed12 -s -split -v | sort | uniq > HepG2_m6a_shSetD2_writerShare_responsive.bed12

bedtools intersect -a HepG2_m6a_shCont_writerShare_responsive.bed12 -b $histonePeak -wa | sort -t $'\t' -k 1,1 -k 2,2n | uniq > HepG2_m6a_shCont_writerShare_responsive_H3K36me3.bed12
bedtools intersect -a HepG2_m6a_shSetD2_writerShare_responsive.bed12 -b $histonePeak -wa | sort -t $'\t' -k 1,1 -k 2,2n | uniq > HepG2_m6a_shSetD2_writerShare_responsive_H3K36me3.bed12
exomePeakToSummit.pl --name shCont_writerShare_m6a -b HepG2_m6a_shCont_writerShare_responsive.bed12 -o $bed6Path/HepG2_m6a_shCont_writerShare_responsive.bed6
exomePeakToSummit.pl --name shSetD2_writerShare_m6a -b HepG2_m6a_shSetD2_writerShare_responsive.bed12 -o $bed6Path/HepG2_m6a_shSetD2_writerShare_responsive.bed6
exomePeakToSummit.pl --name shCont_writerShare_m6a -b HepG2_m6a_shCont_writerShare_responsive_H3K36me3.bed12 -o $bed6Path/HepG2_m6a_shCont_writerShare_responsive_H3K36me3.bed6
exomePeakToSummit.pl --name shSetD2_writerShare_m6a -b HepG2_m6a_shSetD2_writerShare_responsive_H3K36me3.bed12 -o $bed6Path/HepG2_m6a_shSetD2_writerShare_responsive_H3K36me3.bed6

awk '{if(FNR>1){if($15<(1/2) || $16<(1/2) || $17<(1/2)) {for(i=1;i<12;i++) {printf("%s\t",$i);} printf("%s", $12); printf "\n";}}}' ./../xls/allPeak/HepG2_m6a_peak_FC.txt > HepG2_m6a_shCont_writerAll_responsive.bed12
bedtools intersect -a HepG2_shSetD2_m6a.bed12 -b HepG2_m6a_shCont_writerAll_responsive.bed12 -s -split -u | sort  > ./temp.bed12
bedtools intersect -a temp.bed12 -b HepG2_m6a_shCont_SetD2_responsive.bed12 -s -split -v | sort | uniq > HepG2_m6a_shSetD2_writerAll_responsive.bed12

bedtools intersect -a HepG2_m6a_shCont_writerAll_responsive.bed12 -b $histonePeak -wa | sort -t $'\t' -k 1,1 -k 2,2n | uniq > HepG2_m6a_shCont_writerAll_responsive_H3K36me3.bed12
bedtools intersect -a HepG2_m6a_shSetD2_writerAll_responsive.bed12 -b $histonePeak -wa | sort -t $'\t' -k 1,1 -k 2,2n | uniq > HepG2_m6a_shSetD2_writerAll_responsive_H3K36me3.bed12
exomePeakToSummit.pl --name shCont_writerAll_m6a -b HepG2_m6a_shCont_writerAll_responsive.bed12 -o $bed6Path/HepG2_m6a_shCont_writerAll_responsive.bed6
exomePeakToSummit.pl --name shSetD2_writerAll_m6a -b HepG2_m6a_shSetD2_writerAll_responsive.bed12 -o $bed6Path/HepG2_m6a_shSetD2_writerAll_responsive.bed6
exomePeakToSummit.pl --name shCont_writerAll_m6a -b HepG2_m6a_shCont_writerAll_responsive_H3K36me3.bed12 -o $bed6Path/HepG2_m6a_shCont_writerAll_responsive_H3K36me3.bed6
exomePeakToSummit.pl --name shSetD2_writerAll_m6a -b HepG2_m6a_shSetD2_writerAll_responsive_H3K36me3.bed12 -o $bed6Path/HepG2_m6a_shSetD2_writerAll_responsive_H3K36me3.bed6

### identify non writer targets
cat HepG2_m6a_shCont_M14_responsive.bed12 HepG2_m6a_shCont_M3_responsive.bed12 HepG2_m6a_shCont_WTAP_responsive.bed12 | sort -t $'\t' -k1,1V -k2,2n | uniq > HepG2_shCont_allWriter_responsive.bed12
awk '{if(FNR>1){if($15>=(1/2)&&$16>=(1/2)&&$17>=(1/2)) {for(i=1;i<12;i++) {printf("%s\t",$i);} printf("%s", $12); printf "\n";}}}' ./../xls/allPeak/HepG2_m6a_peak_FC.txt > HepG2_m6a_shCont_nonWriter_responsive.bed12
#bed2exclude.pl -split -s -a HepG2_shCont_m6a.bed12 -b HepG2_shCont_allWriter_responsive.bed12 -o ./HepG2_m6a_shCont_nonWriter_responsive.bed12
bedtools intersect -a HepG2_shSetD2_m6a.bed12 -b HepG2_m6a_shCont_nonWriter_responsive.bed12 -s -split -u | sort  > ./temp.bed12
bedtools intersect -a temp.bed12 -b HepG2_m6a_shCont_SetD2_responsive.bed12 -s -split -v | sort | uniq > HepG2_m6a_shSetD2_nonWriter_responsive.bed12

bedtools intersect -a HepG2_m6a_shCont_nonWriter_responsive.bed12 -b $histonePeak -wa | sort -t $'\t' -k 1,1 -k 2,2n | uniq > HepG2_m6a_shCont_nonWriter_responsive_H3K36me3.bed12
bedtools intersect -a HepG2_m6a_shSetD2_nonWriter_responsive.bed12 -b $histonePeak -wa | sort -t $'\t' -k 1,1 -k 2,2n | uniq > HepG2_m6a_shSetD2_nonWriter_responsive_H3K36me3.bed12
exomePeakToSummit.pl --name shCont_non_writer_m6a -b HepG2_m6a_shCont_nonWriter_responsive.bed12 -o $bed6Path/HepG2_m6a_shCont_nonWriter_responsive.bed6
exomePeakToSummit.pl --name shSetD2_non_writer_m6a -b HepG2_m6a_shSetD2_nonWriter_responsive.bed12 -o $bed6Path/HepG2_m6a_shSetD2_nonWriter_responsive.bed6
exomePeakToSummit.pl --name shCont_nonWriter_m6a -b HepG2_m6a_shCont_nonWriter_responsive_H3K36me3.bed12 -o $bed6Path/HepG2_m6a_shCont_nonWriter_responsive_H3K36me3.bed6
exomePeakToSummit.pl --name shSetD2_nonWriter_m6a -b HepG2_m6a_shSetD2_nonWriter_responsive_H3K36me3.bed12 -o $bed6Path/HepG2_m6a_shSetD2_nonWriter_responsive_H3K36me3.bed6

rm -f ./temp.bed12

###draw bin distribution based on count and bed12 for writer targets
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shCont_M14_responsive.bed12 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shCont_M14_responsive.bed12bin
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shSetD2_M14_responsive.bed12 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shSetD2_M14_responsive.bed12bin
paste HepG2_m6a_shCont_M14_responsive.bed12bin HepG2_m6a_shSetD2_M14_responsive.bed12bin | cut -f 1,2,3,6 > HepG2_m6a_M14_responsive.bed12bin
sed  -i '1i region\tbin\tshCont\tshSetD2' HepG2_m6a_M14_responsive.bed12bin
rm -f HepG2_m6a_shCont_M14_responsive.bed12bin HepG2_m6a_shSetD2_M14_responsive.bed12bin

bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shCont_M3_responsive.bed12 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shCont_M3_responsive.bed12bin
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shSetD2_M3_responsive.bed12 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shSetD2_M3_responsive.bed12bin
paste HepG2_m6a_shCont_M3_responsive.bed12bin HepG2_m6a_shSetD2_M3_responsive.bed12bin | cut -f 1,2,3,6 > HepG2_m6a_M3_responsive.bed12bin
sed  -i '1i region\tbin\tshCont\tshSetD2' HepG2_m6a_M3_responsive.bed12bin
rm -f HepG2_m6a_shCont_M3_responsive.bed12bin HepG2_m6a_shSetD2_M3_responsive.bed12bin

bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shCont_WTAP_responsive.bed12 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shCont_WTAP_responsive.bed12bin
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shSetD2_WTAP_responsive.bed12 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shSetD2_WTAP_responsive.bed12bin
paste HepG2_m6a_shCont_WTAP_responsive.bed12bin HepG2_m6a_shSetD2_WTAP_responsive.bed12bin | cut -f 1,2,3,6 > HepG2_m6a_WTAP_responsive.bed12bin
sed  -i '1i region\tbin\tshCont\tshSetD2' HepG2_m6a_WTAP_responsive.bed12bin
rm -f HepG2_m6a_shCont_WTAP_responsive.bed12bin HepG2_m6a_shSetD2_WTAP_responsive.bed12bin

bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shCont_writerShare_responsive.bed12 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shCont_writerShare_responsive.bed12bin
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shSetD2_writerShare_responsive.bed12 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shSetD2_writerShare_responsive.bed12bin
paste HepG2_m6a_shCont_writerShare_responsive.bed12bin HepG2_m6a_shSetD2_writerShare_responsive.bed12bin | cut -f 1,2,3,6 > HepG2_m6a_writerShare_responsive.bed12bin
sed  -i '1i region\tbin\tshCont\tshSetD2' HepG2_m6a_writerShare_responsive.bed12bin
rm -f HepG2_m6a_shCont_writerShare_responsive.bed12bin HepG2_m6a_shSetD2_writerShare_responsive.bed12bin

bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shCont_writerAll_responsive.bed12 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shCont_writerAll_responsive.bed12bin
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shSetD2_writerAll_responsive.bed12 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shSetD2_writerAll_responsive.bed12bin
paste HepG2_m6a_shCont_writerAll_responsive.bed12bin HepG2_m6a_shSetD2_writerAll_responsive.bed12bin | cut -f 1,2,3,6 > HepG2_m6a_writerAll_responsive.bed12bin
sed  -i '1i region\tbin\tshCont\tshSetD2' HepG2_m6a_writerAll_responsive.bed12bin
rm -f HepG2_m6a_shCont_writerAll_responsive.bed12bin HepG2_m6a_shSetD2_writerAll_responsive.bed12bin

bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shCont_nonWriter_responsive.bed12 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shCont_nonWriter_responsive.bed12bin
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shSetD2_nonWriter_responsive.bed12 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shSetD2_nonWriter_responsive.bed12bin
paste HepG2_m6a_shCont_nonWriter_responsive.bed12bin HepG2_m6a_shSetD2_nonWriter_responsive.bed12bin | cut -f 1,2,3,6 > HepG2_m6a_nonWriter_responsive.bed12bin
sed  -i '1i region\tbin\tshCont\tshSetD2' HepG2_m6a_nonWriter_responsive.bed12bin
rm -f HepG2_m6a_shCont_nonWriter_responsive.bed12bin HepG2_m6a_shSetD2_nonWriter_responsive.bed12bin

##

###draw bin distribution based on count and bed12, overlapped with H3K36me3 peaks
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shCont_M14_responsive_H3K36me3.bed12 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shCont_M14_responsive_H3K36me3.bed12bin
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shSetD2_M14_responsive_H3K36me3.bed12 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shSetD2_M14_responsive_H3K36me3.bed12bin
paste HepG2_m6a_shCont_M14_responsive_H3K36me3.bed12bin HepG2_m6a_shSetD2_M14_responsive_H3K36me3.bed12bin | cut -f 1,2,3,6 > HepG2_m6a_M14_responsive_H3K36me3.bed12bin
sed  -i '1i region\tbin\tshCont\tshSetD2' HepG2_m6a_M14_responsive_H3K36me3.bed12bin
rm -f HepG2_m6a_shCont_M14_responsive_H3K36me3.bed12bin HepG2_m6a_shSetD2_M14_responsive_H3K36me3.bed12bin

bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shCont_M3_responsive_H3K36me3.bed12 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shCont_M3_responsive_H3K36me3.bed12bin
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shSetD2_M3_responsive_H3K36me3.bed12 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shSetD2_M3_responsive_H3K36me3.bed12bin
paste HepG2_m6a_shCont_M3_responsive_H3K36me3.bed12bin HepG2_m6a_shSetD2_M3_responsive_H3K36me3.bed12bin | cut -f 1,2,3,6 > HepG2_m6a_M3_responsive_H3K36me3.bed12bin
sed  -i '1i region\tbin\tshCont\tshSetD2' HepG2_m6a_M3_responsive_H3K36me3.bed12bin
rm -f HepG2_m6a_shCont_M3_responsive_H3K36me3.bed12bin HepG2_m6a_shSetD2_M3_responsive_H3K36me3.bed12bin

bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shCont_WTAP_responsive_H3K36me3.bed12 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shCont_WTAP_responsive_H3K36me3.bed12bin
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shSetD2_WTAP_responsive_H3K36me3.bed12 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shSetD2_WTAP_responsive_H3K36me3.bed12bin
paste HepG2_m6a_shCont_WTAP_responsive_H3K36me3.bed12bin HepG2_m6a_shSetD2_WTAP_responsive_H3K36me3.bed12bin | cut -f 1,2,3,6 > HepG2_m6a_WTAP_responsive_H3K36me3.bed12bin
sed  -i '1i region\tbin\tshCont\tshSetD2' HepG2_m6a_WTAP_responsive_H3K36me3.bed12bin
rm -f HepG2_m6a_shCont_WTAP_responsive_H3K36me3.bed12bin HepG2_m6a_shSetD2_WTAP_responsive_H3K36me3.bed12bin

bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shCont_writerShare_responsive_H3K36me3.bed12 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shCont_writerShare_responsive_H3K36me3.bed12bin
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shSetD2_writerShare_responsive_H3K36me3.bed12 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shSetD2_writerShare_responsive_H3K36me3.bed12bin
paste HepG2_m6a_shCont_writerShare_responsive_H3K36me3.bed12bin HepG2_m6a_shSetD2_writerShare_responsive_H3K36me3.bed12bin | cut -f 1,2,3,6 > HepG2_m6a_writerShare_responsive_H3K36me3.bed12bin
sed  -i '1i region\tbin\tshCont\tshSetD2' HepG2_m6a_writerShare_responsive_H3K36me3.bed12bin
rm -f HepG2_m6a_shCont_writerShare_responsive_H3K36me3.bed12bin HepG2_m6a_shSetD2_writerShare_responsive_H3K36me3.bed12bin

bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shCont_writerAll_responsive_H3K36me3.bed12 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shCont_writerAll_responsive_H3K36me3.bed12bin
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shSetD2_writerAll_responsive_H3K36me3.bed12 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shSetD2_writerAll_responsive_H3K36me3.bed12bin
paste HepG2_m6a_shCont_writerAll_responsive_H3K36me3.bed12bin HepG2_m6a_shSetD2_writerAll_responsive_H3K36me3.bed12bin | cut -f 1,2,3,6 > HepG2_m6a_writerAll_responsive_H3K36me3.bed12bin
sed  -i '1i region\tbin\tshCont\tshSetD2' HepG2_m6a_writerAll_responsive_H3K36me3.bed12bin
rm -f HepG2_m6a_shCont_writerAll_responsive_H3K36me3.bed12bin HepG2_m6a_shSetD2_writerAll_responsive_H3K36me3.bed12bin

bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shCont_nonWriter_responsive_H3K36me3.bed12 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shCont_nonWriter_responsive_H3K36me3.bed12bin
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shSetD2_nonWriter_responsive_H3K36me3.bed12 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shSetD2_nonWriter_responsive_H3K36me3.bed12bin
paste HepG2_m6a_shCont_nonWriter_responsive_H3K36me3.bed12bin HepG2_m6a_shSetD2_nonWriter_responsive_H3K36me3.bed12bin | cut -f 1,2,3,6 > HepG2_m6a_nonWriter_responsive_H3K36me3.bed12bin
sed  -i '1i region\tbin\tshCont\tshSetD2' HepG2_m6a_nonWriter_responsive_H3K36me3.bed12bin
rm -f HepG2_m6a_shCont_nonWriter_responsive_H3K36me3.bed12bin HepG2_m6a_shSetD2_nonWriter_responsive_H3K36me3.bed12bin

### writer-responsive + setd2-responsive
cat HepG2_m6a_shCont_M14_responsive.bed12 HepG2_m6a_shCont_SetD2_responsive.bed12 | sort | uniq -d > HepG2_m6a_M14_SetD2_responsive.bed12
exomePeakToSummit.pl --name M14_SetD2_responsive_m6a -b HepG2_m6a_M14_SetD2_responsive.bed12 -o $bed6Path/HepG2_m6a_M14_SetD2_responsive.bed6
bedBinDistribution.pl -strand -smooth move -span 15 -t count --input HepG2_m6a_M14_SetD2_responsive.bed12 -bed6 $mRNAAnnotation -o ./HepG2_m6a_M14_SetD2_responsive.bed12bin

cat HepG2_m6a_shCont_M3_responsive.bed12 HepG2_m6a_shCont_SetD2_responsive.bed12 | sort | uniq -d > HepG2_m6a_M3_SetD2_responsive.bed12
exomePeakToSummit.pl --name M3_SetD2_responsive_m6a -b HepG2_m6a_M3_SetD2_responsive.bed12 -o $bed6Path/HepG2_m6a_M3_SetD2_responsive.bed6
bedBinDistribution.pl -strand -smooth move -span 15 -t count --input HepG2_m6a_M3_SetD2_responsive.bed12 -bed6 $mRNAAnnotation -o ./HepG2_m6a_M3_SetD2_responsive.bed12bin

cat HepG2_m6a_shCont_WTAP_responsive.bed12 HepG2_m6a_shCont_SetD2_responsive.bed12 | sort | uniq -d > HepG2_m6a_WTAP_SetD2_responsive.bed12
exomePeakToSummit.pl --name WTAP_SetD2_responsive_m6a -b HepG2_m6a_WTAP_SetD2_responsive.bed12 -o $bed6Path/HepG2_m6a_WTAP_SetD2_responsive.bed6
bedBinDistribution.pl -strand -smooth move -span 15 -t count --input HepG2_m6a_WTAP_SetD2_responsive.bed12 -bed6 $mRNAAnnotation -o ./HepG2_m6a_WTAP_SetD2_responsive.bed12bin

cat HepG2_m6a_shCont_writerShare_responsive.bed12 HepG2_m6a_shCont_SetD2_responsive.bed12 | sort | uniq -d > HepG2_m6a_writerShare_SetD2_responsive.bed12
exomePeakToSummit.pl --name writerShare_SetD2_responsive_m6a -b HepG2_m6a_writerShare_SetD2_responsive.bed12 -o $bed6Path/HepG2_m6a_writerShare_SetD2_responsive.bed6
bedBinDistribution.pl -strand -smooth move -span 15 -t count --input HepG2_m6a_writerShare_SetD2_responsive.bed12 -bed6 $mRNAAnnotation -o ./HepG2_m6a_writerShare_SetD2_responsive.bed12bin

cat HepG2_m6a_shCont_nonWriter_responsive.bed12 HepG2_m6a_shCont_SetD2_responsive.bed12 | sort | uniq -d > HepG2_m6a_nonWriter_SetD2_responsive.bed12
exomePeakToSummit.pl --name nonWriter_SetD2_responsive_m6a -b HepG2_m6a_nonWriter_SetD2_responsive.bed12 -o $bed6Path/HepG2_m6a_nonWriter_SetD2_responsive.bed6
bedBinDistribution.pl -strand -smooth move -span 15 -t count --input HepG2_m6a_nonWriter_SetD2_responsive.bed12 -bed6 $mRNAAnnotation -o ./HepG2_m6a_nonWriter_SetD2_responsive.bed12bin

paste HepG2_m6a_M14_SetD2_responsive.bed12bin HepG2_m6a_M3_SetD2_responsive.bed12bin HepG2_m6a_WTAP_SetD2_responsive.bed12bin \
  HepG2_m6a_writerShare_SetD2_responsive.bed12bin HepG2_m6a_nonWriter_SetD2_responsive.bed12bin |\
  cut -f 1,2,3,6,9,12,15 > HepG2_m6a_writer_SetD2_responsive.bed12bin
sed  -i '1i region\tbin\tM14_SetD2\tM3_SetD2\tWTAP_SetD2\twriterShare_SetD2\tnonWriter_SetD2' HepG2_m6a_writer_SetD2_responsive.bed12bin
rm -f HepG2_m6a_M14_SetD2_responsive.bed12bin HepG2_m6a_M3_SetD2_responsive.bed12bin HepG2_m6a_WTAP_SetD2_responsive.bed12bin \
  HepG2_m6a_writerShare_SetD2_responsive.bed12bin HepG2_m6a_nonWriter_SetD2_responsive.bed12bin

### writer-responsive + setd2-responsive + H3K36me3
bedtools intersect -a HepG2_m6a_M14_SetD2_responsive.bed12 -b $histonePeak -wa | sort -t $'\t' -k 1,1 -k 2,2n | uniq > HepG2_m6a_M14_SetD2_responsive_H3K36me3.bed12
exomePeakToSummit.pl --name writerShare_SetD2_responsive_m6a -b HepG2_m6a_M14_SetD2_responsive_H3K36me3.bed12 -o $bed6Path/HepG2_m6a_M14_SetD2_responsive_H3K36me3.bed6
bedBinDistribution.pl -strand -smooth move -span 15 -t count --input HepG2_m6a_M14_SetD2_responsive_H3K36me3.bed12 -bed6 $mRNAAnnotation -o ./HepG2_m6a_M14_SetD2_responsive_H3K36me3.bed12bin

bedtools intersect -a HepG2_m6a_M3_SetD2_responsive.bed12 -b $histonePeak -wa | sort -t $'\t' -k 1,1 -k 2,2n | uniq > HepG2_m6a_M3_SetD2_responsive_H3K36me3.bed12
exomePeakToSummit.pl --name writerShare_SetD2_responsive_m6a -b HepG2_m6a_M3_SetD2_responsive_H3K36me3.bed12 -o $bed6Path/HepG2_m6a_M3_SetD2_responsive_H3K36me3.bed6
bedBinDistribution.pl -strand -smooth move -span 15 -t count --input HepG2_m6a_M3_SetD2_responsive_H3K36me3.bed12 -bed6 $mRNAAnnotation -o ./HepG2_m6a_M3_SetD2_responsive_H3K36me3.bed12bin

bedtools intersect -a HepG2_m6a_WTAP_SetD2_responsive.bed12 -b $histonePeak -wa | sort -t $'\t' -k 1,1 -k 2,2n | uniq > HepG2_m6a_WTAP_SetD2_responsive_H3K36me3.bed12
exomePeakToSummit.pl --name writerShare_SetD2_responsive_m6a -b HepG2_m6a_WTAP_SetD2_responsive_H3K36me3.bed12 -o $bed6Path/HepG2_m6a_WTAP_SetD2_responsive_H3K36me3.bed6
bedBinDistribution.pl -strand -smooth move -span 15 -t count --input HepG2_m6a_WTAP_SetD2_responsive_H3K36me3.bed12 -bed6 $mRNAAnnotation -o ./HepG2_m6a_WTAP_SetD2_responsive_H3K36me3.bed12bin

bedtools intersect -a HepG2_m6a_writerShare_SetD2_responsive.bed12 -b $histonePeak -wa | sort -t $'\t' -k 1,1 -k 2,2n | uniq > HepG2_m6a_writerShare_SetD2_responsive_H3K36me3.bed12
exomePeakToSummit.pl --name writerShare_SetD2_responsive_m6a -b HepG2_m6a_writerShare_SetD2_responsive_H3K36me3.bed12 -o $bed6Path/HepG2_m6a_writerShare_SetD2_responsive_H3K36me3.bed6
bedBinDistribution.pl -strand -smooth move -span 15 -t count --input HepG2_m6a_writerShare_SetD2_responsive_H3K36me3.bed12 -bed6 $mRNAAnnotation -o ./HepG2_m6a_writerShare_SetD2_responsive_H3K36me3.bed12bin

bedtools intersect -a HepG2_m6a_nonWriter_SetD2_responsive.bed12 -b $histonePeak -wa | sort -t $'\t' -k 1,1 -k 2,2n | uniq > HepG2_m6a_nonWriter_SetD2_responsive_H3K36me3.bed12
exomePeakToSummit.pl --name writerShare_SetD2_responsive_m6a -b HepG2_m6a_nonWriter_SetD2_responsive_H3K36me3.bed12 -o $bed6Path/HepG2_m6a_nonWriter_SetD2_responsive_H3K36me3.bed6
bedBinDistribution.pl -strand -smooth move -span 15 -t count --input HepG2_m6a_nonWriter_SetD2_responsive_H3K36me3.bed12 -bed6 $mRNAAnnotation -o ./HepG2_m6a_nonWriter_SetD2_responsive_H3K36me3.bed12bin

paste HepG2_m6a_M14_SetD2_responsive_H3K36me3.bed12bin HepG2_m6a_M3_SetD2_responsive_H3K36me3.bed12bin HepG2_m6a_WTAP_SetD2_responsive_H3K36me3.bed12bin \
  HepG2_m6a_writerShare_SetD2_responsive_H3K36me3.bed12bin HepG2_m6a_nonWriter_SetD2_responsive_H3K36me3.bed12bin |\
  cut -f 1,2,3,6,9,12,15 > HepG2_m6a_writer_SetD2_responsive_H3K36me3.bed12bin
sed  -i '1i region\tbin\tM14_SetD2\tM3_SetD2\tWTAP_SetD2\twriterShare_SetD2\tnonWriter_SetD2' HepG2_m6a_writer_SetD2_responsive_H3K36me3.bed12bin
rm -f HepG2_m6a_M14_SetD2_responsive_H3K36me3.bed12bin HepG2_m6a_M3_SetD2_responsive_H3K36me3.bed12bin HepG2_m6a_WTAP_SetD2_responsive_H3K36me3.bed12bin \
  HepG2_m6a_writerShare_SetD2_responsive_H3K36me3.bed12bin HepG2_m6a_nonWriter_SetD2_responsive_H3K36me3.bed12bin


###draw all bed12 peak distribution
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_shCont_m6a.bed12 -bed6 $mRNAAnnotation -o ./HepG2_shCont_m6a.bed12bin
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_shSetD2_m6a.bed12 -bed6 $mRNAAnnotation -o ./HepG2_shSetD2_m6a.bed12bin
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_shM14_m6a.bed12 -bed6 $mRNAAnnotation -o ./HepG2_shM14_m6a.bed12bin
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_shM3_m6a.bed12 -bed6 $mRNAAnnotation -o ./HepG2_shM3_m6a.bed12bin
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_shWTAP_m6a.bed12 -bed6 $mRNAAnnotation -o ./HepG2_shWTAP_m6a.bed12bin
paste HepG2_shCont_m6a.bed12bin HepG2_shSetD2_m6a.bed12bin HepG2_shM14_m6a.bed12bin HepG2_shM3_m6a.bed12bin HepG2_shWTAP_m6a.bed12bin | cut -f 1,2,3,6,9,12,15 > HepG2_m6a.bed12bin
sed  -i '1i region\tbin\tshCont\tshSetD2\tshM14\tshM3\tshWTAP' HepG2_m6a.bed12bin
rm -f HepG2_shCont_m6a.bed12bin HepG2_shSetD2_m6a.bed12bin HepG2_shM14_m6a.bed12bin HepG2_shM3_m6a.bed12bin HepG2_shWTAP_m6a.bed12bin

###draw all bed12 peak + H3K36me3 distribution
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_shCont_m6a_H3K36me3.bed12 -bed6 $mRNAAnnotation -o ./HepG2_shCont_m6a_H3K36me3.bed12bin
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_shSetD2_m6a_H3K36me3.bed12 -bed6 $mRNAAnnotation -o ./HepG2_shSetD2_m6a_H3K36me3.bed12bin
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_shM14_m6a_H3K36me3.bed12 -bed6 $mRNAAnnotation -o ./HepG2_shM14_m6a_H3K36me3.bed12bin
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_shM3_m6a_H3K36me3.bed12 -bed6 $mRNAAnnotation -o ./HepG2_shM3_m6a_H3K36me3.bed12bin
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_shWTAP_m6a_H3K36me3.bed12 -bed6 $mRNAAnnotation -o ./HepG2_shWTAP_m6a_H3K36me3.bed12bin
paste HepG2_shCont_m6a_H3K36me3.bed12bin HepG2_shSetD2_m6a_H3K36me3.bed12bin HepG2_shM14_m6a_H3K36me3.bed12bin HepG2_shM3_m6a_H3K36me3.bed12bin HepG2_shWTAP_m6a_H3K36me3.bed12bin | cut -f 1,2,3,6,9,12,15 > HepG2_m6a_H3K36me3.bed12bin
sed  -i '1i region\tbin\tshCont\tshSetD2\tshM14\tshM3\tshWTAP' HepG2_m6a_H3K36me3.bed12bin
rm -f HepG2_shCont_m6a_H3K36me3.bed12bin HepG2_shSetD2_m6a_H3K36me3.bed12bin HepG2_shM14_m6a_H3K36me3.bed12bin HepG2_shM3_m6a_H3K36me3.bed12bin HepG2_shWTAP_m6a_H3K36me3.bed12bin

#### gene type counts and region counts of SetD2 responsive m6a sites
geneDistribution.pl -strand --input ./HepG2_m6a_shCont_SetD2_responsive.bed12 -bed6 $allAnnotation -o ./HepG2_m6a_shCont_SetD2_responsive.gene
sed  -i '1i geneType\tpeakNumber' ./HepG2_m6a_shCont_SetD2_responsive.gene
regionDistribution.pl -strand -size 200 -f '5utr,cds,stopCodon,3utr' --input ./HepG2_m6a_shCont_SetD2_responsive.bed12 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shCont_SetD2_responsive.region
sed  -i '1i region\tpeakNumber\tenrichment' ./HepG2_m6a_shCont_SetD2_responsive.region

regionDistribution.pl -strand -size 200 -f '5utr,cds,stopCodon,3utr' --input ./HepG2_m6a_M14_SetD2_responsive.bed12 -bed6 $mRNAAnnotation -o ./HepG2_m6a_M14_SetD2_responsive.region
regionDistribution.pl -strand -size 200 -f '5utr,cds,stopCodon,3utr' --input ./HepG2_m6a_M3_SetD2_responsive.bed12 -bed6 $mRNAAnnotation -o ./HepG2_m6a_M3_SetD2_responsive.region
regionDistribution.pl -strand -size 200 -f '5utr,cds,stopCodon,3utr' --input ./HepG2_m6a_WTAP_SetD2_responsive.bed12 -bed6 $mRNAAnnotation -o ./HepG2_m6a_WTAP_SetD2_responsive.region
regionDistribution.pl -strand -size 200 -f '5utr,cds,stopCodon,3utr' --input ./HepG2_m6a_writerShare_SetD2_responsive.bed12 -bed6 $mRNAAnnotation -o ./HepG2_m6a_writerShare_SetD2_responsive.region
regionDistribution.pl -strand -size 200 -f '5utr,cds,stopCodon,3utr' --input ./HepG2_m6a_nonWriter_SetD2_responsive.bed12 -bed6 $mRNAAnnotation -o ./HepG2_m6a_nonWriter_SetD2_responsive.region

paste HepG2_m6a_M14_SetD2_responsive.region HepG2_m6a_M3_SetD2_responsive.region HepG2_m6a_WTAP_SetD2_responsive.region \
  HepG2_m6a_writerShare_SetD2_responsive.region HepG2_m6a_nonWriter_SetD2_responsive.region |\
  cut -f 1,2,3,5,6,8,9,11,12,14,15 > HepG2_m6a_writers_SetD2_responsive_enrichment.region
sed  -i '1i region\tM14_SetD2_peak\tM14_SetD2_enrichment\tM3_SetD2_peak\tM3_SetD2_enrichment\tWTAP_SetD2_peak\tWTAP_SetD2_enrichment\tshare_peak\tshare_enrichment\tnonWriter_peak\tnonWriter_enrichment' ./HepG2_m6a_writers_SetD2_responsive_enrichment.region
rm -f HepG2_m6a_M14_SetD2_responsive.region HepG2_m6a_M3_SetD2_responsive.region HepG2_m6a_WTAP_SetD2_responsive.region \
  HepG2_m6a_writerShare_SetD2_responsive.region HepG2_m6a_nonWriter_SetD2_responsive.region

###draw bin distribution based on count and bed6 for writer targets
cd $bed6Path
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shCont_M14_responsive.bed6 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shCont_M14_responsive.bed6bin
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shSetD2_M14_responsive.bed6 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shSetD2_M14_responsive.bed6bin
paste HepG2_m6a_shCont_M14_responsive.bed6bin HepG2_m6a_shSetD2_M14_responsive.bed6bin | cut -f 1,2,3,6 > HepG2_m6a_M14_responsive.bed6bin
sed  -i '1i region\tbin\tshCont\tshSetD2' HepG2_m6a_M14_responsive.bed6bin
rm -f HepG2_m6a_shCont_M14_responsive.bed6bin HepG2_m6a_shSetD2_M14_responsive.bed6bin

bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shCont_M3_responsive.bed6 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shCont_M3_responsive.bed6bin
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shSetD2_M3_responsive.bed6 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shSetD2_M3_responsive.bed6bin
paste HepG2_m6a_shCont_M3_responsive.bed6bin HepG2_m6a_shSetD2_M3_responsive.bed6bin | cut -f 1,2,3,6 > HepG2_m6a_M3_responsive.bed6bin
sed  -i '1i region\tbin\tshCont\tshSetD2' HepG2_m6a_M3_responsive.bed6bin
rm -f HepG2_m6a_shCont_M3_responsive.bed6bin HepG2_m6a_shSetD2_M3_responsive.bed6bin

bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shCont_WTAP_responsive.bed6 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shCont_WTAP_responsive.bed6bin
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shSetD2_WTAP_responsive.bed6 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shSetD2_WTAP_responsive.bed6bin
paste HepG2_m6a_shCont_WTAP_responsive.bed6bin HepG2_m6a_shSetD2_WTAP_responsive.bed6bin | cut -f 1,2,3,6 > HepG2_m6a_WTAP_responsive.bed6bin
sed  -i '1i region\tbin\tshCont\tshSetD2' HepG2_m6a_WTAP_responsive.bed6bin
rm -f HepG2_m6a_shCont_WTAP_responsive.bed6bin HepG2_m6a_shSetD2_WTAP_responsive.bed6bin

bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shCont_writerShare_responsive.bed6 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shCont_writerShare_responsive.bed6bin
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shSetD2_writerShare_responsive.bed6 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shSetD2_writerShare_responsive.bed6bin
paste HepG2_m6a_shCont_writerShare_responsive.bed6bin HepG2_m6a_shSetD2_writerShare_responsive.bed6bin | cut -f 1,2,3,6 > HepG2_m6a_writerShare_responsive.bed6bin
sed  -i '1i region\tbin\tshCont\tshSetD2' HepG2_m6a_writerShare_responsive.bed6bin
rm -f HepG2_m6a_shCont_writerShare_responsive.bed6bin HepG2_m6a_shSetD2_writerShare_responsive.bed6bin

bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shCont_nonWriter_responsive.bed6 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shCont_nonWriter_responsive.bed6bin
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shSetD2_nonWriter_responsive.bed6 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shSetD2_nonWriter_responsive.bed6bin
paste HepG2_m6a_shCont_nonWriter_responsive.bed6bin HepG2_m6a_shSetD2_nonWriter_responsive.bed6bin | cut -f 1,2,3,6 > HepG2_m6a_nonWriter_responsive.bed6bin
sed  -i '1i region\tbin\tshCont\tshSetD2' HepG2_m6a_nonWriter_responsive.bed6bin
rm -f HepG2_m6a_shCont_nonWriter_responsive.bed6bin HepG2_m6a_shSetD2_nonWriter_responsive.bed6bin

###draw bin distribution based on count and bed6, overlapped with H3K36me3 peaks
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shCont_M14_responsive_H3K36me3.bed6 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shCont_M14_responsive_H3K36me3.bed6bin
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shSetD2_M14_responsive_H3K36me3.bed6 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shSetD2_M14_responsive_H3K36me3.bed6bin
paste HepG2_m6a_shCont_M14_responsive_H3K36me3.bed6bin HepG2_m6a_shSetD2_M14_responsive_H3K36me3.bed6bin | cut -f 1,2,3,6 > HepG2_m6a_M14_responsive_H3K36me3.bed6bin
sed  -i '1i region\tbin\tshCont\tshSetD2' HepG2_m6a_M14_responsive_H3K36me3.bed6bin
rm -f HepG2_m6a_shCont_M14_responsive_H3K36me3.bed6bin HepG2_m6a_shSetD2_M14_responsive_H3K36me3.bed6bin

bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shCont_M3_responsive_H3K36me3.bed6 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shCont_M3_responsive_H3K36me3.bed6bin
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shSetD2_M3_responsive_H3K36me3.bed6 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shSetD2_M3_responsive_H3K36me3.bed6bin
paste HepG2_m6a_shCont_M3_responsive_H3K36me3.bed6bin HepG2_m6a_shSetD2_M3_responsive_H3K36me3.bed6bin | cut -f 1,2,3,6 > HepG2_m6a_M3_responsive_H3K36me3.bed6bin
sed  -i '1i region\tbin\tshCont\tshSetD2' HepG2_m6a_M3_responsive_H3K36me3.bed6bin
rm -f HepG2_m6a_shCont_M3_responsive_H3K36me3.bed6bin HepG2_m6a_shSetD2_M3_responsive_H3K36me3.bed6bin

bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shCont_WTAP_responsive_H3K36me3.bed6 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shCont_WTAP_responsive_H3K36me3.bed6bin
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shSetD2_WTAP_responsive_H3K36me3.bed6 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shSetD2_WTAP_responsive_H3K36me3.bed6bin
paste HepG2_m6a_shCont_WTAP_responsive_H3K36me3.bed6bin HepG2_m6a_shSetD2_WTAP_responsive_H3K36me3.bed6bin | cut -f 1,2,3,6 > HepG2_m6a_WTAP_responsive_H3K36me3.bed6bin
sed  -i '1i region\tbin\tshCont\tshSetD2' HepG2_m6a_WTAP_responsive_H3K36me3.bed6bin
rm -f HepG2_m6a_shCont_WTAP_responsive_H3K36me3.bed6bin HepG2_m6a_shSetD2_WTAP_responsive_H3K36me3.bed6bin

bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shCont_writerShare_responsive_H3K36me3.bed6 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shCont_writerShare_responsive_H3K36me3.bed6bin
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shSetD2_writerShare_responsive_H3K36me3.bed6 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shSetD2_writerShare_responsive_H3K36me3.bed6bin
paste HepG2_m6a_shCont_writerShare_responsive_H3K36me3.bed6bin HepG2_m6a_shSetD2_writerShare_responsive_H3K36me3.bed6bin | cut -f 1,2,3,6 > HepG2_m6a_writerShare_responsive_H3K36me3.bed6bin
sed  -i '1i region\tbin\tshCont\tshSetD2' HepG2_m6a_writerShare_responsive_H3K36me3.bed6bin
rm -f HepG2_m6a_shCont_writerShare_responsive_H3K36me3.bed6bin HepG2_m6a_shSetD2_writerShare_responsive_H3K36me3.bed6bin

bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shCont_writerAll_responsive_H3K36me3.bed6 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shCont_writerAll_responsive_H3K36me3.bed6bin
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shSetD2_writerAll_responsive_H3K36me3.bed6 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shSetD2_writerAll_responsive_H3K36me3.bed6bin
paste HepG2_m6a_shCont_writerAll_responsive_H3K36me3.bed6bin HepG2_m6a_shSetD2_writerAll_responsive_H3K36me3.bed6bin | cut -f 1,2,3,6 > HepG2_m6a_writerAll_responsive_H3K36me3.bed6bin
sed  -i '1i region\tbin\tshCont\tshSetD2' HepG2_m6a_writerAll_responsive_H3K36me3.bed6bin
rm -f HepG2_m6a_shCont_writerAll_responsive_H3K36me3.bed6bin HepG2_m6a_shSetD2_writerAll_responsive_H3K36me3.bed6bin

bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shCont_nonWriter_responsive_H3K36me3.bed6 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shCont_nonWriter_responsive_H3K36me3.bed6bin
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shSetD2_nonWriter_responsive_H3K36me3.bed6 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shSetD2_nonWriter_responsive_H3K36me3.bed6bin
paste HepG2_m6a_shCont_nonWriter_responsive_H3K36me3.bed6bin HepG2_m6a_shSetD2_nonWriter_responsive_H3K36me3.bed6bin | cut -f 1,2,3,6 > HepG2_m6a_nonWriter_responsive_H3K36me3.bed6bin
sed  -i '1i region\tbin\tshCont\tshSetD2' HepG2_m6a_nonWriter_responsive_H3K36me3.bed6bin
rm -f HepG2_m6a_shCont_nonWriter_responsive_H3K36me3.bed6bin HepG2_m6a_shSetD2_nonWriter_responsive_H3K36me3.bed6bin

### writer-responsive + setd2-responsive
bedBinDistribution.pl -strand -smooth move -span 15 -t count --input HepG2_m6a_M14_SetD2_responsive.bed6 -bed6 $mRNAAnnotation -o ./HepG2_m6a_M14_SetD2_responsive.bed6bin
bedBinDistribution.pl -strand -smooth move -span 15 -t count --input HepG2_m6a_M3_SetD2_responsive.bed6 -bed6 $mRNAAnnotation -o ./HepG2_m6a_M3_SetD2_responsive.bed6bin
bedBinDistribution.pl -strand -smooth move -span 15 -t count --input HepG2_m6a_WTAP_SetD2_responsive.bed6 -bed6 $mRNAAnnotation -o ./HepG2_m6a_WTAP_SetD2_responsive.bed6bin
bedBinDistribution.pl -strand -smooth move -span 15 -t count --input HepG2_m6a_writerShare_SetD2_responsive.bed6 -bed6 $mRNAAnnotation -o ./HepG2_m6a_writerShare_SetD2_responsive.bed6bin
bedBinDistribution.pl -strand -smooth move -span 15 -t count --input HepG2_m6a_nonWriter_SetD2_responsive.bed6 -bed6 $mRNAAnnotation -o ./HepG2_m6a_nonWriter_SetD2_responsive.bed6bin

paste HepG2_m6a_M14_SetD2_responsive.bed6bin HepG2_m6a_M3_SetD2_responsive.bed6bin HepG2_m6a_WTAP_SetD2_responsive.bed6bin \
  HepG2_m6a_writerShare_SetD2_responsive.bed6bin HepG2_m6a_nonWriter_SetD2_responsive.bed6bin |\
  cut -f 1,2,3,6,9,12,15 > HepG2_m6a_writer_SetD2_responsive.bed6bin
sed  -i '1i region\tbin\tM14_SetD2\tM3_SetD2\tWTAP_SetD2\twriterShare_SetD2\tnonWriter_SetD2' HepG2_m6a_writer_SetD2_responsive.bed6bin
rm -f HepG2_m6a_M14_SetD2_responsive.bed6bin HepG2_m6a_M3_SetD2_responsive.bed6bin HepG2_m6a_WTAP_SetD2_responsive.bed6bin \
  HepG2_m6a_writerShare_SetD2_responsive.bed6bin HepG2_m6a_nonWriter_SetD2_responsive.bed6bin

### writer-responsive + setd2-responsive + H3K36me3
bedBinDistribution.pl -strand -smooth move -span 15 -t count --input HepG2_m6a_M14_SetD2_responsive_H3K36me3.bed6 -bed6 $mRNAAnnotation -o ./HepG2_m6a_M14_SetD2_responsive_H3K36me3.bed6bin
bedBinDistribution.pl -strand -smooth move -span 15 -t count --input HepG2_m6a_M3_SetD2_responsive_H3K36me3.bed6 -bed6 $mRNAAnnotation -o ./HepG2_m6a_M3_SetD2_responsive_H3K36me3.bed6bin
bedBinDistribution.pl -strand -smooth move -span 15 -t count --input HepG2_m6a_WTAP_SetD2_responsive_H3K36me3.bed6 -bed6 $mRNAAnnotation -o ./HepG2_m6a_WTAP_SetD2_responsive_H3K36me3.bed6bin
bedBinDistribution.pl -strand -smooth move -span 15 -t count --input HepG2_m6a_writerShare_SetD2_responsive_H3K36me3.bed6 -bed6 $mRNAAnnotation -o ./HepG2_m6a_writerShare_SetD2_responsive_H3K36me3.bed6bin
bedBinDistribution.pl -strand -smooth move -span 15 -t count --input HepG2_m6a_nonWriter_SetD2_responsive_H3K36me3.bed6 -bed6 $mRNAAnnotation -o ./HepG2_m6a_nonWriter_SetD2_responsive_H3K36me3.bed6bin

paste HepG2_m6a_M14_SetD2_responsive_H3K36me3.bed6bin HepG2_m6a_M3_SetD2_responsive_H3K36me3.bed6bin HepG2_m6a_WTAP_SetD2_responsive_H3K36me3.bed6bin \
  HepG2_m6a_writerShare_SetD2_responsive_H3K36me3.bed6bin HepG2_m6a_nonWriter_SetD2_responsive_H3K36me3.bed6bin |\
  cut -f 1,2,3,6,9,12,15 > HepG2_m6a_writer_SetD2_responsive_H3K36me3.bed6bin
sed  -i '1i region\tbin\tM14_SetD2\tM3_SetD2\tWTAP_SetD2\twriterShare_SetD2\tnonWriter_SetD2' HepG2_m6a_writer_SetD2_responsive_H3K36me3.bed6bin
rm -f HepG2_m6a_M14_SetD2_responsive_H3K36me3.bed6bin HepG2_m6a_M3_SetD2_responsive_H3K36me3.bed6bin HepG2_m6a_WTAP_SetD2_responsive_H3K36me3.bed6bin \
  HepG2_m6a_writerShare_SetD2_responsive_H3K36me3.bed6bin HepG2_m6a_nonWriter_SetD2_responsive_H3K36me3.bed6bin


###draw all bed6 peak distribution
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_shCont_m6a.bed6 -bed6 $mRNAAnnotation -o ./HepG2_shCont_m6a.bed6bin
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_shSetD2_m6a.bed6 -bed6 $mRNAAnnotation -o ./HepG2_shSetD2_m6a.bed6bin
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_shM14_m6a.bed6 -bed6 $mRNAAnnotation -o ./HepG2_shM14_m6a.bed6bin
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_shM3_m6a.bed6 -bed6 $mRNAAnnotation -o ./HepG2_shM3_m6a.bed6bin
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_shWTAP_m6a.bed6 -bed6 $mRNAAnnotation -o ./HepG2_shWTAP_m6a.bed6bin
paste HepG2_shCont_m6a.bed6bin HepG2_shSetD2_m6a.bed6bin HepG2_shM14_m6a.bed6bin HepG2_shM3_m6a.bed6bin HepG2_shWTAP_m6a.bed6bin | cut -f 1,2,3,6,9,12,15 > HepG2_m6a.bed6bin
sed  -i '1i region\tbin\tshCont\tshSetD2\tshM14\tshM3\tshWTAP' HepG2_m6a.bed6bin
rm -f HepG2_shCont_m6a.bed6bin HepG2_shSetD2_m6a.bed6bin HepG2_shM14_m6a.bed6bin HepG2_shM3_m6a.bed6bin HepG2_shWTAP_m6a.bed6bin

bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_shCont_m6a_H3K36me3.bed6 -bed6 $mRNAAnnotation -o ./HepG2_shCont_m6a_H3K36me3.bed6bin
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_shSetD2_m6a_H3K36me3.bed6 -bed6 $mRNAAnnotation -o ./HepG2_shSetD2_m6a_H3K36me3.bed6bin
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_shM14_m6a_H3K36me3.bed6 -bed6 $mRNAAnnotation -o ./HepG2_shM14_m6a_H3K36me3.bed6bin
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_shM3_m6a_H3K36me3.bed6 -bed6 $mRNAAnnotation -o ./HepG2_shM3_m6a_H3K36me3.bed6bin
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_shWTAP_m6a_H3K36me3.bed6 -bed6 $mRNAAnnotation -o ./HepG2_shWTAP_m6a_H3K36me3.bed6bin
paste HepG2_shCont_m6a_H3K36me3.bed6bin HepG2_shSetD2_m6a_H3K36me3.bed6bin HepG2_shM14_m6a_H3K36me3.bed6bin HepG2_shM3_m6a_H3K36me3.bed6bin HepG2_shWTAP_m6a_H3K36me3.bed6bin | cut -f 1,2,3,6,9,12,15 > HepG2_m6a_H3K36me3.bed6bin
sed  -i '1i region\tbin\tshCont\tshSetD2\tshM14\tshM3\tshWTAP' HepG2_m6a_H3K36me3.bed6bin
rm -f HepG2_shCont_m6a_H3K36me3.bed6bin HepG2_shSetD2_m6a_H3K36me3.bed6bin HepG2_shM14_m6a_H3K36me3.bed6bin HepG2_shM3_m6a_H3K36me3.bed6bin HepG2_shWTAP_m6a_H3K36me3.bed6bin

#### gene type counts and region counts of SetD2 responsive m6a sites
geneDistribution.pl -strand --input ./HepG2_m6a_shCont_SetD2_responsive.bed6 -bed6 $allAnnotation -o ./HepG2_m6a_shCont_SetD2_responsive.gene
sed  -i '1i geneType\tpeakNumber' ./HepG2_m6a_shCont_SetD2_responsive.gene
regionDistribution.pl -strand -size 200 -f '5utr,cds,stopCodon,3utr' --input ./HepG2_m6a_shCont_SetD2_responsive.bed6 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shCont_SetD2_responsive.region
sed  -i '1i region\tpeakNumber\tenrichment' ./HepG2_m6a_shCont_SetD2_responsive.region
