#!/bin/sh

cd /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq

## xls
rm -rf /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/Hela/venn
mkdir /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/Hela/venn
cd /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/Hela/venn
### peak fold_enrichment
bedtools intersect -a ../bed6/Hela_macs_shCont_peaks.bed -b ../bed6/Hela_macs_shSetD2_peaks.bed | sort -t $'\t' -k 1,1V -k 2,2n > Hela_macs_shCont_shSetD2_interset.bed
bedtools merge -i Hela_macs_shCont_shSetD2_interset.bed | awk 'BEGIN{OFS="\t";FS="\t";count=1};{name="mergePeak="FNR;print $0,name,0,"+";}' > Hela_macs_shCont_shSetD2_interset_merge.bed

bedtools intersect -a ../bed6/Hela_macs_shCont_peaks.bed -b Hela_macs_shCont_shSetD2_interset_merge.bed -wa -wb > Hela_macs_shCont_interset_merge.txt
bedtools intersect -a ../bed6/Hela_macs_shSetD2_peaks.bed -b Hela_macs_shCont_shSetD2_interset_merge.bed -wa -wb > Hela_macs_shSetD2_interset_merge.txt
awk '
BEGIN{OFS="\t";FS="\t"};
ARGIND==1{hashPeakName[$4]=$10;}
ARGIND==2{if($4 in hashPeakName) {print hashPeakName[$4];}else{print "peakShCont="FNR;}}
' \
Hela_macs_shCont_interset_merge.txt ../bed6/Hela_macs_shCont_peaks.bed | sort | uniq > Hela_macs_shCont_peaks_name.txt

awk '
BEGIN{OFS="\t";FS="\t"};
ARGIND==1{hashPeakName[$4]=$10;}
ARGIND==2{if($4 in hashPeakName) {print hashPeakName[$4];}else{print "peakShSetD2="FNR;}}
' \
Hela_macs_shSetD2_interset_merge.txt ../bed6/Hela_macs_shSetD2_peaks.bed | sort | uniq > Hela_macs_shSetD2_peaks_name.txt

rm -f Hela_macs_shCont_shSetD2_interset.bed Hela_macs_shCont_shSetD2_interset_merge.bed Hela_macs_shCont_interset_merge.txt Hela_macs_shSetD2_interset_merge.txt
