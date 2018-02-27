#!/bin/sh

cd /data/zhoukr/hhl_setd2_m6a/HepG2_ChIP-seq

## xls
rm -rf /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/venn
mkdir /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/venn
cd /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/venn
### peak fold_enrichment
bedtools intersect -a ../bed6/HepG2_macs_shCont_peaks.bed -b ../bed6/HepG2_macs_shSetD2_peaks.bed | sort -t $'\t' -k 1,1V -k 2,2n > HepG2_macs_shCont_shSetD2_interset.bed
bedtools merge -i HepG2_macs_shCont_shSetD2_interset.bed | awk 'BEGIN{OFS="\t";FS="\t";count=1};{name="mergePeak="FNR;print $0,name,0,"+";}' > HepG2_macs_shCont_shSetD2_interset_merge.bed

bedtools intersect -a ../bed6/HepG2_macs_shCont_peaks.bed -b HepG2_macs_shCont_shSetD2_interset_merge.bed -wa -wb > HepG2_macs_shCont_interset_merge.txt
bedtools intersect -a ../bed6/HepG2_macs_shSetD2_peaks.bed -b HepG2_macs_shCont_shSetD2_interset_merge.bed -wa -wb > HepG2_macs_shSetD2_interset_merge.txt
awk '
BEGIN{OFS="\t";FS="\t"};
ARGIND==1{hashPeakName[$4]=$10;}
ARGIND==2{if($4 in hashPeakName) {print hashPeakName[$4];}else{print "peakShCont="FNR;}}
' \
HepG2_macs_shCont_interset_merge.txt ../bed6/HepG2_macs_shCont_peaks.bed | sort | uniq > HepG2_macs_shCont_peaks_name.txt

awk '
BEGIN{OFS="\t";FS="\t"};
ARGIND==1{hashPeakName[$4]=$10;}
ARGIND==2{if($4 in hashPeakName) {print hashPeakName[$4];}else{print "peakShSetD2="FNR;}}
' \
HepG2_macs_shSetD2_interset_merge.txt ../bed6/HepG2_macs_shSetD2_peaks.bed | sort | uniq > HepG2_macs_shSetD2_peaks_name.txt

rm -f HepG2_macs_shCont_shSetD2_interset.bed HepG2_macs_shCont_shSetD2_interset_merge.bed HepG2_macs_shCont_interset_merge.txt HepG2_macs_shSetD2_interset_merge.txt
