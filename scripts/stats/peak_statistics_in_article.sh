#!/bin/sh

targetPath=/data/zhoukr/hhl_setd2_m6a/analysis/stats
heHelaPath=/data/zhoukr/hhl_setd2_m6a/others/hechuan/Hela/tophat

echo -e "###ChIP-seq peak statistics(PeakNum)"
echo -e "cell-line\tshCont\tshSetD2\tshM14\tshM3\tshWTAP"

### HepG2, ChIP-seq
cd /data/zhoukr/hhl_setd2_m6a/HepG2_ChIP-seq
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x shCont/macs1.4-merge/shCont_peaks.xls -o $targetPath/shCont_peaks.bed
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x shSetD2/macs1.4-merge/shSetD2_peaks.xls -o $targetPath/shSetD2_peaks.bed
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x shM14/macs1.4/shM14_rep1_peaks.xls -o $targetPath/shM14_peaks.bed
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x shM3/macs1.4/shM3_rep1_peaks.xls -o $targetPath/shM3_peaks.bed
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x shWTAP/macs1.4/shWTAP_rep1_peaks.xls -o $targetPath/shWTAP_peaks.bed

shContPeakNum=`cat $targetPath/shCont_peaks.bed | wc -l`
shSetD2PeakNum=`cat $targetPath/shSetD2_peaks.bed | wc -l`
shM14PeakNum=`cat $targetPath/shM14_peaks.bed | wc -l`
shM3PeakNum=`cat $targetPath/shM3_peaks.bed | wc -l`
shWTAPPeakNum=`cat $targetPath/shWTAP_peaks.bed | wc -l`

echo -e "HepG2\t${shContPeakNum}\t${shSetD2PeakNum}\t${shM14PeakNum}\t${shM3PeakNum}\t${shWTAPPeakNum}"



### Hela, ChIP-seq
cd /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x shCont/macs1.4/rep2/shCont_rep2_peaks.xls -o $targetPath/shCont_peaks.bed
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x shSetD2/macs1.4/rep2/shSetD2_rep2_peaks.xls -o $targetPath/shSetD2_peaks.bed
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x shM14/macs1.4/shM14_rep1_peaks.xls -o $targetPath/shM14_peaks.bed
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x shM3/macs1.4/shM3_rep1_peaks.xls -o $targetPath/shM3_peaks.bed
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x shWTAP/macs1.4/shWTAP_rep1_peaks.xls -o $targetPath/shWTAP_peaks.bed

shContPeakNum=`cat $targetPath/shCont_peaks.bed | wc -l`
shSetD2PeakNum=`cat $targetPath/shSetD2_peaks.bed | wc -l`
shM14PeakNum=`cat $targetPath/shM14_peaks.bed | wc -l`
shM3PeakNum=`cat $targetPath/shM3_peaks.bed | wc -l`
shWTAPPeakNum=`cat $targetPath/shWTAP_peaks.bed | wc -l`

echo -e "Hela\t${shContPeakNum}\t${shSetD2PeakNum}\t${shM14PeakNum}\t${shM3PeakNum}\t${shWTAPPeakNum}"

### mESCs, ChIP-seq
cd /data/zhoukr/hhl_setd2_m6a/mouse_ESCs_ChIP-seq/
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x Ctrl_D0/macs1.4-merge/shCont_peaks.xls -name shCont_macs -o $targetPath/shCont_peaks.bed
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x SetD2-KD_D0/macs1.4-merge/shSetD2_peaks.xls -name shSetD2_macs -o $targetPath/shSetD2_peaks.bed
shContPeakNum=`cat $targetPath/shCont_peaks.bed | wc -l`
shSetD2PeakNum=`cat $targetPath/shSetD2_peaks.bed | wc -l`
echo -e "mESCs-D0\t${shContPeakNum}\t${shSetD2PeakNum}"

### mMEFs, ChIP-seq
cd /data/zhoukr/hhl_setd2_m6a/mouse_MEFs_ChIP-seq/
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x SetD2-WT_NA/macs1.4-merge/shCont_peaks.xls -name shCont_macs -o $targetPath/shCont_peaks.bed
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x SetD2-KO_NA/macs1.4-merge/shSetD2_peaks.xls -name shSetD2_macs -o $targetPath/shSetD2_peaks.bed
shContPeakNum=`cat $targetPath/shCont_peaks.bed | wc -l`
shSetD2PeakNum=`cat $targetPath/shSetD2_peaks.bed | wc -l`
echo -e "mEFs\t${shContPeakNum}\t${shSetD2PeakNum}"



echo -e "\n\n###m6A-seq peak statistics(PeakNum)"
echo -e "cell-line\tshCont\tshSetD2\tshM14\tshM3\tshWTAP"
### HepG2, m6A-seq
cd /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq

exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shCont/shCont_m6A_rep1-2-3/peak.xls -o $targetPath/shCont_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shSetD2/shSetD2_m6A_rep1-2-3/peak.xls -o $targetPath/shSetD2_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shM14/shM14_m6A_rep1-2-3/peak.xls -o $targetPath/shM14_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shM3/shM3_m6A_rep1-2-3/peak.xls -o $targetPath/shM3_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shWTAP/shWTAP_m6A_rep1-2-3/peak.xls -o $targetPath/shWTAP_peaks.bed

shContPeakNum=`cat $targetPath/shCont_peaks.bed | wc -l`
shSetD2PeakNum=`cat $targetPath/shSetD2_peaks.bed | wc -l`
shM14PeakNum=`cat $targetPath/shM14_peaks.bed | wc -l`
shM3PeakNum=`cat $targetPath/shM3_peaks.bed | wc -l`
shWTAPPeakNum=`cat $targetPath/shWTAP_peaks.bed | wc -l`
echo -e "HepG2\t${shContPeakNum}\t${shSetD2PeakNum}\t${shM14PeakNum}\t${shM3PeakNum}\t${shWTAPPeakNum}"

### Hela, m6A-seq
cd /data/zhoukr/hhl_setd2_m6a/Hela_m6A-seq
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shCont/shCont_m6A_rep3/peak.xls -o $targetPath/shCont_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shSetD2/shSetD2_m6A_rep3/peak.xls -o $targetPath/shSetD2_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shM14/shM14_m6A_rep3/peak.xls -o $targetPath/shM14_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shM3/shM3_m6A_rep3/peak.xls -o $targetPath/shM3_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shWTAP/shWTAP_m6A_rep3/peak.xls -o $targetPath/shWTAP_peaks.bed

shContPeakNum=`cat $targetPath/shCont_peaks.bed | wc -l`
shSetD2PeakNum=`cat $targetPath/shSetD2_peaks.bed | wc -l`
shM14PeakNum=`cat $targetPath/shM14_peaks.bed | wc -l`
shM3PeakNum=`cat $targetPath/shM3_peaks.bed | wc -l`
shWTAPPeakNum=`cat $targetPath/shWTAP_peaks.bed | wc -l`
echo -e "Hela\t${shContPeakNum}\t${shSetD2PeakNum}\t${shM14PeakNum}\t${shM3PeakNum}\t${shWTAPPeakNum}"

### mESCs, m6A-seq
cd /data/zhoukr/hhl_setd2_m6a/mouse_ESCs_m6A-seq
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -x Ctrl_D0/Ctrl_D0_m6A/peak.xls -o $targetPath/shCont_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -x SetD2-KD_D0/SetD2-KD_D0_m6A/peak.xls -o $targetPath/shSetD2_peaks.bed

shContPeakNum=`cat $targetPath/shCont_peaks.bed | wc -l`
shSetD2PeakNum=`cat $targetPath/shSetD2_peaks.bed | wc -l`
echo -e "mESCs-D0\t${shContPeakNum}\t${shSetD2PeakNum}\tNA\tNA\tNA"

exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -x Ctrl_D6/Ctrl_D6_m6A/peak.xls -o $targetPath/shCont_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -x SetD2-KD_D6/SetD2-KD_D6_m6A/peak.xls -o $targetPath/shSetD2_peaks.bed
shContPeakNum=`cat $targetPath/shCont_peaks.bed | wc -l`
shSetD2PeakNum=`cat $targetPath/shSetD2_peaks.bed | wc -l`
echo -e "mESCs-D6\t${shContPeakNum}\t${shSetD2PeakNum}\tNA\tNA\tNA"

### mMEFs, m6A-seq
cd /data/zhoukr/hhl_setd2_m6a/mouse_MEFs_m6A-seq
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -x SetD2-WT_NA/SetD2-WT_NA_m6A/peak.xls -o $targetPath/shCont_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -x SetD2-KO_NA/SetD2-KO_NA_m6A/peak.xls -o $targetPath/shSetD2_peaks.bed
shContPeakNum=`cat $targetPath/shCont_peaks.bed | wc -l`
shSetD2PeakNum=`cat $targetPath/shSetD2_peaks.bed | wc -l`
echo -e "mEFs\t${shContPeakNum}\t${shSetD2PeakNum}\tNA\tNA\tNA"


echo -e "\n\n###protein-responsive peaks statistics (m6A-seq)"
echo -e "cell-line\tSETD2-responsive\tMETTL14-responsive\tMETTL3-responsive\tWTAP-responsive\tnon-writer-responsive"
### HepG2, m6A-seq
cd /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/bed12
Setd2_responsive=`cat HepG2_m6a_shCont_SetD2_responsive.bed12 | wc -l`
M14_responsive=`cat HepG2_m6a_shCont_M14_responsive.bed12 | wc -l`
M3_responsive=`cat HepG2_m6a_shCont_M3_responsive.bed12 | wc -l`
WTAP_responsive=`cat HepG2_m6a_shCont_WTAP_responsive.bed12 | wc -l`
nonWriter_responsive=`cat HepG2_m6a_shCont_nonWriter_responsive.bed12 | wc -l`
echo -e "HepG2\t${Setd2_responsive}\t${M14_responsive}\t${M3_responsive}\t${WTAP_responsive}\t${nonWriter_responsive}"

cd /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/Hela/bed12
Setd2_responsive=`cat Hela_m6a_shCont_SetD2_responsive.bed12 | wc -l`
M14_responsive=`cat Hela_m6a_shCont_M14_responsive.bed12 | wc -l`
M3_responsive=`cat Hela_m6a_shCont_M3_responsive.bed12 | wc -l`
WTAP_responsive=`cat Hela_m6a_shCont_WTAP_responsive.bed12 | wc -l`
nonWriter_responsive=`cat Hela_m6a_shCont_nonWriter_responsive.bed12 | wc -l`
echo -e "Hela\t${Setd2_responsive}\t${M14_responsive}\t${M3_responsive}\t${WTAP_responsive}\t${nonWriter_responsive}"

cd /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/mESCs/bed12
D0Setd2_responsive=`cat mESCs_m6a_shCont_SetD2_responsive-D0.bed12 | wc -l`
D6Setd2_responsive=`cat mESCs_m6a_shCont_SetD2_responsive-D6.bed12 | wc -l`
echo -e "mESCs-D0\t${D0Setd2_responsive}\tNA\tNA\tNA\tNA"
echo -e "mESCs-D6\t${D6Setd2_responsive}\tNA\tNA\tNA\tNA"

cd /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/mMEFs/bed12
Setd2_responsive=`cat mMEFs_m6a_shCont_SetD2_responsive.bed12 | wc -l`
echo -e "mEFs\t${Setd2_responsive}\tNA\tNA\tNA\tNA"

histone=/data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/bed6/HepG2_macs_shCont_peaks.bed
histoneResonse=/data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/bed6/HepG2_shCont_SetD2_responsive.bed
m6A=/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/bed12/HepG2_shCont_m6a.bed12
m6aResponse=/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/bed12/HepG2_m6a_shCont_SetD2_responsive.bed12
scatter=/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/scatter/HepG2_m6a_peak_shSetD2_FoldEnrichment.txt
awk '{if(FNR>1){if(($14+0.1)/($13+0.1)>2) {for(i=1;i<12;i++) {printf("%s\t",$i);} printf("%s", $12); printf "\n";}}}' $scatter > ./HepG2m6AHyper.bed
overlapped=`bedtools intersect -a $histone -b $m6A -u | wc -l`
echo -e '##H3K36me3 overlapped with m6A'
echo -e "H3K36me3 overlapped with m6a(loci where H3K36me3 and m6A are both present)\t$overlapped"
overlapped=`bedtools intersect -a $histone -b $m6A -u | bedtools intersect -a stdin -b $histoneResonse -u | wc -l`
echo -e "hypomethylated H3K36me3 peaks in loci where H3K36me3 and m6A are both present\t$overlapped"
overlapped=`bedtools intersect -a $histone -b $m6A -u | bedtools intersect -a stdin -b $histoneResonse -u | bedtools intersect -a stdin -b $m6aResponse -u | wc -l`
echo -e "hypomethylated H3K36me3 peaks in loci where H3K36me3 and m6A are both present overlapped with hypomethylated m6A\t$overlapped"
hypoHyper=`bedtools intersect -a $histone -b $m6A -u | bedtools intersect -a stdin -b $histoneResonse -u | bedtools intersect -a stdin -b ./HepG2m6AHyper.bed -u | wc -l`
echo -e "hypomethylated H3K36me3 peaks in loci where H3K36me3 and m6A are both present overlapped with hypermethylated m6A\t$hypoHyper"

cd $targetPath/
rm -f ./*.bed
