#!/bin/sh

targetPath=/data/zhoukr/hhl_setd2_m6a/analysis/stats
heHelaPath=/data/zhoukr/hhl_setd2_m6a/others/hechuan/Hela/tophat

echo -e "###ChIP-seq peak statistics(PeakNum)"
echo -e "sample\tshCont\tshSetD2\tshM14\tshM3\tshWTAP"

### HepG2, ChIP-seq
cd /data/zhoukr/hhl_setd2_m6a/HepG2_ChIP-seq
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x shCont/macs1.4/rep1/shCont_rep1_peaks.xls -o $targetPath/shCont_peaks.bed
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x shSetD2/macs1.4/rep1/shSetD2_rep1_peaks.xls -o $targetPath/shSetD2_peaks.bed
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x shM14/macs1.4/shM14_rep1_peaks.xls -o $targetPath/shM14_peaks.bed
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x shM3/macs1.4/shM3_rep1_peaks.xls -o $targetPath/shM3_peaks.bed
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x shWTAP/macs1.4/shWTAP_rep1_peaks.xls -o $targetPath/shWTAP_peaks.bed

shContPeakNum=`cat $targetPath/shCont_peaks.bed | wc -l`
shSetD2PeakNum=`cat $targetPath/shSetD2_peaks.bed | wc -l`
shM14PeakNum=`cat $targetPath/shM14_peaks.bed | wc -l`
shM3PeakNum=`cat $targetPath/shM3_peaks.bed | wc -l`
shWTAPPeakNum=`cat $targetPath/shWTAP_peaks.bed | wc -l`

echo -e "HepG2_rep1\t${shContPeakNum}\t${shSetD2PeakNum}\t${shM14PeakNum}\t${shM3PeakNum}\t${shWTAPPeakNum}"

macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x shCont/macs1.4/rep2/shCont_rep2_peaks.xls -o $targetPath/shCont_peaks.bed
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x shSetD2/macs1.4/rep2/shSetD2_rep2_peaks.xls -o $targetPath/shSetD2_peaks.bed

shContPeakNum=`cat $targetPath/shCont_peaks.bed | wc -l`
shSetD2PeakNum=`cat $targetPath/shSetD2_peaks.bed | wc -l`

echo -e "HepG2_rep2\t${shContPeakNum}\t${shSetD2PeakNum}\tNA\tNA\tNA"

macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x shCont/macs1.4-merge/shCont_peaks.xls -o $targetPath/shCont_peaks.bed
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x shSetD2/macs1.4-merge/shSetD2_peaks.xls -o $targetPath/shSetD2_peaks.bed

shContPeakNum=`cat $targetPath/shCont_peaks.bed | wc -l`
shSetD2PeakNum=`cat $targetPath/shSetD2_peaks.bed | wc -l`

echo -e "HepG2_rep1-2\t${shContPeakNum}\t${shSetD2PeakNum}\tNA\tNA\tNA"


### Hela, ChIP-seq
cd /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x shCont/macs1.4/rep1/shCont_rep1_peaks.xls -o $targetPath/shCont_peaks.bed
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x shSetD2/macs1.4/rep1/shSetD2_rep1_peaks.xls -o $targetPath/shSetD2_peaks.bed
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x shM14/macs1.4/shM14_rep1_peaks.xls -o $targetPath/shM14_peaks.bed
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x shM3/macs1.4/shM3_rep1_peaks.xls -o $targetPath/shM3_peaks.bed
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x shWTAP/macs1.4/shWTAP_rep1_peaks.xls -o $targetPath/shWTAP_peaks.bed

shContPeakNum=`cat $targetPath/shCont_peaks.bed | wc -l`
shSetD2PeakNum=`cat $targetPath/shSetD2_peaks.bed | wc -l`
shM14PeakNum=`cat $targetPath/shM14_peaks.bed | wc -l`
shM3PeakNum=`cat $targetPath/shM3_peaks.bed | wc -l`
shWTAPPeakNum=`cat $targetPath/shWTAP_peaks.bed | wc -l`

echo -e "Hela_rep1\t${shContPeakNum}\t${shSetD2PeakNum}\t${shM14PeakNum}\t${shM3PeakNum}\t${shWTAPPeakNum}"

macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x shCont/macs1.4/rep2/shCont_rep2_peaks.xls -o $targetPath/shCont_peaks.bed
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x shSetD2/macs1.4/rep2/shSetD2_rep2_peaks.xls -o $targetPath/shSetD2_peaks.bed

shContPeakNum=`cat $targetPath/shCont_peaks.bed | wc -l`
shSetD2PeakNum=`cat $targetPath/shSetD2_peaks.bed | wc -l`

echo -e "Hela_rep2\t${shContPeakNum}\t${shSetD2PeakNum}"

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
echo -e "mMEFs\t${shContPeakNum}\t${shSetD2PeakNum}"



echo -e "\n\n###m6A-seq peak statistics(PeakNum)"
echo -e "samples\tshCont\tshSetD2\tshM14\tshM3\tshWTAP"
### HepG2, m6A-seq
cd /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shCont/shCont_m6A_rep1/peak.xls -o $targetPath/shCont_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shSetD2/shSetD2_m6A_rep1/peak.xls -o $targetPath/shSetD2_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shM14/shM14_m6A_rep1/peak.xls -o $targetPath/shM14_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shM3/shM3_m6A_rep1/peak.xls -o $targetPath/shM3_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shWTAP/shWTAP_m6A_rep1/peak.xls -o $targetPath/shWTAP_peaks.bed

shContPeakNum=`cat $targetPath/shCont_peaks.bed | wc -l`
shSetD2PeakNum=`cat $targetPath/shSetD2_peaks.bed | wc -l`
shM14PeakNum=`cat $targetPath/shM14_peaks.bed | wc -l`
shM3PeakNum=`cat $targetPath/shM3_peaks.bed | wc -l`
shWTAPPeakNum=`cat $targetPath/shWTAP_peaks.bed | wc -l`
echo -e "HepG2_rep1\t${shContPeakNum}\t${shSetD2PeakNum}\t${shM14PeakNum}\t${shM3PeakNum}\t${shWTAPPeakNum}"

exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shCont/shCont_m6A_rep2/peak.xls -o $targetPath/shCont_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shSetD2/shSetD2_m6A_rep2/peak.xls -o $targetPath/shSetD2_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shM14/shM14_m6A_rep2/peak.xls -o $targetPath/shM14_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shM3/shM3_m6A_rep2/peak.xls -o $targetPath/shM3_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shWTAP/shWTAP_m6A_rep2/peak.xls -o $targetPath/shWTAP_peaks.bed

shContPeakNum=`cat $targetPath/shCont_peaks.bed | wc -l`
shSetD2PeakNum=`cat $targetPath/shSetD2_peaks.bed | wc -l`
shM14PeakNum=`cat $targetPath/shM14_peaks.bed | wc -l`
shM3PeakNum=`cat $targetPath/shM3_peaks.bed | wc -l`
shWTAPPeakNum=`cat $targetPath/shWTAP_peaks.bed | wc -l`
echo -e "HepG2_rep2\t${shContPeakNum}\t${shSetD2PeakNum}\t${shM14PeakNum}\t${shM3PeakNum}\t${shWTAPPeakNum}"

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
echo -e "HepG2_rep3\t${shContPeakNum}\t${shSetD2PeakNum}\t${shM14PeakNum}\t${shM3PeakNum}\t${shWTAPPeakNum}"

exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shCont/shCont_m6A_rep1-2/peak.xls -o $targetPath/shCont_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shSetD2/shSetD2_m6A_rep1-2/peak.xls -o $targetPath/shSetD2_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shM14/shM14_m6A_rep1-2/peak.xls -o $targetPath/shM14_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shM3/shM3_m6A_rep1-2/peak.xls -o $targetPath/shM3_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shWTAP/shWTAP_m6A_rep1-2/peak.xls -o $targetPath/shWTAP_peaks.bed

shContPeakNum=`cat $targetPath/shCont_peaks.bed | wc -l`
shSetD2PeakNum=`cat $targetPath/shSetD2_peaks.bed | wc -l`
shM14PeakNum=`cat $targetPath/shM14_peaks.bed | wc -l`
shM3PeakNum=`cat $targetPath/shM3_peaks.bed | wc -l`
shWTAPPeakNum=`cat $targetPath/shWTAP_peaks.bed | wc -l`
echo -e "HepG2_rep1-2\t${shContPeakNum}\t${shSetD2PeakNum}\t${shM14PeakNum}\t${shM3PeakNum}\t${shWTAPPeakNum}"

exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shCont/shCont_m6A_rep1-3/peak.xls -o $targetPath/shCont_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shSetD2/shSetD2_m6A_rep1-3/peak.xls -o $targetPath/shSetD2_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shM14/shM14_m6A_rep1-3/peak.xls -o $targetPath/shM14_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shM3/shM3_m6A_rep1-3/peak.xls -o $targetPath/shM3_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shWTAP/shWTAP_m6A_rep1-3/peak.xls -o $targetPath/shWTAP_peaks.bed

shContPeakNum=`cat $targetPath/shCont_peaks.bed | wc -l`
shSetD2PeakNum=`cat $targetPath/shSetD2_peaks.bed | wc -l`
shM14PeakNum=`cat $targetPath/shM14_peaks.bed | wc -l`
shM3PeakNum=`cat $targetPath/shM3_peaks.bed | wc -l`
shWTAPPeakNum=`cat $targetPath/shWTAP_peaks.bed | wc -l`
echo -e "HepG2_rep1-3\t${shContPeakNum}\t${shSetD2PeakNum}\t${shM14PeakNum}\t${shM3PeakNum}\t${shWTAPPeakNum}"

exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shCont/shCont_m6A_rep2-3/peak.xls -o $targetPath/shCont_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shSetD2/shSetD2_m6A_rep2-3/peak.xls -o $targetPath/shSetD2_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shM14/shM14_m6A_rep2-3/peak.xls -o $targetPath/shM14_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shM3/shM3_m6A_rep2-3/peak.xls -o $targetPath/shM3_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shWTAP/shWTAP_m6A_rep2-3/peak.xls -o $targetPath/shWTAP_peaks.bed

shContPeakNum=`cat $targetPath/shCont_peaks.bed | wc -l`
shSetD2PeakNum=`cat $targetPath/shSetD2_peaks.bed | wc -l`
shM14PeakNum=`cat $targetPath/shM14_peaks.bed | wc -l`
shM3PeakNum=`cat $targetPath/shM3_peaks.bed | wc -l`
shWTAPPeakNum=`cat $targetPath/shWTAP_peaks.bed | wc -l`
echo -e "HepG2_rep2-3\t${shContPeakNum}\t${shSetD2PeakNum}\t${shM14PeakNum}\t${shM3PeakNum}\t${shWTAPPeakNum}"

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
echo -e "HepG2_rep1-2-3\t${shContPeakNum}\t${shSetD2PeakNum}\t${shM14PeakNum}\t${shM3PeakNum}\t${shWTAPPeakNum}"

### Hela, m6A-seq
cd /data/zhoukr/hhl_setd2_m6a/Hela_m6A-seq
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shCont/shCont_m6A_rep1/peak.xls -o $targetPath/shCont_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shSetD2/shSetD2_m6A_rep1/peak.xls -o $targetPath/shSetD2_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shM14/shM14_m6A_rep1/peak.xls -o $targetPath/shM14_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shM3/shM3_m6A_rep1/peak.xls -o $targetPath/shM3_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shWTAP/shWTAP_m6A_rep1/peak.xls -o $targetPath/shWTAP_peaks.bed

shContPeakNum=`cat $targetPath/shCont_peaks.bed | wc -l`
shSetD2PeakNum=`cat $targetPath/shSetD2_peaks.bed | wc -l`
shM14PeakNum=`cat $targetPath/shM14_peaks.bed | wc -l`
shM3PeakNum=`cat $targetPath/shM3_peaks.bed | wc -l`
shWTAPPeakNum=`cat $targetPath/shWTAP_peaks.bed | wc -l`
echo -e "Hela_rep1\t${shContPeakNum}\t${shSetD2PeakNum}\t${shM14PeakNum}\t${shM3PeakNum}\t${shWTAPPeakNum}"

exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shCont/shCont_m6A_rep2/peak.xls -o $targetPath/shCont_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shSetD2/shSetD2_m6A_rep2/peak.xls -o $targetPath/shSetD2_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shM14/shM14_m6A_rep2/peak.xls -o $targetPath/shM14_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shM3/shM3_m6A_rep2/peak.xls -o $targetPath/shM3_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shWTAP/shWTAP_m6A_rep2/peak.xls -o $targetPath/shWTAP_peaks.bed

shContPeakNum=`cat $targetPath/shCont_peaks.bed | wc -l`
shSetD2PeakNum=`cat $targetPath/shSetD2_peaks.bed | wc -l`
shM14PeakNum=`cat $targetPath/shM14_peaks.bed | wc -l`
shM3PeakNum=`cat $targetPath/shM3_peaks.bed | wc -l`
shWTAPPeakNum=`cat $targetPath/shWTAP_peaks.bed | wc -l`
echo -e "Hela_rep2\t${shContPeakNum}\t${shSetD2PeakNum}\t${shM14PeakNum}\t${shM3PeakNum}\t${shWTAPPeakNum}"

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
echo -e "Hela_rep3\t${shContPeakNum}\t${shSetD2PeakNum}\t${shM14PeakNum}\t${shM3PeakNum}\t${shWTAPPeakNum}"

exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shCont/shCont_m6A_rep1-2/peak.xls -o $targetPath/shCont_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shSetD2/shSetD2_m6A_rep1-2/peak.xls -o $targetPath/shSetD2_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shM14/shM14_m6A_rep1-2/peak.xls -o $targetPath/shM14_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shM3/shM3_m6A_rep1-2/peak.xls -o $targetPath/shM3_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shWTAP/shWTAP_m6A_rep1-2/peak.xls -o $targetPath/shWTAP_peaks.bed

shContPeakNum=`cat $targetPath/shCont_peaks.bed | wc -l`
shSetD2PeakNum=`cat $targetPath/shSetD2_peaks.bed | wc -l`
shM14PeakNum=`cat $targetPath/shM14_peaks.bed | wc -l`
shM3PeakNum=`cat $targetPath/shM3_peaks.bed | wc -l`
shWTAPPeakNum=`cat $targetPath/shWTAP_peaks.bed | wc -l`
echo -e "Hela_rep1-2\t${shContPeakNum}\t${shSetD2PeakNum}\t${shM14PeakNum}\t${shM3PeakNum}\t${shWTAPPeakNum}"

exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shCont/shCont_m6A_rep1-3/peak.xls -o $targetPath/shCont_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shSetD2/shSetD2_m6A_rep1-3/peak.xls -o $targetPath/shSetD2_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shM14/shM14_m6A_rep1-3/peak.xls -o $targetPath/shM14_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shM3/shM3_m6A_rep1-3/peak.xls -o $targetPath/shM3_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shWTAP/shWTAP_m6A_rep1-3/peak.xls -o $targetPath/shWTAP_peaks.bed

shContPeakNum=`cat $targetPath/shCont_peaks.bed | wc -l`
shSetD2PeakNum=`cat $targetPath/shSetD2_peaks.bed | wc -l`
shM14PeakNum=`cat $targetPath/shM14_peaks.bed | wc -l`
shM3PeakNum=`cat $targetPath/shM3_peaks.bed | wc -l`
shWTAPPeakNum=`cat $targetPath/shWTAP_peaks.bed | wc -l`
echo -e "Hela_rep1-3\t${shContPeakNum}\t${shSetD2PeakNum}\t${shM14PeakNum}\t${shM3PeakNum}\t${shWTAPPeakNum}"

exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shCont/shCont_m6A_rep2-3/peak.xls -o $targetPath/shCont_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shSetD2/shSetD2_m6A_rep2-3/peak.xls -o $targetPath/shSetD2_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shM14/shM14_m6A_rep2-3/peak.xls -o $targetPath/shM14_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shM3/shM3_m6A_rep2-3/peak.xls -o $targetPath/shM3_peaks.bed
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shWTAP/shWTAP_m6A_rep2-3/peak.xls -o $targetPath/shWTAP_peaks.bed

shContPeakNum=`cat $targetPath/shCont_peaks.bed | wc -l`
shSetD2PeakNum=`cat $targetPath/shSetD2_peaks.bed | wc -l`
shM14PeakNum=`cat $targetPath/shM14_peaks.bed | wc -l`
shM3PeakNum=`cat $targetPath/shM3_peaks.bed | wc -l`
shWTAPPeakNum=`cat $targetPath/shWTAP_peaks.bed | wc -l`
echo -e "Hela_rep2-3\t${shContPeakNum}\t${shSetD2PeakNum}\t${shM14PeakNum}\t${shM3PeakNum}\t${shWTAPPeakNum}"

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

echo -e "Hela_rep1-2-3\t${shContPeakNum}\t${shSetD2PeakNum}\t${shM14PeakNum}\t${shM3PeakNum}\t${shWTAPPeakNum}"

exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x shCont/shCont_m6A_rep1-1-3/peak.xls -o $targetPath/shCont_peaks.bed
shContPeakNum=`cat $targetPath/shCont_peaks.bed | wc -l`
echo -e "Hela_rep1-1-3\t${shContPeakNum}\tNA\tNA\tNA\tNA"

exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x $heHelaPath/shCont_new/shCont_m6A_rep1-2/peak.xls -o $targetPath/shCont_peaks.bed
shContPeakNum=`cat $targetPath/shCont_peaks.bed | wc -l`
echo -e "Hechuan_Hela_new_Rep1-2\t${shContPeakNum}\tNA\tNA\tNA\tNA"
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x $heHelaPath/shCont_old/shCont_m6A_rep1-2/peak.xls -o $targetPath/shCont_peaks.bed
shContPeakNum=`cat $targetPath/shCont_peaks.bed | wc -l`
echo -e "Hechuan_Hela_new_Rep1-2\t${shContPeakNum}\tNA\tNA\tNA\tNA"
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x $heHelaPath/shCont/shCont_m6A_rep1-2-3-4/peak.xls -o $targetPath/shCont_peaks.bed
shContPeakNum=`cat $targetPath/shCont_peaks.bed | wc -l`
echo -e "Hechuan_Hela_old_Rep1-2\t${shContPeakNum}\tNA\tNA\tNA\tNA"


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
echo -e "mMEFs\t${shContPeakNum}\t${shSetD2PeakNum}\tNA\tNA\tNA"


cd $targetPath/
rm -f ./*.bed
