#!/bin/sh

cd /data/zhoukr/hhl_setd2_m6a/mouse_MEFs_ChIP-seq
xlsPath=/data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/mMEFs/xls
## xls
rm -rf $xlsPath
mkdir $xlsPath
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -original -x SetD2-WT_NA/macs1.4-merge/shCont_peaks.xls -o $xlsPath/mMEFs_shCont_chip.xls
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -original -x SetD2-KO_NA/macs1.4-merge/shSetD2_peaks.xls -o $xlsPath/mMEFs_shSetD2_chip.xls

### peak fold_enrichment
cd $xlsPath
mkdir allPeak
# all peak
macsPeakFC.pl -cutoff 0.1 -FE -x mMEFs_shCont_chip.xls mMEFs_shSetD2_chip.xls -o ./allPeak/mMEFs_chip_peak_FoldEnrichment.txt
