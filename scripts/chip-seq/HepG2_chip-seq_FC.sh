#!/bin/sh

cd /data/zhoukr/hhl_setd2_m6a/HepG2_ChIP-seq

## xls
rm -rf /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/xls
mkdir /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/xls
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -original -x shCont/macs1.4-merge/shCont_peaks.xls -o /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/xls/HepG2_shCont_chip.xls
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -original -x shSetD2/macs1.4-merge/shSetD2_peaks.xls -o /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/xls/HepG2_shSetD2_chip.xls
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -original -x shM14/macs1.4/shM14_rep1_peaks.xls -o /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/xls/HepG2_shM14_chip.xls
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -original -x shM3/macs1.4/shM3_rep1_peaks.xls -o /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/xls/HepG2_shM3_chip.xls
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -original -x shWTAP/macs1.4/shWTAP_rep1_peaks.xls -o /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/xls/HepG2_shWTAP_chip.xls

### peak fold_enrichment
cd /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/xls
mkdir allPeak
mkdir overlappedPeak
# all peak
macsPeakFC.pl -cutoff 0.1 -FE -x HepG2_shCont_chip.xls HepG2_shSetD2_chip.xls HepG2_shM14_chip.xls HepG2_shM3_chip.xls HepG2_shWTAP_chip.xls -o ./allPeak/HepG2_chip_peak_FoldEnrichment.txt
macsPeakFC.pl -cutoff 0.1 -FE -F 0.5 -x HepG2_shCont_chip.xls HepG2_shSetD2_chip.xls HepG2_shM14_chip.xls HepG2_shM3_chip.xls HepG2_shWTAP_chip.xls -o ./allPeak/HepG2_chip_peak_-F0.5_FoldEnrichment.txt
macsPeakFC.pl -cutoff 0.1 -FE -log 2 -x HepG2_shCont_chip.xls HepG2_shSetD2_chip.xls HepG2_shM14_chip.xls HepG2_shM3_chip.xls HepG2_shWTAP_chip.xls -o ./allPeak/HepG2_chip_peak_FoldEnrichment_log2.txt
macsPeakFC.pl -cutoff 0.1 -FE -log 2 -F 0.5 -x HepG2_shCont_chip.xls HepG2_shSetD2_chip.xls HepG2_shM14_chip.xls HepG2_shM3_chip.xls HepG2_shWTAP_chip.xls -o ./allPeak/HepG2_chip_peak_-F0.5_FoldEnrichment_log2.txt
macsPeakFC.pl -cutoff 0.1 -x HepG2_shCont_chip.xls HepG2_shSetD2_chip.xls HepG2_shM14_chip.xls HepG2_shM3_chip.xls HepG2_shWTAP_chip.xls -o ./allPeak/HepG2_chip_peak_FC.txt
macsPeakFC.pl -cutoff 0.1 -F 0.5 -x HepG2_shCont_chip.xls HepG2_shSetD2_chip.xls HepG2_shM14_chip.xls HepG2_shM3_chip.xls HepG2_shWTAP_chip.xls -o ./allPeak/HepG2_chip_peak_-F0.5_FC.txt
macsPeakFC.pl -cutoff 0.1 -log 2 -x HepG2_shCont_chip.xls HepG2_shSetD2_chip.xls HepG2_shM14_chip.xls HepG2_shM3_chip.xls HepG2_shWTAP_chip.xls -o ./allPeak/HepG2_chip_peak_FC_log2.txt
macsPeakFC.pl -cutoff 0.1 -log 2 -F 0.5 -x HepG2_shCont_chip.xls HepG2_shSetD2_chip.xls HepG2_shM14_chip.xls HepG2_shM3_chip.xls HepG2_shWTAP_chip.xls -o ./allPeak/HepG2_chip_peak_-F0.5_FC_log2.txt

# overlapped peak
macsPeakFC.pl -cutoff 0.1 -overlap 0 1 1 1 1 -FE -x HepG2_shCont_chip.xls HepG2_shSetD2_chip.xls HepG2_shM14_chip.xls HepG2_shM3_chip.xls HepG2_shWTAP_chip.xls -o ./overlappedPeak/HepG2_chip_peak_FoldEnrichment.txt
macsPeakFC.pl -cutoff 0.1 -overlap 0 1 1 1 1 -FE -F 0.5 -x HepG2_shCont_chip.xls HepG2_shSetD2_chip.xls HepG2_shM14_chip.xls HepG2_shM3_chip.xls HepG2_shWTAP_chip.xls -o ./overlappedPeak/HepG2_chip_peak_-F0.5_FoldEnrichment.txt
macsPeakFC.pl -cutoff 0.1 -overlap 0 1 1 1 1 -FE -log 2 -x HepG2_shCont_chip.xls HepG2_shSetD2_chip.xls HepG2_shM14_chip.xls HepG2_shM3_chip.xls HepG2_shWTAP_chip.xls -o ./overlappedPeak/HepG2_chip_peak_FoldEnrichment_log2.txt
macsPeakFC.pl -cutoff 0.1 -overlap 0 1 1 1 1 -FE -log 2 -F 0.5 -x HepG2_shCont_chip.xls HepG2_shSetD2_chip.xls HepG2_shM14_chip.xls HepG2_shM3_chip.xls HepG2_shWTAP_chip.xls -o ./overlappedPeak/HepG2_chip_peak_-F0.5_FoldEnrichment_log2.txt
macsPeakFC.pl -cutoff 0.1 -overlap 0 1 1 1 1 -x HepG2_shCont_chip.xls HepG2_shSetD2_chip.xls HepG2_shM14_chip.xls HepG2_shM3_chip.xls HepG2_shWTAP_chip.xls -o ./overlappedPeak/HepG2_chip_peak_FC.txt
macsPeakFC.pl -cutoff 0.1 -overlap 0 1 1 1 1 -F 0.5 -x HepG2_shCont_chip.xls HepG2_shSetD2_chip.xls HepG2_shM14_chip.xls HepG2_shM3_chip.xls HepG2_shWTAP_chip.xls -o ./overlappedPeak/HepG2_chip_peak_-F0.5_FC.txt
macsPeakFC.pl -cutoff 0.1 -overlap 0 1 1 1 1 -log 2 -x HepG2_shCont_chip.xls HepG2_shSetD2_chip.xls HepG2_shM14_chip.xls HepG2_shM3_chip.xls HepG2_shWTAP_chip.xls -o ./overlappedPeak/HepG2_chip_peak_FC_log2.txt
macsPeakFC.pl -cutoff 0.1 -overlap 0 1 1 1 1 -log 2 -F 0.5 -x HepG2_shCont_chip.xls HepG2_shSetD2_chip.xls HepG2_shM14_chip.xls HepG2_shM3_chip.xls HepG2_shWTAP_chip.xls -o ./overlappedPeak/HepG2_chip_peak_-F0.5_FC_log2.txt
