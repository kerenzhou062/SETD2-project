#!/bin/sh

cd /data/zhoukr/hhl_setd2_m6a/Hela_m6A-seq

xlsPath=/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/Hela/boxplot/xls
## xls
rm -rf $xlsPath
mkdir -p $xlsPath
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 0 -original -format xls -x shCont/shCont_m6A_rep3/peak.xls -o $xlsPath/Hela_shCont_m6a.xls
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 0 -original -format xls -x shSetD2/shSetD2_m6A_rep3/peak.xls -o $xlsPath/Hela_shSetD2_m6a.xls
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 0 -original -format xls -x shM14/shM14_m6A_rep3/peak.xls -o $xlsPath/Hela_shM14_m6a.xls
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 0 -original -format xls -x shM3/shM3_m6A_rep3/peak.xls -o $xlsPath/Hela_shM3_m6a.xls
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 0 -original -format xls -x shWTAP/shWTAP_m6A_rep3/peak.xls -o $xlsPath/Hela_shWTAP_m6a.xls

### peak fold_enrichment
cd $xlsPath
# all peak
exomePeakFC.pl -cutoff 0.1 -FE -x Hela_shCont_m6a.xls Hela_shCont_m6a.xls Hela_shSetD2_m6a.xls Hela_shM14_m6a.xls Hela_shM3_m6a.xls Hela_shWTAP_m6a.xls \
  -o ./allPeak/Hela_m6a_peak_FoldEnrichment.txt
sed -i '1i chr\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont_ref\tshCont\tshSetD2\tshM14\tshM3\tshWTAP' ./allPeak/Hela_m6a_peak_FoldEnrichment.txt

exomePeakFC.pl -cutoff 0.1 -FE -log 2 -x Hela_shCont_m6a.xls Hela_shCont_m6a.xls Hela_shSetD2_m6a.xls Hela_shM14_m6a.xls Hela_shM3_m6a.xls Hela_shWTAP_m6a.xls \
  -o ./allPeak/Hela_m6a_peak_FoldEnrichment_log2.txt
sed -i '1i chr\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont_ref\tshCont\tshSetD2\tshM14\tshM3\tshWTAP' ./allPeak/Hela_m6a_peak_FoldEnrichment_log2.txt


# overlapped peak
exomePeakFC.pl -cutoff 0.1 -overlap 0 1 1 1 1 1 -FE -x Hela_shCont_m6a.xls Hela_shCont_m6a.xls Hela_shSetD2_m6a.xls Hela_shM14_m6a.xls Hela_shM3_m6a.xls Hela_shWTAP_m6a.xls \
  -o ./overlappedPeak/Hela_m6a_peak_FoldEnrichment.txt
sed -i '1i chr\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont_ref\tshCont\tshSetD2\tshM14\tshM3\tshWTAP' ./overlappedPeak/Hela_m6a_peak_FoldEnrichment.txt

exomePeakFC.pl -cutoff 0.1 -overlap 0 1 1 1 1 1 -FE -log 2 -x Hela_shCont_m6a.xls Hela_shCont_m6a.xls Hela_shSetD2_m6a.xls Hela_shM14_m6a.xls Hela_shM3_m6a.xls Hela_shWTAP_m6a.xls \
  -o ./overlappedPeak/Hela_m6a_peak_FoldEnrichment_log2.txt
sed -i '1i chr\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont_ref\tshCont\tshSetD2\tshM14\tshM3\tshWTAP' ./overlappedPeak/Hela_m6a_peak_FoldEnrichment_log2.txt
