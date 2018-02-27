#!/bin/sh

cd /data/zhoukr/hhl_setd2_m6a/Hela_m6A-seq

xlsPath=/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/Hela/scatter
expPath=/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/Hela/RSEM
## xls
rm -rf $xlsPath
mkdir $xlsPath
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -format xls -x shCont/shCont_m6A_rep3/peak.xls -o $xlsPath/Hela_shCont_m6a.xls
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -format xls -x shSetD2/shSetD2_m6A_rep3/peak.xls -o $xlsPath/Hela_shSetD2_m6a.xls
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -format xls -x shM14/shM14_m6A_rep3/peak.xls -o $xlsPath/Hela_shM14_m6a.xls
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -format xls -x shM3/shM3_m6A_rep3/peak.xls -o $xlsPath/Hela_shM3_m6a.xls
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -format xls -x shWTAP/shWTAP_m6A_rep3/peak.xls -o $xlsPath/Hela_shWTAP_m6a.xls

### peak fold_enrichment
cd $xlsPath
# all peak
exome2PeaksFC.pl -cutoff 0.1 -nonOverlapA -nonOverlapB -FE -x Hela_shCont_m6a.xls Hela_shSetD2_m6a.xls \
  -o ./Hela_m6a_peak_shSetD2_FoldEnrichment.txt
sed -i '1i chr\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont\tshSetD2' ./Hela_m6a_peak_shSetD2_FoldEnrichment.txt

exome2PeaksFC.pl -cutoff 0.1 -nonOverlapA -nonOverlapB -FE -x Hela_shCont_m6a.xls Hela_shM14_m6a.xls \
  -o ./Hela_m6a_peak_shM14_FoldEnrichment.txt
sed -i '1i chr\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont\tshM14' ./Hela_m6a_peak_shM14_FoldEnrichment.txt

exome2PeaksFC.pl -cutoff 0.1 -nonOverlapA -nonOverlapB -FE -x Hela_shCont_m6a.xls Hela_shM3_m6a.xls \
  -o ./Hela_m6a_peak_shM3_FoldEnrichment.txt
sed -i '1i chr\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont\tshM3' ./Hela_m6a_peak_shM3_FoldEnrichment.txt

exome2PeaksFC.pl -cutoff 0.1 -nonOverlapA -nonOverlapB -FE -x Hela_shCont_m6a.xls Hela_shWTAP_m6a.xls \
  -o ./Hela_m6a_peak_shWTAP_FoldEnrichment.txt
sed -i '1i chr\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont\tshWTAP' ./Hela_m6a_peak_shWTAP_FoldEnrichment.txt


