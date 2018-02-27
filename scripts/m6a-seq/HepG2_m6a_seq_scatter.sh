#!/bin/sh

cd /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq

xlsPath=/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/scatter
expPath=/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/HepG2/RSEM
## xls
rm -rf $xlsPath
mkdir $xlsPath
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -format xls -x shCont/shCont_m6A_rep1-2-3/peak.xls -o $xlsPath/HepG2_shCont_m6a.xls
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -format xls -x shSetD2/shSetD2_m6A_rep1-2-3/peak.xls -o $xlsPath/HepG2_shSetD2_m6a.xls
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -format xls -x shM14/shM14_m6A_rep1-2-3/peak.xls -o $xlsPath/HepG2_shM14_m6a.xls
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -format xls -x shM3/shM3_m6A_rep1-2-3/peak.xls -o $xlsPath/HepG2_shM3_m6a.xls
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -format xls -x shWTAP/shWTAP_m6A_rep1-2-3/peak.xls -o $xlsPath/HepG2_shWTAP_m6a.xls

### peak fold_enrichment
cd $xlsPath
# all peak
exome2PeaksFC.pl -cutoff 0.1 -nonOverlapA -nonOverlapB -FE -x HepG2_shCont_m6a.xls HepG2_shSetD2_m6a.xls \
  -o ./HepG2_m6a_peak_shSetD2_FoldEnrichment.txt
sed -i '1i chr\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont\tshSetD2' ./HepG2_m6a_peak_shSetD2_FoldEnrichment.txt

exome2PeaksFC.pl -cutoff 0.1 -nonOverlapA -nonOverlapB -FE -x HepG2_shCont_m6a.xls HepG2_shM14_m6a.xls \
  -o ./HepG2_m6a_peak_shM14_FoldEnrichment.txt
sed -i '1i chr\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont\tshM14' ./HepG2_m6a_peak_shM14_FoldEnrichment.txt

exome2PeaksFC.pl -cutoff 0.1 -nonOverlapA -nonOverlapB -FE -x HepG2_shCont_m6a.xls HepG2_shM3_m6a.xls \
  -o ./HepG2_m6a_peak_shM3_FoldEnrichment.txt
sed -i '1i chr\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont\tshM3' ./HepG2_m6a_peak_shM3_FoldEnrichment.txt

exome2PeaksFC.pl -cutoff 0.1 -nonOverlapA -nonOverlapB -FE -x HepG2_shCont_m6a.xls HepG2_shWTAP_m6a.xls \
  -o ./HepG2_m6a_peak_shWTAP_FoldEnrichment.txt
sed -i '1i chr\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont\tshWTAP' ./HepG2_m6a_peak_shWTAP_FoldEnrichment.txt

mkdir cumulative
cd cumulative

awk '{print log($15+1)/log(2)}' ../HepG2_shCont_m6a.xls > HepG2_shCont_m6a.txt
awk '{print log($15+1)/log(2)}' ../HepG2_shSetD2_m6a.xls > HepG2_shSetD2_m6a.txt

cumulativePlot.pl -header false -interval 0.01 --input HepG2_shCont_m6a.txt HepG2_shSetD2_m6a.txt -o HepG2_m6A_shCont_shSetD2_cumulativePlot.txt
