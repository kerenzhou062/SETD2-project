#!/bin/sh

xlsPath=/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/mouse/scatter
rm -rf $xlsPath
mkdir -p $xlsPath

cd /data/zhoukr/hhl_setd2_m6a/mouse_ESCs_m6A-seq
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -format xls -x Ctrl_D0/Ctrl_D0_m6A/peak.xls -o $xlsPath/mESCs_shCont_D0_m6a.xls
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -format xls -x Ctrl_D6/Ctrl_D6_m6A/peak.xls -o $xlsPath/mESCs_shCont_D6_m6a.xls
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -format xls -x SetD2-KD_D0/SetD2-KD_D0_m6A/peak.xls -o $xlsPath/mESCs_shSetD2_D0_m6a.xls
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -format xls -x SetD2-KD_D6/SetD2-KD_D6_m6A/peak.xls -o $xlsPath/mESCs_shSetD2_D6_m6a.xls

cd /data/zhoukr/hhl_setd2_m6a/mouse_MEFs_m6A-seq
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -format xls -x SetD2-WT_NA/SetD2-WT_NA_m6A/peak.xls -o $xlsPath/mMEFs_shCont_m6a.xls
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -format xls -x SetD2-KO_NA/SetD2-KO_NA_m6A/peak.xls -o $xlsPath/mMEFs_shSetD2_m6a.xls
### peak fold_enrichment
cd $xlsPath
# all peak
exome2PeaksFC.pl -cutoff 0.1 -nonOverlapA -nonOverlapB -FE -x mESCs_shCont_D0_m6a.xls mESCs_shCont_D6_m6a.xls \
  -o ./mESCs_m6a_peak_shContD6vsD0_FoldEnrichment.txt
sed -i '1i chr\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont-D0\tshCont-D6' ./mESCs_m6a_peak_shContD6vsD0_FoldEnrichment.txt

exome2PeaksFC.pl -cutoff 0.1 -nonOverlapA -nonOverlapB -FE -x mESCs_shCont_D0_m6a.xls mESCs_shSetD2_D6_m6a.xls \
  -o ./mESCs_m6a_peak_shSetD2D6vsShContD0_FoldEnrichment.txt
sed -i '1i chr\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont-D0\tshSetD2-D6' ./mESCs_m6a_peak_shContD6vsD0_FoldEnrichment.txt

exome2PeaksFC.pl -cutoff 0.1 -nonOverlapA -nonOverlapB -FE -x mESCs_shCont_D0_m6a.xls mESCs_shSetD2_D6_m6a.xls \
  -o ./mESCs_m6a_peak_shSetD2D6vsD0_FoldEnrichment.txt
sed -i '1i chr\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshSetD2-D0\tshSetD2-D6' ./mESCs_m6a_peak_shSetD2D6vsD0_FoldEnrichment.txt

exome2PeaksFC.pl -cutoff 0.1 -nonOverlapA -nonOverlapB -FE -x mESCs_shCont_D6_m6a.xls mESCs_shSetD2_D6_m6a.xls \
  -o ./mESCs_m6a_peak_shSetD2D6vsShContD6_FoldEnrichment.txt
sed -i '1i chr\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont-D6\tshSetD2-D6' ./mESCs_m6a_peak_shSetD2D6vsShContD6_FoldEnrichment.txt

exome2PeaksFC.pl -cutoff 0.1 -nonOverlapA -nonOverlapB -FE -x mMEFs_shCont_m6a.xls mMEFs_shSetD2_m6a.xls \
  -o ./mEFs_m6a_peak_shSetD2vsShCont_FoldEnrichment.txt
sed -i '1i chr\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont\tshSetD2' ./mEFs_m6a_peak_shSetD2vsShCont_FoldEnrichment.txt
