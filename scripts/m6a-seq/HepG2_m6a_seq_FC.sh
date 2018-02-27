#!/bin/sh

cd /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq

xlsPath=/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/xls
expPath=/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/HepG2/RSEM
## xls
rm -rf $xlsPath
mkdir $xlsPath
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -format xls -x shCont/shCont_m6A_rep1-2-3/peak.xls -o $xlsPath/HepG2_shCont_m6a.xls
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -format xls -x shSetD2/shSetD2_m6A_rep1-2-3/peak.xls -o $xlsPath/HepG2_shSetD2_m6a.xls
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -format xls -x shM14/shM14_m6A_rep1-2-3/peak.xls -o $xlsPath/HepG2_shM14_m6a.xls
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -format xls -x shM3/shM3_m6A_rep1-2-3/peak.xls -o $xlsPath/HepG2_shM3_m6a.xls
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -format xls -x shWTAP/shWTAP_m6A_rep1-2-3/peak.xls -o $xlsPath/HepG2_shWTAP_m6a.xls

##### filter with gene expression
#exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -format xls -x shCont/shCont_m6A_rep1-2-3/peak.xls -o $xlsPath/temp.xls
#exomePeakExpFilter.pl -sample mean -cutoff 1 -expF $expPath/HepG2_RNA-seq_RSEM_shCont_input.txt -peakF $xlsPath/temp.xls -o $xlsPath/HepG2_shCont_m6a.xls
#
#exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -format xls -x shSetD2/shSetD2_m6A_rep1-2-3/peak.xls -o $xlsPath/temp.xls
#exomePeakExpFilter.pl -sample mean -cutoff 1 -expF $expPath/HepG2_RNA-seq_RSEM_shSetD2_input.txt -peakF $xlsPath/temp.xls -o $xlsPath/HepG2_shSetD2_m6a.xls
#
#exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -format xls -x shM14/shM14_m6A_rep1-2-3/peak.xls -o $xlsPath/temp.xls
#exomePeakExpFilter.pl -sample mean -cutoff 1 -expF $expPath/HepG2_RNA-seq_RSEM_shM14_input.txt -peakF $xlsPath/temp.xls -o $xlsPath/HepG2_shM14_m6a.xls
#
#exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -format xls -x shM3/shM3_m6A_rep1-2-3/peak.xls -o $xlsPath/temp.xls
#exomePeakExpFilter.pl -sample mean -cutoff 1 -expF $expPath/HepG2_RNA-seq_RSEM_shM3_input.txt -peakF $xlsPath/temp.xls -o $xlsPath/HepG2_shM3_m6a.xls
#
#exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -format xls -x shWTAP/shWTAP_m6A_rep1-2-3/peak.xls -o $xlsPath/temp.xls
#exomePeakExpFilter.pl -sample mean -cutoff 1 -expF $expPath/HepG2_RNA-seq_RSEM_shWTAP_input.txt -peakF $xlsPath/temp.xls -o $xlsPath/HepG2_shWTAP_m6a.xls
#
#rm -f $xlsPath/temp.xls

### peak fold_enrichment
cd $xlsPath
mkdir allPeak
mkdir overlappedPeak
# all peak
exomePeakFC.pl -cutoff 0.1 -FE -x HepG2_shCont_m6a.xls HepG2_shSetD2_m6a.xls HepG2_shM14_m6a.xls HepG2_shM3_m6a.xls HepG2_shWTAP_m6a.xls \
  -o ./allPeak/HepG2_m6a_peak_FoldEnrichment.txt
sed -i '1i chr\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont\tshSetD2\tshM14\tshM3\tshWTAP' ./allPeak/HepG2_m6a_peak_FoldEnrichment.txt

exomePeakFC.pl -cutoff 0.1 -FE -F 0.5 -x HepG2_shCont_m6a.xls HepG2_shSetD2_m6a.xls HepG2_shM14_m6a.xls HepG2_shM3_m6a.xls HepG2_shWTAP_m6a.xls \
  -o ./allPeak/HepG2_m6a_peak_-F0.5_FoldEnrichment.txt
sed -i '1i chr\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont\tshSetD2\tshM14\tshM3\tshWTAP' ./allPeak/HepG2_m6a_peak_-F0.5_FoldEnrichment.txt

exomePeakFC.pl -cutoff 0.1 -FE -log 2 -x HepG2_shCont_m6a.xls HepG2_shSetD2_m6a.xls HepG2_shM14_m6a.xls HepG2_shM3_m6a.xls HepG2_shWTAP_m6a.xls \
  -o ./allPeak/HepG2_m6a_peak_FoldEnrichment_log2.txt
sed -i '1i chr\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont\tshSetD2\tshM14\tshM3\tshWTAP' ./allPeak/HepG2_m6a_peak_FoldEnrichment_log2.txt

exomePeakFC.pl -cutoff 0.1 -FE -log 2 -F 0.5 -x HepG2_shCont_m6a.xls HepG2_shSetD2_m6a.xls HepG2_shM14_m6a.xls HepG2_shM3_m6a.xls HepG2_shWTAP_m6a.xls \
  -o ./allPeak/HepG2_m6a_peak_-F0.5_FoldEnrichment_log2.txt
sed -i '1i chr\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont\tshSetD2\tshM14\tshM3\tshWTAP' ./allPeak/HepG2_m6a_peak_-F0.5_FoldEnrichment_log2.txt

exomePeakFC.pl -cutoff 0.1 -x HepG2_shCont_m6a.xls HepG2_shSetD2_m6a.xls HepG2_shM14_m6a.xls HepG2_shM3_m6a.xls HepG2_shWTAP_m6a.xls \
  -o ./allPeak/HepG2_m6a_peak_FC.txt
sed -i '1i chr\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont\tshSetD2\tshM14\tshM3\tshWTAP' ./allPeak/HepG2_m6a_peak_FC.txt

exomePeakFC.pl -cutoff 0.1 -F 0.5 -x HepG2_shCont_m6a.xls HepG2_shSetD2_m6a.xls HepG2_shM14_m6a.xls HepG2_shM3_m6a.xls HepG2_shWTAP_m6a.xls \
  -o ./allPeak/HepG2_m6a_peak_-F0.5_FC.txt
sed -i '1i chr\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont\tshSetD2\tshM14\tshM3\tshWTAP' ./allPeak/HepG2_m6a_peak_-F0.5_FC.txt

exomePeakFC.pl -cutoff 0.1 -log 2 -x HepG2_shCont_m6a.xls HepG2_shSetD2_m6a.xls HepG2_shM14_m6a.xls HepG2_shM3_m6a.xls HepG2_shWTAP_m6a.xls \
  -o ./allPeak/HepG2_m6a_peak_FC_log2.txt
sed -i '1i chr\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont\tshSetD2\tshM14\tshM3\tshWTAP' ./allPeak/HepG2_m6a_peak_FC_log2.txt

exomePeakFC.pl -cutoff 0.1 -log 2 -F 0.5 -x HepG2_shCont_m6a.xls HepG2_shSetD2_m6a.xls HepG2_shM14_m6a.xls HepG2_shM3_m6a.xls HepG2_shWTAP_m6a.xls \
  -o ./allPeak/HepG2_m6a_peak_-F0.5_FC_log2.txt
sed -i '1i chr\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont\tshSetD2\tshM14\tshM3\tshWTAP' ./allPeak/HepG2_m6a_peak_-F0.5_FC_log2.txt



# overlapped peak
exomePeakFC.pl -cutoff 0.1 -overlap 0 1 1 1 1 -FE -x HepG2_shCont_m6a.xls HepG2_shSetD2_m6a.xls HepG2_shM14_m6a.xls HepG2_shM3_m6a.xls HepG2_shWTAP_m6a.xls \
  -o ./overlappedPeak/HepG2_m6a_peak_FoldEnrichment.txt
sed -i '1i chr\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont\tshSetD2\tshM14\tshM3\tshWTAP' ./overlappedPeak/HepG2_m6a_peak_FoldEnrichment.txt

exomePeakFC.pl -cutoff 0.1 -overlap 0 1 1 1 1 -FE -F 0.5 -x HepG2_shCont_m6a.xls HepG2_shSetD2_m6a.xls HepG2_shM14_m6a.xls HepG2_shM3_m6a.xls HepG2_shWTAP_m6a.xls \
  -o ./overlappedPeak/HepG2_m6a_peak_-F0.5_FoldEnrichment.txt
sed -i '1i chr\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont\tshSetD2\tshM14\tshM3\tshWTAP' ./overlappedPeak/HepG2_m6a_peak_-F0.5_FoldEnrichment.txt

exomePeakFC.pl -cutoff 0.1 -overlap 0 1 1 1 1 -FE -log 2 -x HepG2_shCont_m6a.xls HepG2_shSetD2_m6a.xls HepG2_shM14_m6a.xls HepG2_shM3_m6a.xls HepG2_shWTAP_m6a.xls \
  -o ./overlappedPeak/HepG2_m6a_peak_FoldEnrichment_log2.txt
sed -i '1i chr\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont\tshSetD2\tshM14\tshM3\tshWTAP' ./overlappedPeak/HepG2_m6a_peak_FoldEnrichment_log2.txt

exomePeakFC.pl -cutoff 0.1 -overlap 0 1 1 1 1 -FE -log 2 -F 0.5 -x HepG2_shCont_m6a.xls HepG2_shSetD2_m6a.xls HepG2_shM14_m6a.xls HepG2_shM3_m6a.xls HepG2_shWTAP_m6a.xls \
  -o ./overlappedPeak/HepG2_m6a_peak_-F0.5_FoldEnrichment_log2.txt
sed -i '1i chr\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont\tshSetD2\tshM14\tshM3\tshWTAP' ./overlappedPeak/HepG2_m6a_peak_-F0.5_FoldEnrichment_log2.txt

exomePeakFC.pl -cutoff 0.1 -overlap 0 1 1 1 1 -x HepG2_shCont_m6a.xls HepG2_shSetD2_m6a.xls HepG2_shM14_m6a.xls HepG2_shM3_m6a.xls HepG2_shWTAP_m6a.xls \
  -o ./overlappedPeak/HepG2_m6a_peak_FC.txt
sed -i '1i chr\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont\tshSetD2\tshM14\tshM3\tshWTAP' ./overlappedPeak/HepG2_m6a_peak_FC.txt

exomePeakFC.pl -cutoff 0.1 -overlap 0 1 1 1 1 -F 0.5 -x HepG2_shCont_m6a.xls HepG2_shSetD2_m6a.xls HepG2_shM14_m6a.xls HepG2_shM3_m6a.xls HepG2_shWTAP_m6a.xls \
  -o ./overlappedPeak/HepG2_m6a_peak_-F0.5_FC.txt
sed -i '1i chr\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont\tshSetD2\tshM14\tshM3\tshWTAP' ./overlappedPeak/HepG2_m6a_peak_-F0.5_FC.txt

exomePeakFC.pl -cutoff 0.1 -overlap 0 1 1 1 1 -log 2 -x HepG2_shCont_m6a.xls HepG2_shSetD2_m6a.xls HepG2_shM14_m6a.xls HepG2_shM3_m6a.xls HepG2_shWTAP_m6a.xls \
  -o ./overlappedPeak/HepG2_m6a_peak_FC_log2.txt
sed -i '1i chr\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont\tshSetD2\tshM14\tshM3\tshWTAP' ./overlappedPeak/HepG2_m6a_peak_FC_log2.txt

exomePeakFC.pl -cutoff 0.1 -overlap 0 1 1 1 1 -log 2 -F 0.5 -x HepG2_shCont_m6a.xls HepG2_shSetD2_m6a.xls HepG2_shM14_m6a.xls HepG2_shM3_m6a.xls HepG2_shWTAP_m6a.xls \
  -o ./overlappedPeak/HepG2_m6a_peak_-F0.5_FC_log2.txt
sed -i '1i chr\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont\tshSetD2\tshM14\tshM3\tshWTAP' ./overlappedPeak/HepG2_m6a_peak_-F0.5_FC_log2.txt

