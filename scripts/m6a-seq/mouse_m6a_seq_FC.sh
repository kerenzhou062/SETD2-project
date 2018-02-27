#!/bin/sh
mouseXlsPath=/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/mouse/xls

cd /data/zhoukr/hhl_setd2_m6a/mouse_ESCs_m6A-seq
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 0 -original -format xls -x Ctrl_D0/Ctrl_D0_m6A/peak.xls -o $mouseXlsPath/mESCs_shCont_D0_m6a.xls
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 0 -original -format xls -x Ctrl_D6/Ctrl_D6_m6A/peak.xls -o $mouseXlsPath/mESCs_shCont_D6_m6a.xls
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 0 -original -format xls -x SetD2-KD_D0/SetD2-KD_D0_m6A/peak.xls -o $mouseXlsPath/mESCs_shSetD2_D0_m6a.xls
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 0 -original -format xls -x SetD2-KD_D6/SetD2-KD_D6_m6A/peak.xls -o $mouseXlsPath/mESCs_shSetD2_D6_m6a.xls

cd /data/zhoukr/hhl_setd2_m6a/mouse_MEFs_m6A-seq
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 0 -original -format xls -x SetD2-WT_NA/SetD2-WT_NA_m6A/peak.xls -o $mouseXlsPath/mMEFs_shCont_m6a.xls
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 0 -original -format xls -x SetD2-KO_NA/SetD2-KO_NA_m6A/peak.xls -o $mouseXlsPath/mMEFs_shSetD2_m6a.xls

cd /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/mouse/xls

cd $mouseXlsPath
rm -rf overlappedPeak
mkdir overlappedPeak

# overlapped peak

exomePeakFC.pl -cutoff 0.1 -overlap 0 1 1 1 1 1 1 -FE -x mESCs_shCont_D0_m6a.xls mESCs_shCont_D0_m6a.xls mESCs_shCont_D6_m6a.xls \
  mESCs_shSetD2_D0_m6a.xls mESCs_shSetD2_D6_m6a.xls \
  mMEFs_shCont_m6a.xls mMEFs_shSetD2_m6a.xls -o ./temp.tmp
awk 'BEGIN{OFS="\t";FS="\t"; bed12="chr\tchromStart\tchromEnd\tid\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts";
  append="mESCs_shCont_D0_ref\tmESCs-shCont-D0\tmESCs-shCont-D6\tmESCs-shSETD2-D0\tmESCs-shSETD2-D6\tmMEFs-shCont\tmMEFs-shSETD2";
  print bed12"\t"append; }
ARGIND==1{split($4,geneIdI,":");split(geneIdI[1],geneIdInfo,".");geneIdi=geneIdInfo[1];geneNamei=geneIdI[2];hashGeneIdName[geneIdi]=geneNamei;}
ARGIND==2{split($4,geneIdH,".");geneIdh=geneIdH[1]; if (geneIdh in hashGeneIdName){geneNameh=hashGeneIdName[geneIdh]}else{geneNameh="-"};$4=$4"\t"geneNameh; print $0}' \
/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.vM9.annotation.all.longest.exon.bed6 ./temp.tmp > ./overlappedPeak/mouse_all_samples_m6a_peak_FE.txt


exomePeakFC.pl -cutoff 0.1 -log 2 -overlap 0 1 1 1 1 1 1 -FE -x mESCs_shCont_D0_m6a.xls mESCs_shCont_D0_m6a.xls mESCs_shCont_D6_m6a.xls \
  mESCs_shSetD2_D0_m6a.xls mESCs_shSetD2_D6_m6a.xls \
  mMEFs_shCont_m6a.xls mMEFs_shSetD2_m6a.xls -o ./temp.tmp
awk 'BEGIN{OFS="\t";FS="\t"; bed12="chr\tchromStart\tchromEnd\tid\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts";
  append="mESCs_shCont_D0_ref\tmESCs-shCont-D0\tmESCs-shCont-D6\tmESCs-shSETD2-D0\tmESCs-shSETD2-D6\tmMEFs-shCont\tmMEFs-shSETD2";
  print bed12"\t"append; }
ARGIND==1{split($4,geneIdI,":");split(geneIdI[1],geneIdInfo,".");geneIdi=geneIdInfo[1];geneNamei=geneIdI[2];hashGeneIdName[geneIdi]=geneNamei;}
ARGIND==2{split($4,geneIdH,".");geneIdh=geneIdH[1]; if (geneIdh in hashGeneIdName){geneNameh=hashGeneIdName[geneIdh]}else{geneNameh="-"};$4=$4"\t"geneNameh; print $0}' \
/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.vM9.annotation.all.longest.exon.bed6 ./temp.tmp > ./overlappedPeak/mouse_all_samples_m6a_peak_FE_log2.txt

rm -f ./temp.tmp
