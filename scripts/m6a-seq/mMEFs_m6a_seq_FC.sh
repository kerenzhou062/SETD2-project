#!/bin/sh
xlsPath=/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/mMEFs/xls
cd /data/zhoukr/hhl_setd2_m6a/mouse_MEFs_m6A-seq

## xls
rm -rf $xlsPath
mkdir $xlsPath
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -format xls -x SetD2-WT_NA/SetD2-WT_NA_m6A/peak.xls -o $xlsPath/mMEFs_shCont_m6a.xls
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -format xls -x SetD2-KO_NA/SetD2-KO_NA_m6A/peak.xls -o $xlsPath/mMEFs_shSetD2_m6a.xls

### peak fold_enrichment
cd /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/mMEFs/xls
mkdir allPeak


exomePeakFC.pl -cutoff 0.1 -FE -x mMEFs_shCont_m6a.xls mMEFs_shSetD2_m6a.xls -o ./temp.tmp
awk 'BEGIN{OFS="\t";FS="\t"; print "chr\tchromStart\tchromEnd\tid\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont\tshSetD2"; }
ARGIND==1{split($4,geneIdI,":");split(geneIdI[1],geneIdInfo,".");geneIdi=geneIdInfo[1];geneNamei=geneIdI[2];hashGeneIdName[geneIdi]=geneNamei;}
ARGIND==2{split($4,geneIdH,".");geneIdh=geneIdH[1]; if (geneIdh in hashGeneIdName){geneNameh=hashGeneIdName[geneIdh]}else{geneNameh="-"};$4=$4"\t"geneNameh; print $0}' \
/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.vM9.annotation.all.longest.exon.bed6 ./temp.tmp > ./allPeak/mMEFs_shCont_vs_shSetD2_m6a_peak_FE.txt

exomePeakFC.pl -cutoff 0.1 -FE -log 2 -x mMEFs_shCont_m6a.xls mMEFs_shSetD2_m6a.xls -o ./temp.tmp
awk 'BEGIN{OFS="\t";FS="\t"; print "chr\tchromStart\tchromEnd\tid\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont\tshSetD2"; }
ARGIND==1{split($4,geneIdI,":");split(geneIdI[1],geneIdInfo,".");geneIdi=geneIdInfo[1];geneNamei=geneIdI[2];hashGeneIdName[geneIdi]=geneNamei;}
ARGIND==2{split($4,geneIdH,".");geneIdh=geneIdH[1]; if (geneIdh in hashGeneIdName){geneNameh=hashGeneIdName[geneIdh]}else{geneNameh="-"};$4=$4"\t"geneNameh; print $0}' \
/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.vM9.annotation.all.longest.exon.bed6 ./temp.tmp > ./allPeak/mMEFs_shCont_vs_shSetD2_m6a_peak_FE_log2.txt

exomePeakFC.pl -cutoff 0.1 -x mMEFs_shCont_m6a.xls mMEFs_shSetD2_m6a.xls -o ./temp.tmp
awk 'BEGIN{OFS="\t";FS="\t"; print "chr\tchromStart\tchromEnd\tid\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont\tshSetD2"; }
ARGIND==1{split($4,geneIdI,":");split(geneIdI[1],geneIdInfo,".");geneIdi=geneIdInfo[1];geneNamei=geneIdI[2];hashGeneIdName[geneIdi]=geneNamei;}
ARGIND==2{split($4,geneIdH,".");geneIdh=geneIdH[1]; if (geneIdh in hashGeneIdName){geneNameh=hashGeneIdName[geneIdh]}else{geneNameh="-"};$4=$4"\t"geneNameh; print $0}' \
/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.vM9.annotation.all.longest.exon.bed6 ./temp.tmp > ./allPeak/mMEFs_shCont_vs_shSetD2_m6a_peak_FC.txt

exomePeakFC.pl -cutoff 0.1 -log 2 -x mMEFs_shCont_m6a.xls mMEFs_shSetD2_m6a.xls -o ./temp.tmp
awk 'BEGIN{OFS="\t";FS="\t"; print "chr\tchromStart\tchromEnd\tid\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont\tshSetD2"; }
ARGIND==1{split($4,geneIdI,":");split(geneIdI[1],geneIdInfo,".");geneIdi=geneIdInfo[1];geneNamei=geneIdI[2];hashGeneIdName[geneIdi]=geneNamei;}
ARGIND==2{split($4,geneIdH,".");geneIdh=geneIdH[1]; if (geneIdh in hashGeneIdName){geneNameh=hashGeneIdName[geneIdh]}else{geneNameh="-"};$4=$4"\t"geneNameh; print $0}' \
/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.vM9.annotation.all.longest.exon.bed6 ./temp.tmp > ./allPeak/mMEFs_shCont_vs_shSetD2_m6a_peak_FC_log2.txt


rm -f ./temp.tmp
