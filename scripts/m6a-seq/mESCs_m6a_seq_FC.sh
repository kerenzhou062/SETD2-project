#!/bin/sh
xlsPath=/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/mESCs/xls
cd /data/zhoukr/hhl_setd2_m6a/mouse_ESCs_m6A-seq

## xls
rm -rf $xlsPath
mkdir $xlsPath
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -format xls -x Ctrl_D0/Ctrl_D0_m6A/peak.xls -o $xlsPath/mESCs_shCont_D0_m6a.xls
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -format xls -x Ctrl_D6/Ctrl_D6_m6A/peak.xls -o $xlsPath/mESCs_shCont_D6_m6a.xls
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -format xls -x SetD2-KD_D0/SetD2-KD_D0_m6A/peak.xls -o $xlsPath/mESCs_shSetD2_D0_m6a.xls
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -format xls -x SetD2-KD_D6/SetD2-KD_D6_m6A/peak.xls -o $xlsPath/mESCs_shSetD2_D6_m6a.xls

### peak fold_enrichment
cd /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/mESCs/xls
mkdir allPeak
mkdir overlappedPeak
mkdir GOSET
# all peak


exomePeakFC.pl -cutoff 0.1 -FE -x mESCs_shCont_D0_m6a.xls mESCs_shCont_D6_m6a.xls mESCs_shSetD2_D0_m6a.xls mESCs_shSetD2_D6_m6a.xls -o ./temp.tmp
awk 'BEGIN{OFS="\t";FS="\t"; print "chr\tchromStart\tchromEnd\tid\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont-D0\tshCont-D6\tshSetD2-D0\tshSetD2-D6"; }
ARGIND==1{split($4,geneIdI,":");split(geneIdI[1],geneIdInfo,".");geneIdi=geneIdInfo[1];geneNamei=geneIdI[2];hashGeneIdName[geneIdi]=geneNamei;}
ARGIND==2{split($4,geneIdH,".");geneIdh=geneIdH[1]; if (geneIdh in hashGeneIdName){geneNameh=hashGeneIdName[geneIdh]}else{geneNameh="-"};$4=$4"\t"geneNameh; print $0}' \
/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.vM9.annotation.all.longest.exon.bed6 ./temp.tmp > ./allPeak/mESCs_all_samples_m6a_peak_FE.txt

exomePeakFC.pl -cutoff 0.1 -FE -log 2 -x mESCs_shCont_D0_m6a.xls mESCs_shCont_D6_m6a.xls mESCs_shSetD2_D0_m6a.xls mESCs_shSetD2_D6_m6a.xls -o ./temp.tmp
awk 'BEGIN{OFS="\t";FS="\t"; print "chr\tchromStart\tchromEnd\tid\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont-D0\tshCont-D6\tshSetD2-D0\tshSetD2-D6"; }
ARGIND==1{split($4,geneIdI,":");split(geneIdI[1],geneIdInfo,".");geneIdi=geneIdInfo[1];geneNamei=geneIdI[2];hashGeneIdName[geneIdi]=geneNamei;}
ARGIND==2{split($4,geneIdH,".");geneIdh=geneIdH[1]; if (geneIdh in hashGeneIdName){geneNameh=hashGeneIdName[geneIdh]}else{geneNameh="-"};$4=$4"\t"geneNameh; print $0}' \
/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.vM9.annotation.all.longest.exon.bed6 ./temp.tmp > ./allPeak/mESCs_all_samples_m6a_peak_FE_log2.txt

exomePeakFC.pl -cutoff 0.1 -x mESCs_shCont_D0_m6a.xls mESCs_shCont_D6_m6a.xls mESCs_shSetD2_D0_m6a.xls mESCs_shSetD2_D6_m6a.xls -o ./temp.tmp
awk 'BEGIN{OFS="\t";FS="\t"; print "chr\tchromStart\tchromEnd\tid\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont-D0\tshCont-D6\tshSetD2-D0\tshSetD2-D6"; }
ARGIND==1{split($4,geneIdI,":");split(geneIdI[1],geneIdInfo,".");geneIdi=geneIdInfo[1];geneNamei=geneIdI[2];hashGeneIdName[geneIdi]=geneNamei;}
ARGIND==2{split($4,geneIdH,".");geneIdh=geneIdH[1]; if (geneIdh in hashGeneIdName){geneNameh=hashGeneIdName[geneIdh]}else{geneNameh="-"};$4=$4"\t"geneNameh; print $0}' \
/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.vM9.annotation.all.longest.exon.bed6 ./temp.tmp > ./allPeak/mESCs_all_samples_m6a_peak_FC.txt

exomePeakFC.pl -cutoff 0.1 -log 2 -x mESCs_shCont_D0_m6a.xls mESCs_shCont_D6_m6a.xls mESCs_shSetD2_D0_m6a.xls mESCs_shSetD2_D6_m6a.xls -o ./temp.tmp
awk 'BEGIN{OFS="\t";FS="\t"; print "chr\tchromStart\tchromEnd\tid\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont-D0\tshCont-D6\tshSetD2-D0\tshSetD2-D6"; }
ARGIND==1{split($4,geneIdI,":");split(geneIdI[1],geneIdInfo,".");geneIdi=geneIdInfo[1];geneNamei=geneIdI[2];hashGeneIdName[geneIdi]=geneNamei;}
ARGIND==2{split($4,geneIdH,".");geneIdh=geneIdH[1]; if (geneIdh in hashGeneIdName){geneNameh=hashGeneIdName[geneIdh]}else{geneNameh="-"};$4=$4"\t"geneNameh; print $0}' \
/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.vM9.annotation.all.longest.exon.bed6 ./temp.tmp > ./allPeak/mESCs_all_samples_m6a_peak_FC_log2.txt

## GO genesets
#exomePeakFC.pl -cutoff 0.1 -log 2 -x mESCs_shCont_D0_m6a.xls mESCs_shSetD2_D0_m6a.xls -o ./temp.tmp
#awk 'BEGIN{OFS="\t";FS="\t"; print "chr\tchromStart\tchromEnd\tid\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont\tshSetD2"; }
#ARGIND==1{split($4,geneIdI,":");split(geneIdI[1],geneIdInfo,".");geneIdi=geneIdInfo[1];geneNamei=geneIdI[2];hashGeneIdName[geneIdi]=geneNamei;}
#ARGIND==2{split($4,geneIdH,".");geneIdh=geneIdH[1]; if (geneIdh in hashGeneIdName){geneNameh=hashGeneIdName[geneIdh]}else{geneNameh="-"};$4=$4"\t"geneNameh; print $0}' \
#/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.vM9.annotation.all.longest.exon.bed6 ./temp.tmp > ./allPeak/mESCs_D0_m6a_peak_FC_log2.txt
#
#exomePeakFC.pl -cutoff 0.1 -log 2 -x mESCs_shCont_D6_m6a.xls mESCs_shSetD2_D6_m6a.xls -o ./temp.tmp
#awk 'BEGIN{OFS="\t";FS="\t"; print "chr\tchromStart\tchromEnd\tid\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont\tshSetD2"; }
#ARGIND==1{split($4,geneIdI,":");split(geneIdI[1],geneIdInfo,".");geneIdi=geneIdInfo[1];geneNamei=geneIdI[2];hashGeneIdName[geneIdi]=geneNamei;}
#ARGIND==2{split($4,geneIdH,".");geneIdh=geneIdH[1]; if (geneIdh in hashGeneIdName){geneNameh=hashGeneIdName[geneIdh]}else{geneNameh="-"};$4=$4"\t"geneNameh; print $0}' \
#/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.vM9.annotation.all.longest.exon.bed6 ./temp.tmp > ./allPeak/mESCs_D6_m6a_peak_FC_log2.txt
#
#exomePeakFC.pl -cutoff 0.1 -log 2 -x mESCs_shCont_D0_m6a.xls mESCs_shSetD2_D6_m6a.xls -o ./temp.tmp
#awk 'BEGIN{OFS="\t";FS="\t"; print "chr\tchromStart\tchromEnd\tid\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont\tshSetD2"; }
#ARGIND==1{split($4,geneIdI,":");split(geneIdI[1],geneIdInfo,".");geneIdi=geneIdInfo[1];geneNamei=geneIdI[2];hashGeneIdName[geneIdi]=geneNamei;}
#ARGIND==2{split($4,geneIdH,".");geneIdh=geneIdH[1]; if (geneIdh in hashGeneIdName){geneNameh=hashGeneIdName[geneIdh]}else{geneNameh="-"};$4=$4"\t"geneNameh; print $0}' \
#/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.vM9.annotation.all.longest.exon.bed6 ./temp.tmp > ./allPeak/mESCs_shSetD2vsshCont_D6vsD0_m6a_peak_FC_log2.txt
#
#awk 'BEGIN{OFS="\t";};{ if(NR>1){ cutoff=log(1/2)/log(2);if($15<cutoff && $14 > 2){ gsub(/\.[0-9]+/,"" ,$4);print $4; } } }' ./allPeak/mESCs_D0_m6a_peak_FC_log2.txt | sort | uniq | awk 'BEGIN{OFS="\t";print "ENSEMBL"; };{ print $1;} ' > ./GOSET/mESCs_D0_shSetD2_hypomethylation_geneID.txt
#awk 'BEGIN{OFS="\t";};{ if(NR>1){ cutoff=log(1/2)/log(2);if($15<cutoff && $14 > 2){ gsub(/\.[0-9]+/,"" ,$4);print $4; } } }' ./allPeak/mESCs_D6_m6a_peak_FC_log2.txt | sort | uniq | awk 'BEGIN{OFS="\t";print "ENSEMBL"; };{ print $1;} ' > ./GOSET/mESCs_D6_shSetD2_hypomethylation_geneID.txt
#awk 'BEGIN{OFS="\t";};{ if(NR>1){ cutoff=log(1/2)/log(2);if($15<cutoff && $14 > 2){ gsub(/\.[0-9]+/,"" ,$4);print $4; } } }' ./allPeak/mESCs_shSetD2vsshCont_D6vsD0_m6a_peak_FC_log2.txt | sort | uniq | awk 'BEGIN{OFS="\t";print "ENSEMBL"; };{ print $1;} ' > ./GOSET/mESCs_shSetD2vsshCont_D6vsD0_shSetD2_hypomethylation_geneID.txt
#
#awk 'BEGIN{OFS="\t";};{ if(NR>1){ cutoff=log(1)/log(2);if($18<cutoff && $16 < log(0.05)){ gsub(/\.[0-9]+/,"" ,$4);print $4; } } }' ./../R/mESCs_shSetD2_D0_diff_peak.xls | sort | uniq | awk 'BEGIN{OFS="\t";print "ENSEMBL"; };{ print $1;} ' > ./GOSET/mESCs_D0_exomePeak_shSetD2_hypomethylation_geneID.txt
#awk 'BEGIN{OFS="\t";};{ if(NR>1){ cutoff=log(1)/log(2);if($18<cutoff && $16 < log(0.05)){ gsub(/\.[0-9]+/,"" ,$4);print $4; } } }' ./../R/mESCs_shSetD2_D6_diff_peak.xls | sort | uniq | awk 'BEGIN{OFS="\t";print "ENSEMBL"; };{ print $1;} ' > ./GOSET/mESCs_D6_exomePeak_shSetD2_hypomethylation_geneID.txt
#
# overlapped peak

exomePeakFC.pl -cutoff 0.1 -overlap 0 1 1 1 -FE -x mESCs_shCont_D0_m6a.xls mESCs_shCont_D6_m6a.xls mESCs_shSetD2_D0_m6a.xls mESCs_shSetD2_D6_m6a.xls -o ./temp.tmp
awk 'BEGIN{OFS="\t";FS="\t"; print "chr\tchromStart\tchromEnd\tid\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont-D0\tshCont-D6\tshSetD2-D0\tshSetD2-D6"; }
ARGIND==1{split($4,geneIdI,":");split(geneIdI[1],geneIdInfo,".");geneIdi=geneIdInfo[1];geneNamei=geneIdI[2];hashGeneIdName[geneIdi]=geneNamei;}
ARGIND==2{split($4,geneIdH,".");geneIdh=geneIdH[1]; if (geneIdh in hashGeneIdName){geneNameh=hashGeneIdName[geneIdh]}else{geneNameh="-"};$4=$4"\t"geneNameh; print $0}' \
/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.vM9.annotation.all.longest.exon.bed6 ./temp.tmp > ./overlappedPeak/mESCs_all_samples_m6a_peak_FE.txt


exomePeakFC.pl -cutoff 0.1 -overlap 0 1 1 1 -x mESCs_shCont_D0_m6a.xls mESCs_shCont_D6_m6a.xls mESCs_shSetD2_D0_m6a.xls mESCs_shSetD2_D6_m6a.xls -o ./temp.tmp
awk 'BEGIN{OFS="\t";FS="\t"; print "chr\tchromStart\tchromEnd\tid\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont-D0\tshCont-D6\tshSetD2-D0\tshSetD2-D6"; }
ARGIND==1{split($4,geneIdI,":");split(geneIdI[1],geneIdInfo,".");geneIdi=geneIdInfo[1];geneNamei=geneIdI[2];hashGeneIdName[geneIdi]=geneNamei;}
ARGIND==2{split($4,geneIdH,".");geneIdh=geneIdH[1]; if (geneIdh in hashGeneIdName){geneNameh=hashGeneIdName[geneIdh]}else{geneNameh="-"};$4=$4"\t"geneNameh; print $0}' \
/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.vM9.annotation.all.longest.exon.bed6 ./temp.tmp > ./overlappedPeak/mESCs_all_samples_m6a_peak_FC.txt

exomePeakFC.pl -cutoff 0.1 -log 2 -overlap 0 1 1 1 -x mESCs_shCont_D0_m6a.xls mESCs_shCont_D6_m6a.xls mESCs_shSetD2_D0_m6a.xls mESCs_shSetD2_D6_m6a.xls -o ./temp.tmp
awk 'BEGIN{OFS="\t";FS="\t"; print "chr\tchromStart\tchromEnd\tid\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont-D0\tshCont-D6\tshSetD2-D0\tshSetD2-D6"; }
ARGIND==1{split($4,geneIdI,":");split(geneIdI[1],geneIdInfo,".");geneIdi=geneIdInfo[1];geneNamei=geneIdI[2];hashGeneIdName[geneIdi]=geneNamei;}
ARGIND==2{split($4,geneIdH,".");geneIdh=geneIdH[1]; if (geneIdh in hashGeneIdName){geneNameh=hashGeneIdName[geneIdh]}else{geneNameh="-"};$4=$4"\t"geneNameh; print $0}' \
/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.vM9.annotation.all.longest.exon.bed6 ./temp.tmp > ./overlappedPeak/mESCs_all_samples_m6a_peak_FC_log2.txt


exomePeakFC.pl -cutoff 0.1 -overlap 0 1 -x mESCs_shCont_D0_m6a.xls mESCs_shSetD2_D0_m6a.xls -o ./temp.tmp
awk 'BEGIN{OFS="\t";FS="\t"; print "chr\tchromStart\tchromEnd\tid\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont\tshSetD2"; }
ARGIND==1{split($4,geneIdI,":");split(geneIdI[1],geneIdInfo,".");geneIdi=geneIdInfo[1];geneNamei=geneIdI[2];hashGeneIdName[geneIdi]=geneNamei;}
ARGIND==2{split($4,geneIdH,".");geneIdh=geneIdH[1]; if (geneIdh in hashGeneIdName){geneNameh=hashGeneIdName[geneIdh]}else{geneNameh="-"};$4=$4"\t"geneNameh; print $0}' \
/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.vM9.annotation.all.longest.exon.bed6 ./temp.tmp > ./overlappedPeak/mESCs_D0_m6a_peak_FC.txt

exomePeakFC.pl -cutoff 0.1 -log 2 -overlap 0 1 -x mESCs_shCont_D0_m6a.xls mESCs_shSetD2_D0_m6a.xls -o ./temp.tmp
awk 'BEGIN{OFS="\t";FS="\t"; print "chr\tchromStart\tchromEnd\tid\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont\tshSetD2"; }
ARGIND==1{split($4,geneIdI,":");split(geneIdI[1],geneIdInfo,".");geneIdi=geneIdInfo[1];geneNamei=geneIdI[2];hashGeneIdName[geneIdi]=geneNamei;}
ARGIND==2{split($4,geneIdH,".");geneIdh=geneIdH[1]; if (geneIdh in hashGeneIdName){geneNameh=hashGeneIdName[geneIdh]}else{geneNameh="-"};$4=$4"\t"geneNameh; print $0}' \
/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.vM9.annotation.all.longest.exon.bed6 ./temp.tmp > ./overlappedPeak/mESCs_D0_m6a_peak_FC_log2.txt

exomePeakFC.pl -cutoff 0.1 -overlap 0 1 -x mESCs_shCont_D6_m6a.xls mESCs_shSetD2_D6_m6a.xls -o ./temp.tmp
awk 'BEGIN{OFS="\t";FS="\t"; print "chr\tchromStart\tchromEnd\tid\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont\tshSetD2"; }
ARGIND==1{split($4,geneIdI,":");split(geneIdI[1],geneIdInfo,".");geneIdi=geneIdInfo[1];geneNamei=geneIdI[2];hashGeneIdName[geneIdi]=geneNamei;}
ARGIND==2{split($4,geneIdH,".");geneIdh=geneIdH[1]; if (geneIdh in hashGeneIdName){geneNameh=hashGeneIdName[geneIdh]}else{geneNameh="-"};$4=$4"\t"geneNameh; print $0}' \
/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.vM9.annotation.all.longest.exon.bed6 ./temp.tmp > ./overlappedPeak/mESCs_D6_m6a_peak_FC.txt

exomePeakFC.pl -cutoff 0.1 -log 2 -overlap 0 1 -x mESCs_shCont_D6_m6a.xls mESCs_shSetD2_D6_m6a.xls -o ./temp.tmp
awk 'BEGIN{OFS="\t";FS="\t"; print "chr\tchromStart\tchromEnd\tid\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont\tshSetD2"; }
ARGIND==1{split($4,geneIdI,":");split(geneIdI[1],geneIdInfo,".");geneIdi=geneIdInfo[1];geneNamei=geneIdI[2];hashGeneIdName[geneIdi]=geneNamei;}
ARGIND==2{split($4,geneIdH,".");geneIdh=geneIdH[1]; if (geneIdh in hashGeneIdName){geneNameh=hashGeneIdName[geneIdh]}else{geneNameh="-"};$4=$4"\t"geneNameh; print $0}' \
/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.vM9.annotation.all.longest.exon.bed6 ./temp.tmp > ./overlappedPeak/mESCs_D6_m6a_peak_FC_log2.txt

exomePeakFC.pl -cutoff 0.1 -overlap 0 1 1 -x mESCs_shCont_D0_m6a.xls mESCs_shSetD2_D0_m6a.xls mESCs_shSetD2_D6_m6a.xls -o ./temp.tmp
awk 'BEGIN{OFS="\t";FS="\t"; print "chr\tchromStart\tchromEnd\tid\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont\tshSetD2"; }
ARGIND==1{split($4,geneIdI,":");split(geneIdI[1],geneIdInfo,".");geneIdi=geneIdInfo[1];geneNamei=geneIdI[2];hashGeneIdName[geneIdi]=geneNamei;}
ARGIND==2{split($4,geneIdH,".");geneIdh=geneIdH[1]; if (geneIdh in hashGeneIdName){geneNameh=hashGeneIdName[geneIdh]}else{geneNameh="-"};$4=$4"\t"geneNameh; print $0}' \
/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.vM9.annotation.all.longest.exon.bed6 ./temp.tmp > ./overlappedPeak/mESCs_D6vsD0_m6a_peak_FC.txt

exomePeakFC.pl -cutoff 0.1 -overlap 0 1 1 -log 2 -x mESCs_shCont_D0_m6a.xls mESCs_shSetD2_D0_m6a.xls mESCs_shSetD2_D6_m6a.xls -o ./temp.tmp
awk 'BEGIN{OFS="\t";FS="\t"; print "chr\tchromStart\tchromEnd\tid\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tshCont\tshSetD2"; }
ARGIND==1{split($4,geneIdI,":");split(geneIdI[1],geneIdInfo,".");geneIdi=geneIdInfo[1];geneNamei=geneIdI[2];hashGeneIdName[geneIdi]=geneNamei;}
ARGIND==2{split($4,geneIdH,".");geneIdh=geneIdH[1]; if (geneIdh in hashGeneIdName){geneNameh=hashGeneIdName[geneIdh]}else{geneNameh="-"};$4=$4"\t"geneNameh; print $0}' \
/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.vM9.annotation.all.longest.exon.bed6 ./temp.tmp > ./overlappedPeak/mESCs_D6vsD0_m6a_peak_FC_log2.txt

rm -f ./temp.tmp
