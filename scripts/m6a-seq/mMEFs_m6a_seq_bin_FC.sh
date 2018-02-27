#!/bin/sh

cd /data/zhoukr/hhl_setd2_m6a/mouse_MEFs_m6A-seq

mRNAAnnotation=/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.vM9.annotation.mRNA.longest.exon.bed6
allAnnotation=/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.vM9.annotation.all.longest.exon.bed6
bed12Path=/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/mMEFs/bed12
bed6Path=/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/mMEFs/bed6

### bed12
rm -rf $bed12Path
mkdir $bed12Path
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x SetD2-WT_NA/SetD2-WT_NA_m6A/peak.xls -o $bed12Path/mMEFs_shCont_m6a.bed12
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x SetD2-KO_NA/SetD2-KO_NA_m6A/peak.xls -o $bed12Path/mMEFs_shSetD2_m6a.bed12

### bed6
rm -rf $bed6Path
mkdir $bed6Path
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -x SetD2-WT_NA/SetD2-WT_NA_m6A/peak.xls -o $bed6Path/mMEFs_shCont_m6a.bed6
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -x SetD2-KO_NA/SetD2-KO_NA_m6A/peak.xls -o $bed6Path/mMEFs_shSetD2_m6a.bed6


###draw all bed12 peak distribution
cd /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/mMEFs/bed12
bedBinDistribution.pl -strand -smooth move -t count --input mMEFs_shCont_m6a.bed12 -bed6 $mRNAAnnotation -o ./mMEFs_shCont_m6a.bed12bin
bedBinDistribution.pl -strand -smooth move -t count --input mMEFs_shSetD2_m6a.bed12 -bed6 $mRNAAnnotation -o ./mMEFs_shSetD2_m6a.bed12bin
paste mMEFs_shCont_m6a.bed12bin mMEFs_shSetD2_m6a.bed12bin | cut -f 1,2,3,6,9,12 > mMEFs_m6a.bed12bin
sed  -i '1i region\tbin\tshCont\tshSetD2' mMEFs_m6a.bed12bin
rm -f mMEFs_shCont_m6a.bed12bin  mMEFs_shSetD2_m6a.bed12bin

###draw all bed6 peak distribution
cd /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/mMEFs/bed6
bedBinDistribution.pl -strand -smooth move -t count --input mMEFs_shCont_m6a.bed6 -bed6 $mRNAAnnotation -o ./mMEFs_shCont_m6a.bed6bin
bedBinDistribution.pl -strand -smooth move -t count --input mMEFs_shSetD2_m6a.bed6 -bed6 $mRNAAnnotation -o ./mMEFs_shSetD2_m6a.bed6bin
paste mMEFs_shCont_m6a.bed6bin mMEFs_shSetD2_m6a.bed6bin  | cut -f 1,2,3,6,9,12 > mMEFs_m6a.bed6bin
sed  -i '1i region\tbin\tshCont\tshSetD2' mMEFs_m6a.bed6bin
rm -f mMEFs_shCont_m6a.bed6bin mMEFs_shSetD2_m6a.bed6bin

#### SetD2 responsive m6a sites
cd /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/mMEFs/bed12
cat ./../xls/allPeak/mMEFs_shCont_vs_shSetD2_m6a_peak_FC_log2.txt | cut -f 1,2,3,4,6,7,8,9,10,11,12,13,14,15 | awk '{if($14<log(1/2)/log(2)) {for(i=1;i<12;i++) {printf("%s\t",$i);} printf("%s", $12); printf "\n";}}' > ./mMEFs_m6a_shCont_SetD2_responsive.bed12

#### gene type counts and region counts of SetD2 responsive m6a sites
geneDistribution.pl -strand --input ./mMEFs_m6a_shCont_SetD2_responsive.bed12 -bed6 $allAnnotation -o ./mMEFs_m6a_shCont_SetD2_responsive.gene
sed  -i '1i geneType\tpeakNumber' ./mMEFs_m6a_shCont_SetD2_responsive.gene


cd /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/mMEFs/bed6

exomePeakToSummit.pl --name SetD2_responsive_m6a -b ./../bed12/mMEFs_m6a_shCont_SetD2_responsive.bed12 -o ./mMEFs_m6a_shCont_SetD2_responsive.bed6

#### gene type counts and region counts of SetD2 responsive m6a sites
geneDistribution.pl -strand --input ./mMEFs_m6a_shCont_SetD2_responsive.bed6 -bed6 $allAnnotation -o ./mMEFs_m6a_shCont_SetD2_responsive.gene
sed  -i '1i geneType\tpeakNumber' ./mMEFs_m6a_shCont_SetD2_responsive.gene
regionDistribution.pl -strand --input ./mMEFs_m6a_shCont_SetD2_responsive.bed6 -bed6 $allAnnotation -o ./mMEFs_m6a_shCont_SetD2_responsive.region
sed  -i '1i region\tpeakNumber\tenrichment' ./mMEFs_m6a_shCont_SetD2_responsive.region
