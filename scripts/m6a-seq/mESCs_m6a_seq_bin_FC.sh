#!/bin/sh

cd /data/zhoukr/hhl_setd2_m6a/mouse_ESCs_m6A-seq

mRNAAnnotation=/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.vM9.annotation.mRNA.longest.exon.bed6
allAnnotation=/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.vM9.annotation.all.longest.exon.bed6
bed12Path=/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/mESCs/bed12
bed6Path=/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/mESCs/bed6

### bed12
rm -rf $bed12Path
mkdir $bed12Path
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x Ctrl_D0/Ctrl_D0_m6A/peak.xls -o $bed12Path/mESCs_shCont_D0_m6a.bed12
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x Ctrl_D6/Ctrl_D6_m6A/peak.xls -o $bed12Path/mESCs_shCont_D6_m6a.bed12
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x SetD2-KD_D0/SetD2-KD_D0_m6A/peak.xls -o $bed12Path/mESCs_shSetD2_D0_m6a.bed12
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x SetD2-KD_D6/SetD2-KD_D6_m6A/peak.xls -o $bed12Path/mESCs_shSetD2_D6_m6a.bed12

### bed6
rm -rf $bed6Path
mkdir $bed6Path
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -x Ctrl_D0/Ctrl_D0_m6A/peak.xls -o $bed6Path/mESCs_shCont_D0_m6a.bed6
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -x Ctrl_D6/Ctrl_D6_m6A/peak.xls -o $bed6Path/mESCs_shCont_D6_m6a.bed6
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -x SetD2-KD_D0/SetD2-KD_D0_m6A/peak.xls -o $bed6Path/mESCs_shSetD2_D0_m6a.bed6
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -x SetD2-KD_D6/SetD2-KD_D6_m6A/peak.xls -o $bed6Path/mESCs_shSetD2_D6_m6a.bed6

###draw all bed12 peak distribution
cd /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/mESCs/bed12
bedBinDistribution.pl -strand -smooth move -t count --input mESCs_shCont_D0_m6a.bed12 -bed6 $mRNAAnnotation -o ./mESCs_shCont_D0_m6a.bed12bin
bedBinDistribution.pl -strand -smooth move -t count --input mESCs_shCont_D6_m6a.bed12 -bed6 $mRNAAnnotation -o ./mESCs_shCont_D6_m6a.bed12bin
bedBinDistribution.pl -strand -smooth move -t count --input mESCs_shSetD2_D0_m6a.bed12 -bed6 $mRNAAnnotation -o ./mESCs_shSetD2_D0_m6a.bed12bin
bedBinDistribution.pl -strand -smooth move -t count --input mESCs_shSetD2_D6_m6a.bed12 -bed6 $mRNAAnnotation -o ./mESCs_shSetD2_D6_m6a.bed12bin
paste mESCs_shCont_D0_m6a.bed12bin mESCs_shCont_D6_m6a.bed12bin mESCs_shSetD2_D0_m6a.bed12bin mESCs_shSetD2_D6_m6a.bed12bin | cut -f 1,2,3,6,9,12 > mESCs_m6a.bed12bin
sed  -i '1i region\tbin\tshCont-D0\tshCont-D6\tshSetD2-D0\tshSetD2-D6' mESCs_m6a.bed12bin
rm -f mESCs_shCont_D0_m6a.bed12bin mESCs_shCont_D6_m6a.bed12bin mESCs_shSetD2_D0_m6a.bed12bin mESCs_shSetD2_D6_m6a.bed12bin

###draw all bed6 peak distribution
cd /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/mESCs/bed6
bedBinDistribution.pl -strand -smooth move -t count --input mESCs_shCont_D0_m6a.bed6 -bed6 $mRNAAnnotation -o ./mESCs_shCont_D0_m6a.bed6bin
bedBinDistribution.pl -strand -smooth move -t count --input mESCs_shCont_D6_m6a.bed6 -bed6 $mRNAAnnotation -o ./mESCs_shCont_D6_m6a.bed6bin
bedBinDistribution.pl -strand -smooth move -t count --input mESCs_shSetD2_D0_m6a.bed6 -bed6 $mRNAAnnotation -o ./mESCs_shSetD2_D0_m6a.bed6bin
bedBinDistribution.pl -strand -smooth move -t count --input mESCs_shSetD2_D6_m6a.bed6 -bed6 $mRNAAnnotation -o ./mESCs_shSetD2_D6_m6a.bed6bin
paste mESCs_shCont_D0_m6a.bed6bin mESCs_shCont_D6_m6a.bed6bin mESCs_shSetD2_D0_m6a.bed6bin mESCs_shSetD2_D6_m6a.bed6bin | cut -f 1,2,3,6,9,12 > mESCs_m6a.bed6bin
sed  -i '1i region\tbin\tshCont-D0\tshCont-D6\tshSetD2-D0\tshSetD2-D6' mESCs_m6a.bed6bin
rm -f mESCs_shCont_D0_m6a.bed6bin mESCs_shCont_D6_m6a.bed6bin mESCs_shSetD2_D0_m6a.bed6bin mESCs_shSetD2_D6_m6a.bed6bin

#### SetD2 responsive m6a sites
cd /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/mESCs/bed12
cat ./../xls/allPeak/mESCs_all_samples_m6a_peak_FC_log2.txt | cut -f 1,2,3,4,6,7,8,9,10,11,12,13,14,15 | awk '{if(FNR>1){if($14<log(1/2)/log(2)) {for(i=1;i<12;i++) {printf("%s\t",$i);} printf("%s", $12); printf "\n";}}}' > ./mESCs_m6a_shCont_SetD2_responsive-D0.bed12
cat ./../xls/allPeak/mESCs_all_samples_m6a_peak_FC_log2.txt | cut -f 1,2,3,4,6,7,8,9,10,11,12,13,14,15 | awk '{if(FNR>1){if($16<log(1/2)/log(2)) {for(i=1;i<12;i++) {printf("%s\t",$i);} printf("%s", $12); printf "\n";}}}' > ./mESCs_m6a_shCont_SetD2_responsive-D6.bed12

#### gene type counts and region counts of SetD2 responsive m6a sites
geneDistribution.pl -strand --input ./mESCs_m6a_shCont_SetD2_responsive-D0.bed12 -bed6 $allAnnotation -o ./mESCs_m6a_shCont_SetD2_responsive-D0.gene
sed  -i '1i geneType\tpeakNumber' ./mESCs_m6a_shCont_SetD2_responsive-D0.gene
regionDistribution.pl -strand --input ./mESCs_m6a_shCont_SetD2_responsive-D0.bed12 -bed6 $allAnnotation -o ./mESCs_m6a_shCont_SetD2_responsive-D0.region
sed  -i '1i region\tpeakNumber\tenrichment' ./mESCs_m6a_shCont_SetD2_responsive-D0.region

geneDistribution.pl -strand --input ./mESCs_m6a_shCont_SetD2_responsive-D6.bed12 -bed6 $allAnnotation -o ./mESCs_m6a_shCont_SetD2_responsive-D6.gene
sed  -i '1i geneType\tpeakNumber' ./mESCs_m6a_shCont_SetD2_responsive-D6.gene
regionDistribution.pl -strand --input ./mESCs_m6a_shCont_SetD2_responsive-D6.bed12 -bed6 $allAnnotation -o ./mESCs_m6a_shCont_SetD2_responsive-D6.region
sed  -i '1i region\tpeakNumber\tenrichment' ./mESCs_m6a_shCont_SetD2_responsive-D6.region


cd /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/mESCs/bed6

exomePeakToSummit.pl --name SetD2_responsive_m6a -b ./../bed12/mESCs_m6a_shCont_SetD2_responsive-D0.bed12 -o ./mESCs_m6a_shCont_SetD2_responsive-D0.bed6
exomePeakToSummit.pl --name SetD2_responsive_m6a -b ./../bed12/mESCs_m6a_shCont_SetD2_responsive-D6.bed12 -o ./mESCs_m6a_shCont_SetD2_responsive-D6.bed6

#### gene type counts and region counts of SetD2 responsive m6a sites
geneDistribution.pl -strand --input ./mESCs_m6a_shCont_SetD2_responsive-D0.bed6 -bed6 $allAnnotation -o ./mESCs_m6a_shCont_SetD2_responsive-D0.gene
sed  -i '1i geneType\tpeakNumber' ./mESCs_m6a_shCont_SetD2_responsive-D0.gene
regionDistribution.pl -strand --input ./mESCs_m6a_shCont_SetD2_responsive-D0.bed6 -bed6 $allAnnotation -o ./mESCs_m6a_shCont_SetD2_responsive-D0.region
sed  -i '1i region\tpeakNumber\tenrichment' ./mESCs_m6a_shCont_SetD2_responsive-D0.region

geneDistribution.pl -strand --input ./mESCs_m6a_shCont_SetD2_responsive-D6.bed6 -bed6 $allAnnotation -o ./mESCs_m6a_shCont_SetD2_responsive-D6.gene
sed  -i '1i geneType\tpeakNumber' ./mESCs_m6a_shCont_SetD2_responsive-D6.gene
regionDistribution.pl -strand --input ./mESCs_m6a_shCont_SetD2_responsive-D6.bed6 -bed6 $allAnnotation -o ./mESCs_m6a_shCont_SetD2_responsive-D6.region
sed  -i '1i region\tpeakNumber\tenrichment' ./mESCs_m6a_shCont_SetD2_responsive-D6.region
