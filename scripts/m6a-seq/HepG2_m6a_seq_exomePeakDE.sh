#!/bin/sh
mRNAAnnotation=/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.mRNA.longest.exon.bed6
allAnnotation=/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.all.longest.exon.bed6
rPath=/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/R
bed12Path=/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/bed12
bed6Path=/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/bed6
low=1.2
modest=1.5
high=2

cd /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/exomePeakDE/
rm -rf ./*

awk -v cutoff="$low" 'BEGIN{cutoff=log(1/cutoff)/log(2);}{if(FNR>1)if($18<cutoff){for(i=1;i<12;i++) {printf("%s\t",$i);} printf("%s", $12); printf "\n";}}' OFS="\t" \
    $rPath/shSetD2/HepG2_shSetD2_sig_diff_peak.xls > HepG2_shSetD2_hypo_peak.bed12

awk -v cutoff="$low" 'BEGIN{cutoff=log(1/cutoff)/log(2);}{if(FNR>1)if($18<cutoff){for(i=1;i<12;i++) {printf("%s\t",$i);} printf("%s", $12); printf "\n";}}' OFS="\t" \
    $rPath/shM14/HepG2_shM14_sig_diff_peak.xls > HepG2_shM14_hypo_peak.bed12

awk -v cutoff="$low" 'BEGIN{cutoff=log(1/cutoff)/log(2);}{if(FNR>1)if($18<cutoff){for(i=1;i<12;i++) {printf("%s\t",$i);} printf("%s", $12); printf "\n";}}' OFS="\t" \
    $rPath/shM3/HepG2_shM3_sig_diff_peak.xls > HepG2_shM3_hypo_peak.bed12

awk -v cutoff="$low" 'BEGIN{cutoff=log(1/cutoff)/log(2);}{if(FNR>1)if($18<cutoff){for(i=1;i<12;i++) {printf("%s\t",$i);} printf("%s", $12); printf "\n";}}' OFS="\t" \
    $rPath/shWTAP/HepG2_shWTAP_sig_diff_peak.xls > HepG2_shWTAP_hypo_peak.bed12
exomePeakToSummit.pl -p 1 -b HepG2_shSetD2_hypo_peak.bed12 -o HepG2_shSetD2_hypo_peak.bed6
exomePeakToSummit.pl -p 1 -b HepG2_shM14_hypo_peak.bed12 -o HepG2_shM14_hypo_peak.bed6
exomePeakToSummit.pl -p 1 -b HepG2_shM3_hypo_peak.bed12 -o HepG2_shM3_hypo_peak.bed6
exomePeakToSummit.pl -p 1 -b HepG2_shWTAP_hypo_peak.bed12 -o HepG2_shWTAP_hypo_peak.bed6

bedBinDistribution.pl -strand -smooth move -t percentage --input $bed6Path/HepG2_shCont_m6a.bed6 -bed6 $mRNAAnnotation -o HepG2_shCont_m6a.bed6bin
bedBinDistribution.pl -strand -smooth move -t percentage --input $bed6Path/HepG2_shSetD2_m6a.bed6 -bed6 $mRNAAnnotation -o HepG2_shSetD2_m6a.bed6bin
bedBinDistribution.pl -strand -smooth move -t percentage --input $bed6Path/HepG2_shM14_m6a.bed6 -bed6 $mRNAAnnotation -o HepG2_shM14_m6a.bed6bin
bedBinDistribution.pl -strand -smooth move -t percentage --input $bed6Path/HepG2_shM3_m6a.bed6 -bed6 $mRNAAnnotation -o HepG2_shM3_m6a.bed6bin
bedBinDistribution.pl -strand -smooth move -t percentage --input $bed6Path/HepG2_shWTAP_m6a.bed6 -bed6 $mRNAAnnotation -o HepG2_shWTAP_m6a.bed6bin
bedBinDistribution.pl -strand -smooth move -t percentage --input HepG2_shSetD2_hypo_peak.bed6 -bed6 $mRNAAnnotation -o HepG2_shSetD2_hypo_peak.bed6bin
bedBinDistribution.pl -strand -smooth move -t percentage --input HepG2_shM14_hypo_peak.bed6 -bed6 $mRNAAnnotation -o HepG2_shM14_hypo_peak.bed6bin
bedBinDistribution.pl -strand -smooth move -t percentage --input HepG2_shM3_hypo_peak.bed6 -bed6 $mRNAAnnotation -o HepG2_shM3_hypo_peak.bed6bin
bedBinDistribution.pl -strand -smooth move -t percentage --input HepG2_shWTAP_hypo_peak.bed6 -bed6 $mRNAAnnotation -o HepG2_shWTAP_hypo_peak.bed6bin

paste HepG2_shCont_m6a.bed6bin HepG2_shSetD2_m6a.bed6bin HepG2_shM14_m6a.bed6bin HepG2_shM3_m6a.bed6bin \
    HepG2_shWTAP_m6a.bed6bin HepG2_shSetD2_hypo_peak.bed6bin HepG2_shM14_hypo_peak.bed6bin HepG2_shM3_hypo_peak.bed6bin HepG2_shWTAP_hypo_peak.bed6bin |\
    cut -f 1,2,3,6,9,12,15,18,21,24,27 > HepG2_all_m6A_hypo.bin
sed  -i '1i region\tbin\tshCont\tshSetD2\tshM14\tshM3\tshWTAP\tshSetD2-hypo\tshM14-hypo\tshM3-hypo\tshWTAP-hypo' HepG2_all_m6A_hypo.bin

bedtools intersect -a HepG2_shSetD2_hypo_peak.bed12 -b HepG2_shM14_hypo_peak.bed12 -wa -f 0.5 | sort | uniq > shSetD2_intersect_shM14.bed12
bedtools intersect -a HepG2_shSetD2_hypo_peak.bed12 -b HepG2_shM3_hypo_peak.bed12 -wa -f 0.5 | sort | uniq > shSetD2_intersect_shM3.bed12
bedtools intersect -a HepG2_shSetD2_hypo_peak.bed12 -b HepG2_shWTAP_hypo_peak.bed12 -wa -f 0.5 | sort | uniq > shSetD2_intersect_shWTAP.bed12
cat shSetD2_intersect_shM14.bed12 shSetD2_intersect_shM3.bed12 shSetD2_intersect_shWTAP.bed12 | sort | uniq > HepG2_shSetD2_writer_intersect.bed12

exomePeakToSummit.pl -p 1 -b HepG2_shSetD2_hypo_peak.bed12 -o HepG2_shSetD2_hypo_peak.bed6
exomePeakToSummit.pl -p 1 -b HepG2_shSetD2_writer_intersect.bed12 -o HepG2_shSetD2_writer_intersect.bed6

bedBinDistribution.pl -strand -smooth move -t percentage --input HepG2_shSetD2_hypo_peak.bed6 -bed6 $mRNAAnnotation -o HepG2_shSetD2_hypo_peak.bed6bin
bedBinDistribution.pl -strand -smooth move -t percentage --input HepG2_shSetD2_writer_intersect.bed6 -bed6 $mRNAAnnotation -o HepG2_shSetD2_writer_intersect.bed6bin
paste HepG2_shSetD2_hypo_peak.bed6bin HepG2_shSetD2_writer_intersect.bed6bin | cut -f 1,2,3,6 > HepG2_shSetD2_hypo.bin
sed  -i '1i region\tbin\tHypo-all\tHypo-writer' HepG2_shSetD2_hypo.bin
rm -f *.bed6bin

bedtools intersect -a HepG2_shM14_hypo_peak.bed12 -b shSetD2_intersect_shM14.bed12 -wa -f 0.5 | sort | uniq > shM14_intersect_shSetD2.bed12
bedtools intersect -a HepG2_shM14_hypo_peak.bed12 -b shSetD2_intersect_shM14.bed12 -v -f 0.5 | sort | uniq > shM14_Non_intersect_shSetD2.bed12

exomePeakToSummit.pl -p 1 -b shM14_intersect_shSetD2.bed12 -o shM14_intersect_shSetD2.bed6
exomePeakToSummit.pl -p 1 -b shM14_Non_intersect_shSetD2.bed12 -o shM14_Non_intersect_shSetD2.bed6

bedBinDistribution.pl -strand -smooth move -t percentage --input shM14_intersect_shSetD2.bed6 -bed6 $mRNAAnnotation -o shM14_intersect_shSetD2.bed6bin
bedBinDistribution.pl -strand -smooth move -t percentage --input shM14_Non_intersect_shSetD2.bed6 -bed6 $mRNAAnnotation -o shM14_Non_intersect_shSetD2.bed6bin
paste shM14_intersect_shSetD2.bed6bin shM14_Non_intersect_shSetD2.bed6bin | cut -f 1,2,3,6 > HepG2_shM14_hypo.bin
sed  -i '1i region\tbin\tshM14_intersect_shSetD2\tshM14_Non_intersect_shSetD2' HepG2_shM14_hypo.bin
rm -f *.bed6bin

bedtools intersect -a HepG2_shM3_hypo_peak.bed12 -b shSetD2_intersect_shM3.bed12 -wa -f 0.5 | sort | uniq > shM3_intersect_shSetD2.bed12
bedtools intersect -a HepG2_shM3_hypo_peak.bed12 -b shSetD2_intersect_shM3.bed12 -v -f 0.5 | sort | uniq > shM3_Non_intersect_shSetD2.bed12

exomePeakToSummit.pl -p 1 -b shM3_intersect_shSetD2.bed12 -o shM3_intersect_shSetD2.bed6
exomePeakToSummit.pl -p 1 -b shM3_Non_intersect_shSetD2.bed12 -o shM3_Non_intersect_shSetD2.bed6

bedBinDistribution.pl -strand -smooth move -t percentage --input shM3_intersect_shSetD2.bed6 -bed6 $mRNAAnnotation -o shM3_intersect_shSetD2.bed6bin
bedBinDistribution.pl -strand -smooth move -t percentage --input shM3_Non_intersect_shSetD2.bed6 -bed6 $mRNAAnnotation -o shM3_Non_intersect_shSetD2.bed6bin
paste shM3_intersect_shSetD2.bed6bin shM3_Non_intersect_shSetD2.bed6bin | cut -f 1,2,3,6 > HepG2_shM3_hypo.bin
sed  -i '1i region\tbin\tshM3_intersect_shSetD2\tshM3_Non_intersect_shSetD2' HepG2_shM3_hypo.bin
rm -f *.bed6bin

bedtools intersect -a HepG2_shWTAP_hypo_peak.bed12 -b shSetD2_intersect_shWTAP.bed12 -wa -f 0.5 | sort | uniq > shWTAP_intersect_shSetD2.bed12
bedtools intersect -a HepG2_shWTAP_hypo_peak.bed12 -b shSetD2_intersect_shWTAP.bed12 -v -f 0.5 | sort | uniq > shWTAP_Non_intersect_shSetD2.bed12

exomePeakToSummit.pl -p 1 -b shWTAP_intersect_shSetD2.bed12 -o shWTAP_intersect_shSetD2.bed6
exomePeakToSummit.pl -p 1 -b shWTAP_Non_intersect_shSetD2.bed12 -o shWTAP_Non_intersect_shSetD2.bed6

bedBinDistribution.pl -strand -smooth move -t percentage --input shWTAP_intersect_shSetD2.bed6 -bed6 $mRNAAnnotation -o shWTAP_intersect_shSetD2.bed6bin
bedBinDistribution.pl -strand -smooth move -t percentage --input shWTAP_Non_intersect_shSetD2.bed6 -bed6 $mRNAAnnotation -o shWTAP_Non_intersect_shSetD2.bed6bin
paste shWTAP_intersect_shSetD2.bed6bin shWTAP_Non_intersect_shSetD2.bed6bin | cut -f 1,2,3,6 > HepG2_shWTAP_hypo.bin
sed  -i '1i region\tbin\tshWTAP_intersect_shSetD2\tshWTAP_Non_intersect_shSetD2' HepG2_shWTAP_hypo.bin
rm -f *.bed6bin
