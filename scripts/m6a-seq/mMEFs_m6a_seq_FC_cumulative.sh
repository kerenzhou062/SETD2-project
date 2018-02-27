#!/bin/sh

histonePeak=/data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/mMEFs/bed6/mMEFs_macs_shCont_peaks.bed

rm -rf /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/mMEFs/xls/m6a-H3K36me3
mkdir /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/mMEFs/xls/m6a-H3K36me3
cd /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/mMEFs/xls/m6a-H3K36me3
awk '{$15=log($15)/log(2);print $15}' ./../mMEFs_shCont_m6a.xls > mMEFs_shCont_m6a.txt
awk '{$15=log($15)/log(2);print $15}' ./../mMEFs_shSetD2_m6a.xls > mMEFs_shSetD2_m6a.txt
cumulativePlot.pl -header false -interval 0.01 --input mMEFs_shCont_m6a.txt mMEFs_shSetD2_m6a.txt -o mMEFs_m6A_all_cumulativePlot.txt
sed -i '1i FE\tshCont\tshSetD2' mMEFs_m6A_all_cumulativePlot.txt

awk '{$4=$4":"$15;for(i=1;i<12;i++) {printf("%s\t",$i);} printf("%s", $12); printf "\n";}' ./../mMEFs_shCont_m6a.xls > mMEFs_shCont_m6a.bed12
awk '{$4=$4":"$15;for(i=1;i<12;i++) {printf("%s\t",$i);} printf("%s", $12); printf "\n";}' ./../mMEFs_shSetD2_m6a.xls > mMEFs_shSetD2_m6a.bed12

bed2overlap.pl -aOver -a mMEFs_shCont_m6a.bed12 -b $histonePeak -o mMEFs_shCont_m6a+H3K36me3.bed12
bed2overlap.pl -aOver -a mMEFs_shSetD2_m6a.bed12 -b $histonePeak -o mMEFs_shSetD2_m6a+H3K36me3.bed12

awk '{split($4,arr,":");print log(arr[2])/log(2);}' mMEFs_shCont_m6a+H3K36me3.bed12 > mMEFs_shCont_m6a+H3K36me3.txt
awk '{split($4,arr,":");print log(arr[2])/log(2);}' mMEFs_shSetD2_m6a+H3K36me3.bed12 > mMEFs_shSetD2_m6a+H3K36me3.txt
cumulativePlot.pl -header false -interval 0.01 --input mMEFs_shCont_m6a+H3K36me3.txt mMEFs_shSetD2_m6a+H3K36me3.txt -o mMEFs_m6A_H3K36me3_cumulativePlot.txt
sed -i '1i FE\tshCont\tshSetD2' mMEFs_m6A_H3K36me3_cumulativePlot.txt

bed2overlap.pl -aOver -split -a mMEFs_shCont_m6a+H3K36me3.bed12 -b mMEFs_shSetD2_m6a+H3K36me3.bed12 -o mMEFs_intersect_shCont_m6a+H3K36me3.bed12
bed2overlap.pl -aOver -split -a mMEFs_shSetD2_m6a+H3K36me3.bed12 -b mMEFs_shCont_m6a+H3K36me3.bed12 -o mMEFs_intersect_shSetD2_m6a+H3K36me3.bed12
awk '{split($4,arr,":");print log(arr[2])/log(2);}' mMEFs_intersect_shCont_m6a+H3K36me3.bed12 > mMEFs_intersect_shCont_m6a+H3K36me3.txt
awk '{split($4,arr,":");print log(arr[2])/log(2);}' mMEFs_intersect_shSetD2_m6a+H3K36me3.bed12 > mMEFs_intersect_shSetD2_m6a+H3K36me3.txt
cumulativePlot.pl -header false -interval 0.01 --input mMEFs_intersect_shCont_m6a+H3K36me3.txt mMEFs_intersect_shSetD2_m6a+H3K36me3.txt -o mMEFs_intersect_m6A_H3K36me3_cumulativePlot.txt
sed -i '1i FE\tshCont\tshSetD2' mMEFs_intersect_m6A_H3K36me3_cumulativePlot.txt

bedtools intersect -a mMEFs_shCont_m6a+H3K36me3.bed12 -b mMEFs_shSetD2_m6a+H3K36me3.bed12 -wao -split | \
  awk '{if($25>0){if($4 in hashA){if($25>hashB[2]){hashA[$4]=$16;hashB[$4]=$25}}else{hashA[$4]=$16;hashB[$4]=$25}}}
  END{print "shCont", "shSETD2"; for (x in hashA){split(x,arrA,":");split(hashA[x],arrB,":");print log(arrA[2])/log(2),log(arrB[2])/log(2)}}' OFS="\t" > mMEFs_intersect_m6A_H3K36me3.txt
