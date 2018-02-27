#!/bin/sh

histonePeak=/data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/mESCs/bed6/mESCs_macs_shCont_peaks.bed
rm -rf /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/mESCs/xls/m6a-H3K36me3
mkdir /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/mESCs/xls/m6a-H3K36me3
cd /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/mESCs/xls/m6a-H3K36me3
awk '{$15=log($15)/log(2);print $15}' ./../mESCs_shCont_D0_m6a.xls > mESCs_shCont_m6a.txt
awk '{$15=log($15)/log(2);print $15}' ./../mESCs_shSetD2_D0_m6a.xls > mESCs_shSetD2_m6a.txt
cumulativePlot.pl -header false -interval 0.01 --input mESCs_shCont_m6a.txt mESCs_shSetD2_m6a.txt -o mESCs_m6A_all_cumulativePlot.txt
sed -i '1i FE\tshCont\tshSetD2' mESCs_m6A_all_cumulativePlot.txt

awk '{$4=$4":"$15;for(i=1;i<12;i++) {printf("%s\t",$i);} printf("%s", $12); printf "\n";}' ./../mESCs_shCont_D0_m6a.xls > mESCs_shCont_m6a.bed12
awk '{$4=$4":"$15;for(i=1;i<12;i++) {printf("%s\t",$i);} printf("%s", $12); printf "\n";}' ./../mESCs_shSetD2_D0_m6a.xls > mESCs_shSetD2_m6a.bed12

bed2overlap.pl -aOver -a mESCs_shCont_m6a.bed12 -b $histonePeak -o mESCs_shCont_m6a+H3K36me3.bed12
bed2overlap.pl -aOver -a mESCs_shSetD2_m6a.bed12 -b $histonePeak -o mESCs_shSetD2_m6a+H3K36me3.bed12

awk '{split($4,arr,":");print log(arr[2])/log(2);}' mESCs_shCont_m6a+H3K36me3.bed12 > mESCs_shCont_m6a+H3K36me3.txt
awk '{split($4,arr,":");print log(arr[2])/log(2);}' mESCs_shSetD2_m6a+H3K36me3.bed12 > mESCs_shSetD2_m6a+H3K36me3.txt
cumulativePlot.pl -header false -interval 0.01 --input mESCs_shCont_m6a+H3K36me3.txt mESCs_shSetD2_m6a+H3K36me3.txt -o mESCs_m6A_H3K36me3_cumulativePlot.txt
sed -i '1i FE\tshCont\tshSetD2' mESCs_m6A_H3K36me3_cumulativePlot.txt

bed2overlap.pl -aOver -split -a mESCs_shCont_m6a+H3K36me3.bed12 -b mESCs_shSetD2_m6a+H3K36me3.bed12 -o mESCs_intersect_shCont_m6a+H3K36me3.bed12
bed2overlap.pl -aOver -split -a mESCs_shSetD2_m6a+H3K36me3.bed12 -b mESCs_shCont_m6a+H3K36me3.bed12 -o mESCs_intersect_shSetD2_m6a+H3K36me3.bed12
awk '{split($4,arr,":");print log(arr[2])/log(2);}' mESCs_intersect_shCont_m6a+H3K36me3.bed12 > mESCs_intersect_shCont_m6a+H3K36me3.txt
awk '{split($4,arr,":");print log(arr[2])/log(2);}' mESCs_intersect_shSetD2_m6a+H3K36me3.bed12 > mESCs_intersect_shSetD2_m6a+H3K36me3.txt
cumulativePlot.pl -header false -interval 0.01 --input mESCs_intersect_shCont_m6a+H3K36me3.txt mESCs_intersect_shSetD2_m6a+H3K36me3.txt -o mESCs_intersect_m6A_H3K36me3_cumulativePlot.txt
sed -i '1i FE\tshCont\tshSetD2' mESCs_intersect_m6A_H3K36me3_cumulativePlot.txt

bedtools intersect -a mESCs_shCont_m6a+H3K36me3.bed12 -b mESCs_shSetD2_m6a+H3K36me3.bed12 -wao -split | \
  awk '{if($25>0){if($4 in hashA){if($25>hashB[2]){hashA[$4]=$16;hashB[$4]=$25}}else{hashA[$4]=$16;hashB[$4]=$25}}}
  END{print "shCont", "shSETD2"; for (x in hashA){split(x,arrA,":");split(hashA[x],arrB,":");print log(arrA[2])/log(2),log(arrB[2])/log(2)}}' OFS="\t" > mESCs_intersect_m6A_H3K36me3.txt

### shSETD2_D6 vs shSETD2_D0
mkdir shSETD2_D6vsD0

awk '{$4=$4":"$15;for(i=1;i<12;i++) {printf("%s\t",$i);} printf("%s", $12); printf "\n";}' ./../mESCs_shCont_D0_m6a.xls > shSETD2_D6vsD0/mESCs_shCont_D0_m6a.bed12
awk '{$4=$4":"$15;for(i=1;i<12;i++) {printf("%s\t",$i);} printf("%s", $12); printf "\n";}' ./../mESCs_shSetD2_D0_m6a.xls > shSETD2_D6vsD0/mESCs_shSetD2_D0_m6a.bed12
awk '{$4=$4":"$15;for(i=1;i<12;i++) {printf("%s\t",$i);} printf("%s", $12); printf "\n";}' ./../mESCs_shCont_D6_m6a.xls > shSETD2_D6vsD0/mESCs_shCont_D6_m6a.bed12
awk '{$4=$4":"$15;for(i=1;i<12;i++) {printf("%s\t",$i);} printf("%s", $12); printf "\n";}' ./../mESCs_shSetD2_D6_m6a.xls > shSETD2_D6vsD0/mESCs_shSetD2_D6_m6a.bed12

cd shSETD2_D6vsD0
bedtools intersect -a mESCs_shCont_D0_m6a.bed12 -b mESCs_shCont_D6_m6a.bed12 -wao -split | \
  awk '{if($25>0){if($4 in hashA){if($25>hashB[2]){hashA[$4]=$16;hashB[$4]=$25}}else{hashA[$4]=$16;hashB[$4]=$25}}}
  END{print "D0", "D6"; for (x in hashA){split(x,arrA,":");split(hashA[x],arrB,":");print log(arrA[2])/log(2),log(arrB[2])/log(2)}}' OFS="\t" > mESCs_shCont_D6vsD0_m6a.txt

bedtools intersect -a mESCs_shSetD2_D0_m6a.bed12 -b mESCs_shSetD2_D6_m6a.bed12 -wao -split | \
  awk '{if($25>0){if($4 in hashA){if($25>hashB[2]){hashA[$4]=$16;hashB[$4]=$25}}else{hashA[$4]=$16;hashB[$4]=$25}}}
  END{print "D0", "D6"; for (x in hashA){split(x,arrA,":");split(hashA[x],arrB,":");print log(arrA[2])/log(2),log(arrB[2])/log(2)}}' OFS="\t" > mESCs_shSetD2_D6vsD0_m6a.txt

bedtools intersect -a mESCs_shCont_D0_m6a.bed12 -b mESCs_shSetD2_D0_m6a.bed12 -wao -split | \
  awk '{if($25>0){if($4 in hashA){if($25>hashB[2]){hashA[$4]=$16;hashB[$4]=$25}}else{hashA[$4]=$16;hashB[$4]=$25}}}
  END{print "shCont", "shSETD2"; for (x in hashA){split(x,arrA,":");split(hashA[x],arrB,":");print log(arrA[2])/log(2),log(arrB[2])/log(2)}}' OFS="\t" > mESCs_D0vsD0_m6a.txt

bedtools intersect -a mESCs_shCont_D6_m6a.bed12 -b mESCs_shSetD2_D6_m6a.bed12 -wao -split | \
  awk '{if($25>0){if($4 in hashA){if($25>hashB[2]){hashA[$4]=$16;hashB[$4]=$25}}else{hashA[$4]=$16;hashB[$4]=$25}}}
  END{print "shCont", "shSETD2"; for (x in hashA){split(x,arrA,":");split(hashA[x],arrB,":");print log(arrA[2])/log(2),log(arrB[2])/log(2)}}' OFS="\t" > mESCs_D6vsD6_m6a.txt

