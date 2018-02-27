#!/bin/sh

histonePeak=/data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/Hela/bed6/Hela_macs_shCont_peaks.bed
cd /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/Hela/xls/allPeak
tatgetPath="./target"
if [ ! -d "$tatgetPath" ]
then
  mkdir "$tatgetPath"
fi

cd $tatgetPath

awk '{if(NR>1){cutoff=log(0.5)/log(2); if ($15 < cutoff) {print $14;}}}' ./../Hela_m6a_peak_FC_log2.txt > Hela_shM14_FC_log2.txt
awk '{if(NR>1){cutoff=log(0.5)/log(2); if ($16 < cutoff) {print $14;}}}' ./../Hela_m6a_peak_FC_log2.txt > Hela_shM3_FC_log2.txt
awk '{if(NR>1){cutoff=log(0.5)/log(2); if ($17 < cutoff) {print $14;}}}' ./../Hela_m6a_peak_FC_log2.txt > Hela_shWTAP_FC_log2.txt
awk '{if(NR>1){cutoff=log(0.5)/log(2); if (($15 < cutoff)&&($16 < cutoff)&&($17 < cutoff)) {print $14;}}}' ./../Hela_m6a_peak_FC_log2.txt > Hela_shareTargets_FC_log2.txt
awk '{if(NR>1){cutoff=log(0.5)/log(2); if (($15 >= cutoff)&&($16 >= cutoff)&&($17 >= cutoff)) {print $14;}}}' ./../Hela_m6a_peak_FC_log2.txt > Hela_nonTargets_FC_log2.txt

cumulativePlot.pl -header false -interval 0.01 --input Hela_shM14_FC_log2.txt Hela_shM3_FC_log2.txt Hela_shWTAP_FC_log2.txt Hela_shareTargets_FC_log2.txt Hela_nonTargets_FC_log2.txt -o ./../Hela_cumulativePlot.txt
sed -i '1i FC\tshM14\tshM3\tshWTAP\tshare\tnonTarget' ./../Hela_cumulativePlot.txt

echo -e "Type\tFC" > Hela_m6A_target_violin_Plot.txt
awk '{print "METTL14", $1}' OFS="\t" Hela_shM14_FC_log2.txt >> Hela_m6A_target_violin_Plot.txt
awk '{print "METTL3", $1}' OFS="\t" Hela_shM3_FC_log2.txt >> Hela_m6A_target_violin_Plot.txt
awk '{print "WTAP", $1}' OFS="\t" Hela_shM3_FC_log2.txt >> Hela_m6A_target_violin_Plot.txt
awk '{print "Share", $1}' OFS="\t" Hela_shareTargets_FC_log2.txt >> Hela_m6A_target_violin_Plot.txt
awk '{print "Non-Target", $1}' OFS="\t" Hela_nonTargets_FC_log2.txt >> Hela_m6A_target_violin_Plot.txt

cd /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/Hela/xls/overlappedPeak
if [ ! -d "$tatgetPath" ]
then
  mkdir "$tatgetPath"
fi

cd $tatgetPath

awk '{if(NR>1){cutoff=log(0.5)/log(2); if ($15 < cutoff) {print $14;}}}' ./../Hela_m6a_peak_FC_log2.txt > Hela_shM14_FC_log2.txt
awk '{if(NR>1){cutoff=log(0.5)/log(2); if ($16 < cutoff) {print $14;}}}' ./../Hela_m6a_peak_FC_log2.txt > Hela_shM3_FC_log2.txt
awk '{if(NR>1){cutoff=log(0.5)/log(2); if ($17 < cutoff) {print $14;}}}' ./../Hela_m6a_peak_FC_log2.txt > Hela_shWTAP_FC_log2.txt
awk '{if(NR>1){cutoff=log(0.5)/log(2); if (($15 < cutoff)&&($16 < cutoff)&&($17 < cutoff)) {print $14;}}}' ./../Hela_m6a_peak_FC_log2.txt > Hela_shareTargets_FC_log2.txt
awk '{if(NR>1){cutoff=log(0.5)/log(2); if (($15 >= cutoff)&&($16 >= cutoff)&&($17 >= cutoff)) {print $14;}}}' ./../Hela_m6a_peak_FC_log2.txt > Hela_nonTargets_FC_log2.txt

cumulativePlot.pl -header false -interval 0.01 --input Hela_shM14_FC_log2.txt Hela_shM3_FC_log2.txt Hela_shWTAP_FC_log2.txt Hela_shareTargets_FC_log2.txt Hela_nonTargets_FC_log2.txt -o ./../Hela_cumulativePlot.txt
sed -i '1i FC\tshM14\tshM3\tshWTAP\tshare\tnonTarget' ./../Hela_cumulativePlot.txt

echo -e "Type\tFC" > Hela_m6A_target_violin_Plot.txt
awk '{print "METTL14", $1}' OFS="\t" Hela_shM14_FC_log2.txt >> Hela_m6A_target_violin_Plot.txt
awk '{print "METTL3", $1}' OFS="\t" Hela_shM3_FC_log2.txt >> Hela_m6A_target_violin_Plot.txt
awk '{print "WTAP", $1}' OFS="\t" Hela_shM3_FC_log2.txt >> Hela_m6A_target_violin_Plot.txt
awk '{print "Share", $1}' OFS="\t" Hela_shareTargets_FC_log2.txt >> Hela_m6A_target_violin_Plot.txt
awk '{print "Non-Target", $1}' OFS="\t" Hela_nonTargets_FC_log2.txt >> Hela_m6A_target_violin_Plot.txt

rm -rf /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/Hela/xls/m6a-H3K36me3
mkdir /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/Hela/xls/m6a-H3K36me3
cd /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/Hela/xls/m6a-H3K36me3
awk '{$15=log($15)/log(2);print $15}' ./../Hela_shCont_m6a.xls > Hela_shCont_m6a.txt
awk '{$15=log($15)/log(2);print $15}' ./../Hela_shSetD2_m6a.xls > Hela_shSetD2_m6a.txt
cumulativePlot.pl -header false -interval 0.01 --input Hela_shCont_m6a.txt Hela_shSetD2_m6a.txt -o Hela_m6A_all_cumulativePlot.txt
sed -i '1i FE\tshCont\tshSetD2' Hela_m6A_all_cumulativePlot.txt

awk '{$4=$4":"$15;for(i=1;i<12;i++) {printf("%s\t",$i);} printf("%s", $12); printf "\n";}' ./../Hela_shCont_m6a.xls > Hela_shCont_m6a.bed12
awk '{$4=$4":"$15;for(i=1;i<12;i++) {printf("%s\t",$i);} printf("%s", $12); printf "\n";}' ./../Hela_shSetD2_m6a.xls > Hela_shSetD2_m6a.bed12

bed2overlap.pl -aOver -split -a Hela_shCont_m6a.bed12 -b Hela_shSetD2_m6a.bed12 -o Hela_intersect_shCont_m6a+shSetD2.bed12
bed2overlap.pl -aOver -split -a Hela_shSetD2_m6a.bed12 -b Hela_shCont_m6a.bed12 -o Hela_intersect_shSetD2_m6a+shCont.bed12

awk '{split($4,arr,":");print log(arr[2])/log(2);}' Hela_intersect_shCont_m6a+shSetD2.bed12 > Hela_intersect_shCont_m6a+shSetD2.txt
awk '{split($4,arr,":");print log(arr[2])/log(2);}' Hela_intersect_shSetD2_m6a+shCont.bed12 > Hela_intersect_shSetD2_m6a+shCont.txt
cumulativePlot.pl -header false -interval 0.01 --input Hela_intersect_shCont_m6a+shSetD2.txt Hela_intersect_shSetD2_m6a+shCont.txt -o Hela_intersect_m6A_cumulativePlot.txt
sed -i '1i FE\tshCont\tshSetD2' Hela_intersect_m6A_cumulativePlot.txt

bedtools intersect -a Hela_intersect_shCont_m6a+shSetD2.bed12 -b Hela_intersect_shSetD2_m6a+shCont.bed12 -wao -split | \
  awk '{if($25>0){if($4 in hashA){if($25>hashB[2]){hashA[$4]=$16;hashB[$4]=$25}}else{hashA[$4]=$16;hashB[$4]=$25}}}
  END{print "shCont", "shSETD2"; for (x in hashA){split(x,arrA,":");split(hashA[x],arrB,":");print log(arrA[2])/log(2),log(arrB[2])/log(2)}}' OFS="\t" > Hela_intersect_m6A.txt

bed2overlap.pl -aOver -a Hela_shCont_m6a.bed12 -b $histonePeak -o Hela_shCont_m6a+H3K36me3.bed12
bed2overlap.pl -aOver -a Hela_shSetD2_m6a.bed12 -b $histonePeak -o Hela_shSetD2_m6a+H3K36me3.bed12

awk '{split($4,arr,":");print log(arr[2])/log(2);}' Hela_shCont_m6a+H3K36me3.bed12 > Hela_shCont_m6a+H3K36me3.txt
awk '{split($4,arr,":");print log(arr[2])/log(2);}' Hela_shSetD2_m6a+H3K36me3.bed12 > Hela_shSetD2_m6a+H3K36me3.txt
cumulativePlot.pl -header false -interval 0.01 --input Hela_shCont_m6a+H3K36me3.txt Hela_shSetD2_m6a+H3K36me3.txt -o Hela_m6A_H3K36me3_cumulativePlot.txt
sed -i '1i FE\tshCont\tshSetD2' Hela_m6A_H3K36me3_cumulativePlot.txt

bed2overlap.pl -aOver -split -a Hela_shCont_m6a+H3K36me3.bed12 -b Hela_shSetD2_m6a+H3K36me3.bed12 -o Hela_intersect_shCont_m6a+H3K36me3.bed12
bed2overlap.pl -aOver -split -a Hela_shSetD2_m6a+H3K36me3.bed12 -b Hela_shCont_m6a+H3K36me3.bed12 -o Hela_intersect_shSetD2_m6a+H3K36me3.bed12
awk '{split($4,arr,":");print log(arr[2])/log(2);}' Hela_intersect_shCont_m6a+H3K36me3.bed12 > Hela_intersect_shCont_m6a+H3K36me3.txt
awk '{split($4,arr,":");print log(arr[2])/log(2);}' Hela_intersect_shSetD2_m6a+H3K36me3.bed12 > Hela_intersect_shSetD2_m6a+H3K36me3.txt
cumulativePlot.pl -header false -interval 0.01 --input Hela_intersect_shCont_m6a+H3K36me3.txt Hela_intersect_shSetD2_m6a+H3K36me3.txt -o Hela_intersect_m6A_H3K36me3_cumulativePlot.txt
sed -i '1i FE\tshCont\tshSetD2' Hela_intersect_m6A_H3K36me3_cumulativePlot.txt

bedtools intersect -a Hela_shCont_m6a+H3K36me3.bed12 -b Hela_shSetD2_m6a+H3K36me3.bed12 -wao -split | \
  awk '{if($25>0){if($4 in hashA){if($25>hashB[2]){hashA[$4]=$16;hashB[$4]=$25}}else{hashA[$4]=$16;hashB[$4]=$25}}}
  END{print "shCont", "shSETD2"; for (x in hashA){split(x,arrA,":");split(hashA[x],arrB,":");print log(arrA[2])/log(2),log(arrB[2])/log(2)}}' OFS="\t" > Hela_intersect_m6A_H3K36me3.txt
