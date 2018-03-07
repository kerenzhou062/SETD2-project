#!/bin/sh
HepG2_H3K36me3=/data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/bed6/HepG2_macs_shCont_peaks.bed
chipseq=/data/zhoukr/hhl_setd2_m6a/others/encode/HepG2
bed12Path=/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/bed12
path=/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/heterochomatin
statsFile=$path/heterochomatin_m6a_stats.txt
genomeSize=/data/zhoukr/reference/genome/genome_size/nature/human_hg19_chrsize.txt
cd $path

awk '{print $1, $2, $3}' $HepG2_H3K36me3 > HepG2_macs_shCont_peaks.bed
awk '{print $1, $2, $3}' $bed12Path/HepG2_shCont_m6a.bed12 > HepG2_shCont_m6a.bed
awk '{print $1, $2, $3}' $chipseq/narrowPeak_H3K9me3.bed > narrowPeak_H3K9me3.bed

echo "##statistic overlapped number between histone modification(ENCODE) and m6A" > $statsFile
#echo "histone\tm6A_overlapped_peak_number" >> $statsFile
#for i in $chipseq/broadPeak*
#do
#    name=${i##*/}
#    name=${name%%.*}
#    num=`bedtools intersect -a $bed12Path/HepG2_shCont_m6a.bed12 -b $i -wa |sort|uniq| wc -l`
#    echo -e "$name\t$num" >> $statsFile
#done

#for i in $chipseq/narrowPeak*
#do
#    name=${i##*/}
#    name=${name%%.*}
#    num=`bedtools intersect -a $bed12Path/HepG2_shCont_m6a.bed12 -b $i -wa |sort|uniq| wc -l`
#    echo -e "$name\t$num" >> $statsFile
#done
#echo ""
#echo "##statistic overlapped number between H3K36me3 in house and m6A" >> $statsFile
#echo "histone\tm6A_overlapped_peak_number" >> $statsFile
#num=`bedtools intersect -a $bed12Path/HepG2_shCont_m6a.bed12 -b $HepG2_H3K36me3 -wa |sort|uniq| wc -l`
#echo -e "HepG2_H3K36me3_in_house\t$num" >> $statsFile

echo -e "A\tB\toverlapped_peak_number_of_A" >> $statsFile
num=`bedtools intersect -a $bed12Path/HepG2_shCont_m6a.bed12 -b $HepG2_H3K36me3 -wa |sort|uniq| wc -l`
echo -e "HepG2_m6A\tHepG2_H3K36me3_in_house\t$num" >> $statsFile

num=`bedtools intersect -a $bed12Path/HepG2_shCont_m6a.bed12 -b $chipseq/narrowPeak_H3K9me3.bed -wa |sort|uniq| wc -l`
echo -e "HepG2_m6A\tH3K9me3_GENCODE\t$num" >> $statsFile

num=`bedtools intersect -a $bed12Path/HepG2_shCont_m6a.bed12 -b $chipseq/narrowPeak_H3K27me3.bed -wa |sort|uniq| wc -l`
echo -e "HepG2_m6A\tH3K27me3_GENCODE\t$num" >> $statsFile

num=`bedtools intersect -a $HepG2_H3K36me3 -b $chipseq/narrowPeak_H3K9me3.bed -wa |sort|uniq| wc -l`
echo -e "HepG2_H3K36me3_in_house\tH3K9me3_GENCODE\t$num" >> $statsFile

num=`bedtools intersect -a $HepG2_H3K36me3 -b $chipseq/narrowPeak_H3K27me3.bed -wa |sort|uniq| wc -l`
echo -e "HepG2_H3K36me3_in_house\tH3K27me3_GENCODE\t$num" >> $statsFile

/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/Rscripts/HepG2/heterochromatin_permutation_test.r > /dev/null 2>&1

#cd circos/heatmap
#circos -conf heatmap.m6A_H3K9me3.conf > /dev/null 2>&1
#rm -f heatmap.m6A_H3K9me3.svg
#
#circos -conf heatmap.m6A_H3K9me3_partial.conf > /dev/null 2>&1
#rm -f heatmap.m6A_H3K9me3.svg heatmap.m6A_H3K9me3_partial.svg

