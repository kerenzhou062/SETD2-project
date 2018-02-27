#!/bin/sh
mRNAAnnotation=/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.mRNA.longest.exon.bed6
allAnnotation=/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.all.longest.exon.bed6
histonePeak=/data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/bed6/HepG2_macs_shCont_peaks.bed
bed12Path=/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/bed12
bed6Path=/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/bed6
bamPath=/data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq
readsPath=/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/reads
rm -rf $readsPath
mkdir $readsPath
cd $readsPath

bedReadsCount.pl --input $bed12Path/HepG2_m6a_shCont_M14_responsive.bed12 -o HepG2_m6a_shCont_M14_responsive_reads.bed12 -bam $bamPath/shCont/HepG2_m6A-seq_shCont_IP_rep1.fastq.sorted.bam $bamPath/shCont/HepG2_m6A-seq_shCont_IP_rep2.fastq.sorted.bam $bamPath/shCont/HepG2_m6A-seq_shCont_IP_rep3.fastq.sorted.bam $bamPath/shCont/HepG2_m6A-seq_shCont_input_rep1.fastq.sorted.bam $bamPath/shCont/HepG2_m6A-seq_shCont_input_rep2.fastq.sorted.bam $bamPath/shCont/HepG2_m6A-seq_shCont_input_rep3.fastq.sorted.bam
bedReadsCount.pl --input $bed12Path/HepG2_m6a_shSetD2_M14_responsive.bed12 -o HepG2_m6a_shSetD2_M14_responsive_reads.bed12 -bam $bamPath/shSetD2/HepG2_m6A-seq_shSetD2_IP_rep1.fastq.sorted.bam $bamPath/shSetD2/HepG2_m6A-seq_shSetD2_IP_rep2.fastq.sorted.bam $bamPath/shSetD2/HepG2_m6A-seq_shSetD2_IP_rep3.fastq.sorted.bam $bamPath/shSetD2/HepG2_m6A-seq_shSetD2_input_rep1.fastq.sorted.bam $bamPath/shSetD2/HepG2_m6A-seq_shSetD2_input_rep2.fastq.sorted.bam $bamPath/shSetD2/HepG2_m6A-seq_shSetD2_input_rep3.fastq.sorted.bam

bedBinDistribution.pl -strand -smooth move -span 10 -t count -method center -score 'NA' --input HepG2_m6a_shCont_M14_responsive_reads.bed12 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shCont_M14_responsive_reads.bed12bin
bedBinDistribution.pl -strand -smooth move -span 10 -t count -method center -score 'NA' --input HepG2_m6a_shSetD2_M14_responsive_reads.bed12 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shSetD2_M14_responsive_reads.bed12bin
paste HepG2_m6a_shCont_M14_responsive_reads.bed12bin HepG2_m6a_shSetD2_M14_responsive_reads.bed12bin | cut -f 1,2,3,6 > HepG2_m6a_shCont_M14_reads.bed12bin
sed  -i '1i region\tbin\tshCont\tshSetD2' HepG2_m6a_shCont_M14_reads.bed12bin
rm -f HepG2_m6a_shCont_M14_responsive_reads.bed12bin HepG2_m6a_shSetD2_M14_responsive_reads.bed12bin

exomePeakToSummit.pl --name HepG2_m6a_shCont_M14_responsive_reads -b HepG2_m6a_shCont_M14_responsive_reads.bed12 --pvalue 1e10 -o HepG2_m6a_shCont_M14_responsive_reads.bed6
exomePeakToSummit.pl --name HepG2_m6a_shSetD2_M14_responsive_reads -b HepG2_m6a_shSetD2_M14_responsive_reads.bed12 --pvalue 1e10 -o HepG2_m6a_shSetD2_M14_responsive_reads.bed6
bedBinDistribution.pl -strand -smooth move -span 10 -t count -method center -score 'NA' --input HepG2_m6a_shCont_M14_responsive_reads.bed6 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shCont_M14_responsive_reads.bed6bin
bedBinDistribution.pl -strand -smooth move -span 10 -t count -method center -score 'NA' --input HepG2_m6a_shSetD2_M14_responsive_reads.bed6 -bed6 $mRNAAnnotation -o ./HepG2_m6a_shSetD2_M14_responsive_reads.bed6bin
paste HepG2_m6a_shCont_M14_responsive_reads.bed6bin HepG2_m6a_shSetD2_M14_responsive_reads.bed6bin | cut -f 1,2,3,6 > HepG2_m6a_shCont_M14_reads.bed6bin
sed  -i '1i region\tbin\tshCont\tshSetD2' HepG2_m6a_shCont_M14_reads.bed6bin
rm -f HepG2_m6a_shCont_M14_responsive_reads.bed6bin HepG2_m6a_shSetD2_M14_responsive_reads.bed6bin

