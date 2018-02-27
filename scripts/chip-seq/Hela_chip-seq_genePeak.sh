#!/bin/sh
script=/data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/
histonePeak=/data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/Hela/bed6/Hela_macs_shCont_peaks.bed
geneGtf=/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.gtf
geneBed=/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.transcript.bed6
output=/data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/Hela/genePeak

rm -rf $output
mkdir $output
cd $output

annotatePeaks.pl $histonePeak none -gtf $geneGtf > annotatePeaksResult.tmp 2> /dev/null

$script/parseHomerAnnotate.pl $geneBed annotatePeaksResult.tmp Hela_shCont_gene_histone.txt
