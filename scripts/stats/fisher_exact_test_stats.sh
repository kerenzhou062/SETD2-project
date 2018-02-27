#!/bin/sh

targetPath=/data/zhoukr/hhl_setd2_m6a/analysis/stats
genomeSize=/data/zhoukr/reference/genome/genome_size/nature/human_hg19_chrsize.txt
chipPath=/data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/bed6
m6aPath=/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/bed12
echo -e "# HepG2, a=H3K36me3 peaks, b=m6A peaks"
histoneBed=$chipPath/HepG2_macs_shCont_peaks.bed
m6ABed=$m6aPath/HepG2_shCont_m6a.bed12

fisherResult=`bedtools fisher -a $histoneBed -b $m6ABed -g $genomeSize`

echo -e "$fisherResult\n"

echo -e "# HepG2, a=H3K36me3 peaks, b=SetD2-responsive m6A peaks"
histoneBed=$chipPath/HepG2_macs_shCont_peaks.bed
m6ABed=$m6aPath/HepG2_m6a_shCont_SetD2_responsive.bed12

fisherResult=`bedtools fisher -a $histoneBed -b $m6ABed -g $genomeSize`

echo -e "$fisherResult\n"

echo -e "# HepG2, a=SetD2-responsive H3K36me3 peaks, b=SetD2-responsive m6A peaks"
histoneBed=$chipPath/HepG2_shCont_SetD2_responsive.bed
m6ABed=$m6aPath/HepG2_m6a_shCont_SetD2_responsive.bed12

fisherResult=`bedtools fisher -a $histoneBed -b $m6ABed -g $genomeSize`

echo -e "$fisherResult\n"

#!/bin/sh

targetPath=/data/zhoukr/hhl_setd2_m6a/analysis/stats
genomeSize=/data/zhoukr/reference/genome/genome_size/nature/human_hg19_chrsize.txt
chipPath=/data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/Hela/bed6
m6aPath=/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/Hela/bed12
echo -e "# Hela, a=H3K36me3 peaks, b=m6A peaks"
histoneBed=$chipPath/Hela_macs_shCont_peaks.bed
m6ABed=$m6aPath/Hela_shCont_m6a.bed12

fisherResult=`bedtools fisher -a $histoneBed -b $m6ABed -g $genomeSize`

echo -e "$fisherResult\n"

echo -e "# Hela, a=H3K36me3 peaks, b=SetD2-responsive m6A peaks"
histoneBed=$chipPath/Hela_macs_shCont_peaks.bed
m6ABed=$m6aPath/Hela_m6a_shCont_SetD2_responsive.bed12

fisherResult=`bedtools fisher -a $histoneBed -b $m6ABed -g $genomeSize`

echo -e "$fisherResult\n"

echo -e "# Hela, a=SetD2-responsive H3K36me3 peaks, b=SetD2-responsive m6A peaks"
histoneBed=$chipPath/Hela_shCont_SetD2_responsive.bed
m6ABed=$m6aPath/Hela_m6a_shCont_SetD2_responsive.bed12

fisherResult=`bedtools fisher -a $histoneBed -b $m6ABed -g $genomeSize`

echo -e "$fisherResult\n"

