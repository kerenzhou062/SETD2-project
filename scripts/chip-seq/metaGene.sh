#!/bin/sh

cd /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq

bedBinDistribution.pl --feature 5PStart,geneBody,3PEnd -span 5 -smooth move -t percentage --input HepG2/bed6//HepG2_macs_shCont_peaks.bed -bed6 /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.mRNA.geneFeature.bed6 -o ./metaGene/HepG2_chip_shCont.bin
bedBinDistribution.pl --feature 5PStart,geneBody,3PEnd -span 5 -smooth move -t percentage --input Hela/bed6/Hela_macs_shCont_peaks.bed -bed6 /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.mRNA.geneFeature.bed6 -o ./metaGene/Hela_chip_shCont.bin
bedBinDistribution.pl --feature 5PStart,geneBody,3PEnd -span 5 -smooth move -t percentage --input mESCs/bed6/mESCs_macs_shCont_peaks.bed -bed6 /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.vM9.mRNA.geneFeature.bed6 -o ./metaGene/mESCs_chip_shCont.bin
bedBinDistribution.pl --feature 5PStart,geneBody,3PEnd -span 5 -smooth move -t percentage --input mMEFs/bed6/mMEFs_macs_shCont_peaks.bed -bed6 /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.vM9.mRNA.geneFeature.bed6 -o ./metaGene/mMEFs_chip_shCont.bin

cd /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/metaGene

paste HepG2_chip_shCont.bin Hela_chip_shCont.bin mESCs_chip_shCont.bin mMEFs_chip_shCont.bin | cut -f 1,2,3,6,9,12 > shCont_chip_all.bin
sed  -i '1i region\tbin\tHepG2\tHela\tmESCs\tmMEFs' shCont_chip_all.bin
rm -f HepG2_chip_shCont.bin Hela_chip_shCont.bin mESCs_chip_shCont.bin mMEFs_chip_shCont.bin
