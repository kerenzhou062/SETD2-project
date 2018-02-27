#!/bin/sh

## motif
cd /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/motif
rm -rf ./*
mkdir ./HepG2_shCont_chip


findMotifsGenome.pl ./../bed6/HepG2_macs_shCont_peaks.bed hg19 ./HepG2_shCont_chip -size given -len 5,6,7 -noknown -p 5 -cache 1000 > ./HepG2_shCont_chip.log 2>&1

