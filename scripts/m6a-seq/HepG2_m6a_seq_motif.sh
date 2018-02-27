#!/bin/sh

## motif
cd /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/motif
rm -rf ./*
mkdir ./HepG2_shCont_m6a
mkdir ./HepG2_shSetD2_m6a
mkdir ./HepG2_shM14_m6a
mkdir ./HepG2_shM3_m6a
mkdir ./HepG2_shWTAP_m6a
mkdir ./HepG2_m6a_shCont_SetD2_responsive

findMotifsGenome.pl ./../bed12/HepG2_shCont_m6a.bed12 hg19 ./HepG2_shCont_m6a -rna -size 200 -len 5,6,7 -noknown -p 5 -cache 1000 > ./HepG2_shCont_m6a.log 2>&1
findMotifsGenome.pl ./../bed12/HepG2_shSetD2_m6a.bed12 hg19 ./HepG2_shSetD2_m6a -rna -size 200 -len 5,6,7 -noknown -p 5 -cache 1000 > ./HepG2_shSetD2_m6a.log 2>&1
findMotifsGenome.pl ./../bed12/HepG2_shM14_m6a.bed12 hg19 ./HepG2_shM14_m6a -rna -size 200 -len 5,6,7 -noknown -p 5 -cache 1000 > ./HepG2_shM14_m6a.log 2>&1
findMotifsGenome.pl ./../bed12/HepG2_shM3_m6a.bed12 hg19 ./HepG2_shM3_m6a -rna -size 200 -len 5,6,7 -noknown -p 5 -cache 1000 > ./HepG2_shM3_m6a.log 2>&1
findMotifsGenome.pl ./../bed12/HepG2_shWTAP_m6a.bed12 hg19 ./HepG2_shWTAP_m6a -rna -size 200 -len 5,6,7 -noknown -p 5 -cache 1000 > ./HepG2_shWTAP_m6a.log 2>&1
findMotifsGenome.pl ./../bed12/HepG2_m6a_shCont_SetD2_responsive.bed12 hg19 ./HepG2_m6a_shCont_SetD2_responsive -rna -size 200 -len 5,6,7 -noknown -p 5 -cache 1000 > ./HepG2_m6a_shCont_SetD2_responsive.log 2>&1

