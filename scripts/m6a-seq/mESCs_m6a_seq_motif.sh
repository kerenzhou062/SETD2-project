#!/bin/sh

## motif
cd /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/mESCs/motif
rm -rf ./*
mkdir ./mESCs_shCont_D0_m6a
mkdir ./mESCs_shCont_D6_m6a
mkdir ./mESCs_shSetD2_D0_m6a
mkdir ./mESCs_shSetD2_D6_m6a
mkdir ./mESCs_m6a_shCont_SetD2_responsive-D0
mkdir ./mESCs_m6a_shCont_SetD2_responsive-D6

findMotifsGenome.pl ./../bed6/mESCs_shCont_D0_m6a.bed6 mm10 ./mESCs_shCont_D0_m6a -rna -size 200 -len 5,6,7 -noknown -p 5 -cache 1000 > ./mESCs_shCont_D0_m6a.log 2>&1
findMotifsGenome.pl ./../bed6/mESCs_shCont_D6_m6a.bed6 mm10 ./mESCs_shCont_D6_m6a -rna -size 200 -len 5,6,7 -noknown -p 5 -cache 1000 > ./mESCs_shCont_D6_m6a.log 2>&1
findMotifsGenome.pl ./../bed6/mESCs_shSetD2_D0_m6a.bed6 mm10 ./mESCs_shSetD2_D0_m6a -rna -size 200 -len 5,6,7 -noknown -p 5 -cache 1000 > ./mESCs_shSetD2_D0_m6a.log 2>&1
findMotifsGenome.pl ./../bed6/mESCs_shSetD2_D6_m6a.bed6 mm10 ./mESCs_shSetD2_D6_m6a -rna -size 200 -len 5,6,7 -noknown -p 5 -cache 1000 > ./mESCs_shSetD2_D6_m6a.log 2>&1

findMotifsGenome.pl ./../bed6/mESCs_m6a_shCont_SetD2_responsive-D0.bed6 mm10 ./mESCs_m6a_shCont_SetD2_responsive-D0 -rna -size 200 -len 5,6,7 -noknown -p 5 -cache 1000 > ./mESCs_m6a_shCont_SetD2_responsive-D0.log 2>&1
findMotifsGenome.pl ./../bed6/mESCs_m6a_shCont_SetD2_responsive-D6.bed6 mm10 ./mESCs_m6a_shCont_SetD2_responsive-D6 -rna -size 200 -len 5,6,7 -noknown -p 5 -cache 1000 > ./mESCs_m6a_shCont_SetD2_responsive-D6.log 2>&1

