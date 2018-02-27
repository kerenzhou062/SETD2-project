#!/bin/sh

## motif
cd /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/Hela/motif
rm -rf ./*
mkdir ./Hela_shCont_m6a
mkdir ./Hela_shSetD2_m6a
mkdir ./Hela_shM14_m6a
mkdir ./Hela_shM3_m6a
mkdir ./Hela_shWTAP_m6a
mkdir ./Hela_m6a_shCont_SetD2_responsive

findMotifsGenome.pl ./../bed12/Hela_shCont_m6a.bed12 hg19 ./Hela_shCont_m6a -rna -size 200 -len 5,6,7 -noknown -p 5 -cache 1000 > ./Hela_shCont_m6a.log 2>&1
findMotifsGenome.pl ./../bed12/Hela_shSetD2_m6a.bed12 hg19 ./Hela_shSetD2_m6a -rna -size 200 -len 5,6,7 -noknown -p 5 -cache 1000 > ./Hela_shSetD2_m6a.log 2>&1
findMotifsGenome.pl ./../bed12/Hela_shM14_m6a.bed12 hg19 ./Hela_shM14_m6a -rna -size 200 -len 5,6,7 -noknown -p 5 -cache 1000 > ./Hela_shM14_m6a.log 2>&1
findMotifsGenome.pl ./../bed12/Hela_shM3_m6a.bed12 hg19 ./Hela_shM3_m6a -rna -size 200 -len 5,6,7 -noknown -p 5 -cache 1000 > ./Hela_shM3_m6a.log 2>&1
findMotifsGenome.pl ./../bed12/Hela_shWTAP_m6a.bed12 hg19 ./Hela_shWTAP_m6a -rna -size 200 -len 5,6,7 -noknown -p 5 -cache 1000 > ./Hela_shWTAP_m6a.log 2>&1
findMotifsGenome.pl ./../bed12/Hela_m6a_shCont_SetD2_responsive.bed12 hg19 ./Hela_m6a_shCont_SetD2_responsive -rna -size 200 -len 5,6,7 -noknown -p 5 -cache 1000 > ./Hela_m6a_shCont_SetD2_responsive.log 2>&1

