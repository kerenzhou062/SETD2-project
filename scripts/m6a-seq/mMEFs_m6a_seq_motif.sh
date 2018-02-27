#!/bin/sh

## motif
cd /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/mMEFs/motif
rm -rf ./*
mkdir ./mMEFs_shCont_m6a
mkdir ./mMEFs_shSetD2_m6a
mkdir ./mMEFs_m6a_shCont_SetD2_responsive

findMotifsGenome.pl ./../bed6/mMEFs_shCont_m6a.bed6 mm10 ./mMEFs_shCont_m6a -rna -size 200 -len 5 -noknown -p 5 -cache 1000 > ./mMEFs_shCont_m6a.log 2>&1
findMotifsGenome.pl ./../bed6/mMEFs_shSetD2_m6a.bed6 mm10 ./mMEFs_shSetD2_m6a -rna -size 200 -len 5 -noknown -p 5 -cache 1000 > ./mMEFs_shSetD2_m6a.log 2>&1
findMotifsGenome.pl ./../bed6/mMEFs_m6a_shCont_SetD2_responsive.bed6 mm10 ./mMEFs_m6a_shCont_SetD2_responsive -rna -size 200 -len 5 -noknown -p 5 -cache 1000 > ./mMEFs_m6a_shCont_SetD2_responsive.log 2>&1

