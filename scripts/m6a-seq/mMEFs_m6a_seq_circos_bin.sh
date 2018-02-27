#!/bin/sh

cd /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/mMEFs/circos
exome2PeaksFC.pl -cutoff 0.1 -nonOverlapA -nonOverlapB -x ../xls/mMEFs_shSetD2_m6a.xls ../xls/mMEFs_shCont_m6a.xls -o temp.xls
awk '{gsub(/chr/,"mm",$1);$13=log(1/$13)/log(2);print $1,$2,$3,$13;}' OFS="\t" temp.xls > mESCs_m6a_shSetD2_D0_FC_circos.txt

rm -f temp.xls
