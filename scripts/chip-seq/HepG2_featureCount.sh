#!/bin/sh
featureCountExp=/data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/featureCountExp.py
annotation=/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.gtf
ipBamFile=/data/zhoukr/hhl_setd2_m6a/HepG2_ChIP-seq/shCont/IP/HepG2_ChIP-seq_shCont_IP.fastq.sorted.bam
inputBamFile=/data/zhoukr/hhl_setd2_m6a/HepG2_ChIP-seq/shCont/input/HepG2_ChIP-seq_shCont_input.fastq.sorted.bam
ouput=/data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/featureCount/

rm -rf $ouput
mkdir -p $ouput
cd $ouput

featureCounts -T 10 -a $annotation -g gene_id -F GTF -t gene -M \
  -o $ouput/HepG2_ChIP-seq_shCont_featureCount.txt $ipBamFile $inputBamFile > /dev/null 2>&1

$featureCountExp -m fpkm -i $ouput/HepG2_ChIP-seq_shCont_featureCount.txt -n shCont_IP shCont_input -o $ouput/HepG2_ChIP-seq_shCont_fpkm.txt

awk 'BEGIN{OFS="\t";FS="\t";}
    {if(FNR>2){if($8 > 0 && $7 > 0){FC=$7/$8;if(FC>1){print $1, FC;}}}}
' HepG2_ChIP-seq_shCont_fpkm.txt > HepG2_ChIP-seq_shCont_fpkm_FC.txt
