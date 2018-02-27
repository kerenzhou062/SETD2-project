#!/bin/sh
cd /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/htseq-count
find . -type d | sort | grep "gene-count" | awk '{if(FNR>1)print}' | xargs -I {} rm -rf {}
mkdir -p gene-count

nohup /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/htseq-count/bash/HepG2_shCont_input_rep1.sh > /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/htseq-count/bash/HepG2_shCont_input_rep1.log 2>&1 &
nohup /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/htseq-count/bash/HepG2_shCont_input_rep2.sh > /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/htseq-count/bash/HepG2_shCont_input_rep2.log 2>&1 &
nohup /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/htseq-count/bash/HepG2_shCont_IP_rep1.sh > /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/htseq-count/bash/HepG2_shCont_IP_rep1.log 2>&1 &
nohup /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/htseq-count/bash/HepG2_shCont_IP_rep2.sh > /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/htseq-count/bash/HepG2_shCont_IP_rep2.log 2>&1 &
nohup /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/htseq-count/bash/HepG2_shSetD2_input_rep1.sh > /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/htseq-count/bash/HepG2_shSetD2_input_rep1.log 2>&1 &
nohup /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/htseq-count/bash/HepG2_shSetD2_input_rep2.sh > /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/htseq-count/bash/HepG2_shSetD2_input_rep2.log 2>&1 &
nohup /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/htseq-count/bash/HepG2_shSetD2_IP_rep1.sh > /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/htseq-count/bash/HepG2_shSetD2_IP_rep1.log 2>&1 &
nohup /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/htseq-count/bash/HepG2_shSetD2_IP_rep2.sh > /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/htseq-count/bash/HepG2_shSetD2_IP_rep2.log 2>&1 &

cd gene-count
htseqCountGroup.pl -f HepG2_shCont_input_htseq-count_rep1.txt HepG2_shCont_input_htseq-count_rep2.txt \
  HepG2_shCont_IP_htseq-count_rep1.txt HepG2_shCont_IP_htseq-count_rep2.txt \
  -o ./../HepG2_chip-seq_shCont_htseq-count.txt
sed -i '1i geneID\tshContInput-1\tshContInput-2\tshContIP-1\tshContIP-2' ./../HepG2_chip-seq_shCont_htseq-count.txt

htseqCountGroup.pl -f HepG2_shSetD2_input_htseq-count_rep1.txt HepG2_shSetD2_input_htseq-count_rep2.txt \
  HepG2_shSetD2_IP_htseq-count_rep1.txt HepG2_shSetD2_IP_htseq-count_rep2.txt \
  -o ./../HepG2_chip-seq_shSetD2_htseq-count.txt
sed -i '1i geneID\tshSetD2Input-1\tshSetD2Input-2\tshSetD2IP-1\tshSetD2IP-2' ./../HepG2_chip-seq_shSetD2_htseq-count.txt

