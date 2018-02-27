#!/bin/sh
/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/HepG2_m6A_shCont_input_htseq-count-gene.sh > HepG2_m6A_shCont_input_htseq-count-gene.log 2>&1 &
/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/HepG2_m6A_shCont_input_htseq-count-transcript.sh > HepG2_m6A_shCont_input_htseq-count-transcript.log 2>&1 &

/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/HepG2_m6A_shCont_IP_htseq-count-gene.sh > HepG2_m6A_shCont_IP_htseq-count-gene.log 2>&1 &
/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/HepG2_m6A_shCont_IP_htseq-count-transcript.sh > HepG2_m6A_shCont_IP_htseq-count-transcript.log 2>&1 &


/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/Hela_m6A_shCont_input_htseq-count-gene.sh > Hela_m6A_shCont_input_htseq-count-gene.log 2>&1 &
/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/Hela_m6A_shCont_input_htseq-count-transcript.sh > Hela_m6A_shCont_input_htseq-count-transcript.log 2>&1 &

/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/Hela_m6A_shCont_IP_htseq-count-gene.sh > Hela_m6A_shCont_IP_htseq-count-gene.log 2>&1 &
/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/Hela_m6A_shCont_IP_htseq-count-transcript.sh > Hela_m6A_shCont_IP_htseq-count-transcript.log 2>&1 &

