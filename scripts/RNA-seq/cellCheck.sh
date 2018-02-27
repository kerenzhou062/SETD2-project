#!/bin/sh

## Hela
perl cellLineCCLECheck.pl Hela/RSEM/shCont_input_rep1.genes.results /data/zhoukr/ccle/CCLE_RNAseq_081117.rpkm.gct Hela_rep1 Hela/cellCheck/Hela_rep1_vs_CCLE.txt
perl cellLineCCLECheck.pl Hela/RSEM/shCont_input_rep2.genes.results /data/zhoukr/ccle/CCLE_RNAseq_081117.rpkm.gct Hela_rep2 Hela/cellCheck/Hela_rep2_vs_CCLE.txt
perl cellLineCCLECheck.pl Hela/RSEM/shCont_input_rep3.genes.results /data/zhoukr/ccle/CCLE_RNAseq_081117.rpkm.gct Hela_rep3 Hela/cellCheck/Hela_rep3_vs_CCLE.txt

## hepG2
perl cellLineCCLECheck.pl HepG2/RSEM/shCont_input_rep1.genes.results /data/zhoukr/ccle/CCLE_RNAseq_081117.rpkm.gct HepG2_rep1 HepG2/cellCheck/HepG2_rep1_vs_CCLE.txt
perl cellLineCCLECheck.pl HepG2/RSEM/shCont_input_rep2.genes.results /data/zhoukr/ccle/CCLE_RNAseq_081117.rpkm.gct HepG2_rep2 HepG2/cellCheck/HepG2_rep2_vs_CCLE.txt
perl cellLineCCLECheck.pl HepG2/RSEM/shCont_input_rep3.genes.results /data/zhoukr/ccle/CCLE_RNAseq_081117.rpkm.gct HepG2_rep3 HepG2/cellCheck/HepG2_rep3_vs_CCLE.txt

