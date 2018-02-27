#!/bin/sh

script=/data/zhoukr/hhl_setd2_m6a/analysis/joint-analysis/jointAnalysisFpkm.pl
HepGPath=/data/zhoukr/hhl_setd2_m6a/analysis/joint-analysis/HepG2/fpkm
HelaPath=/data/zhoukr/hhl_setd2_m6a/analysis/joint-analysis/Hela/fpkm
rm -rf $HepGPath
mkdir -p $HepGPath
rm -rf $HelaPath
mkdir -p $HelaPath

##HepG2 foldEnrichment
rsemFile=/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/HepG2/RSEM/HepG2_RNA-seq_RSEM_shCont_input.txt
histoneFile=/data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/genePeak/HepG2_shCont_gene_histone.txt
m6aFile=/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/xls/genePeak/HepG2_shCont_gene_m6a.txt
output=$HepGPath/HepG2_join_exp-histone-m6A_foldEnrichment

$script $rsemFile mean $histoneFile $m6aFile $output

##Hela foldEnrichment
rsemFile=/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/Hela/RSEM/Hela_RNA-seq_RSEM_shCont_input.txt
histoneFile=/data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/Hela/genePeak/Hela_shCont_gene_histone.txt
m6aFile=/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/Hela/xls/genePeak/Hela_shCont_gene_m6a.txt
output=$HelaPath/Hela_join_exp-histone-m6A_foldEnrichment

$script $rsemFile rep3 $histoneFile $m6aFile $output

##HepG2 fpkm
rsemFile=/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/HepG2/RSEM/HepG2_RNA-seq_RSEM_shCont_input.txt
histoneFile=/data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/featureCount/HepG2_ChIP-seq_shCont_fpkm_FC.txt
m6aFile=/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/HepG2/m6A_Level/HepG2_m6a_RSEM_shCont_FC.txt
output=$HepGPath/HepG2_join_exp-histone-m6A_fpkm

$script $rsemFile mean $histoneFile $m6aFile $output

##Hela fpkm
rsemFile=/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/Hela/RSEM/Hela_RNA-seq_RSEM_shCont_input.txt
histoneFile=/data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/Hela/featureCount/Hela_ChIP-seq_shCont_fpkm_FC.txt
m6aFile=/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/Hela/m6A_Level/Hela_m6a_RSEM_shCont_FC.txt
output=$HelaPath/Hela_join_exp-histone-m6A_fpkm

$script $rsemFile rep3 $histoneFile $m6aFile $output

