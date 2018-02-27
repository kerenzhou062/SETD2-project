#!/bin/sh
script=/data/zhoukr/hhl_setd2_m6a/analysis/joint-analysis/jointAnalysisDEGseq.pl
HepGPath=/data/zhoukr/hhl_setd2_m6a/analysis/joint-analysis/HepG2/DEGseq
HelaPath=/data/zhoukr/hhl_setd2_m6a/analysis/joint-analysis/Hela/DEGseq
rm -rf $HepGPath
mkdir -p $HepGPath
rm -rf $HelaPath
mkdir -p $HelaPath

##HepG2 DEGseq
chipPath=/data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/DEGseq/
m6aPath=/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/HepG2/DEGseq/
rsemFile=/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/HepG2/RSEM/HepG2_RNA-seq_RSEM_shCont_input.txt
histoneFile=$chipPath/HepG2_ChIP-seq_DEGseq_intersected_FC.txt
m6aFile=$m6aPath/HepG2_RNA-seq_DEGseq_intersected_FC.txt


$script $rsemFile mean $histoneFile $HepGPath/HepG2_join_exp-histone_DEGseq
$script $rsemFile mean $m6aFile $HepGPath/HepG2_join_exp-m6A_DEGseq

cat $chipPath/HepG2_ChIP-seq_shCont_DEGseq_FC.txt $m6aPath/HepG2_RNA-seq_shCont_DEGseq_FC.txt | \
  cut -f 1 | sort | uniq -c | awk '$1 == 2 { print $2 }' > $HepGPath/HepG2_shCont_intersected_gene.txt

echo -e "GeneID\tH3K36me3\tm6A" > $HepGPath/HepG2_shCont_m6a_H3K36me3_log2.txt

awk 'BEGIN{OFS="\t";FS="\t";}
ARGIND==1{
  hashGene[$1]=1;
}
ARGIND==2{
  hashGeneExp[$1,1]=$2;
}
ARGIND==3{
  hashGeneExp[$1,2]=$2;
}
END{
  for (gene in hashGene) {
    print gene, hashGeneExp[gene,1], hashGeneExp[gene,2];
  }
}
' $HepGPath/HepG2_shCont_intersected_gene.txt $chipPath/HepG2_ChIP-seq_shCont_DEGseq_FC.txt \
  $m6aPath/HepG2_RNA-seq_shCont_DEGseq_FC.txt >> $HepGPath/HepG2_shCont_m6a_H3K36me3_log2.txt

##Hela DEGseq
chipPath=/data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/Hela/DEGseq/
m6aPath=/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/Hela/DEGseq/
rsemFile=/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/Hela/RSEM/Hela_RNA-seq_RSEM_shCont_input.txt
histoneFile=$chipPath/Hela_ChIP-seq_DEGseq_intersected_FC.txt
m6aFile=$m6aPath/Hela_RNA-seq_DEGseq_intersected_FC.txt


$script $rsemFile rep3 $histoneFile $HelaPath/Hela_join_exp-histone_DEGseq
$script $rsemFile rep3 $m6aFile $HelaPath/Hela_join_exp-m6A_DEGseq

cat $chipPath/Hela_ChIP-seq_shCont_DEGseq_FC.txt $m6aPath/Hela_RNA-seq_shCont_DEGseq_FC.txt | \
  cut -f 1 | sort | uniq -c | awk '$1 == 2 { print $2 }' > $HelaPath/Hela_shCont_intersected_gene.txt

echo -e "GeneID\tH3K36me3\tm6A" > $HelaPath/Hela_shCont_m6a_H3K36me3_log2.txt
awk 'BEGIN{OFS="\t";FS="\t";}
ARGIND==1{
  hashGene[$1]=1;
}
ARGIND==2{
  hashGeneExp[$1,1]=$2;
}
ARGIND==3{
  hashGeneExp[$1,2]=$2;
}
END{
  for (gene in hashGene) {
    print gene, hashGeneExp[gene,1], hashGeneExp[gene,2];
  }
}
' $HelaPath/Hela_shCont_intersected_gene.txt $chipPath/Hela_ChIP-seq_shCont_DEGseq_FC.txt \
  $m6aPath/Hela_RNA-seq_shCont_DEGseq_FC.txt >> $HelaPath/Hela_shCont_m6a_H3K36me3_log2.txt

