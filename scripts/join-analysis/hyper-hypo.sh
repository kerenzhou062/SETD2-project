#!/bin/sh
##HepG2
cd /data/zhoukr/hhl_setd2_m6a/analysis/joint-analysis/HepG2
genomeSize=/data/zhoukr/reference/genome/genome_size/nature/human_hg19_chrsize.txt

chipPath=/data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/xls/
macs2PeaksFC.pl -cutoff 0.1 -nonOverlapA -nonOverlapB -x $chipPath/HepG2_shCont_chip.xls $chipPath/HepG2_shSetD2_chip.xls -o temp.txt
awk 'BEGIN{OFS="\t";FS="\t";cutoff=0.5;}
{
  if($7<cutoff){
    print $1,$2,$3,$4,$5,$6;
  }
}
' temp.txt > chip_hypo.bed

awk 'BEGIN{OFS="\t";FS="\t";cutoff=0.5;}
{
  if($7>1/cutoff){
    print $1,$2,$3,$4,$5,$6;
  }
}
' temp.txt > chip_hyper.bed

m6aPath=/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/xls/
exome2PeaksFC.pl -cutoff 0.1 -nonOverlapA -nonOverlapB -x $m6aPath/HepG2_shCont_m6a.xls $m6aPath/HepG2_shSetD2_m6a.xls -o temp.txt
awk 'BEGIN{OFS="\t";FS="\t";cutoff=0.5;}
{
  if($13<cutoff){
    print $1,$2,$3,$4"-"$5,$6;
  }
}
' temp.txt > m6a_hypo.bed

awk 'BEGIN{OFS="\t";FS="\t";cutoff=0.5;}
{
  if($13>1/cutoff){
    print $1,$2,$3,$4"-"$5,$6;
  }
}
' temp.txt > m6a_hyper.bed

rm -f temp.txt

chipHypovsM6aHypo=`bedtools fisher -a chip_hypo.bed -b m6a_hypo.bed -g $genomeSize | awk '{if(FNR==3){split($0,arr,": ");print arr[2];}}'`
chipHypovsM6aHyper=`bedtools fisher -a chip_hypo.bed -b m6a_hyper.bed -g $genomeSize | awk '{if(FNR==3){split($0,arr,": ");print arr[2];}}'`
chipHypervsM6aHypo=`bedtools fisher -a chip_hyper.bed -b m6a_hypo.bed -g $genomeSize | awk '{if(FNR==3){split($0,arr,": ");print arr[2];}}'`
chipHypervsM6aHyper=`bedtools fisher -a chip_hyper.bed -b m6a_hyper.bed -g $genomeSize | awk '{if(FNR==3){split($0,arr,": ");print arr[2];}}'`

echo -e "Type\tH3K36me3-hyper\tH3K36me3-hypo" > hyper-hypo_stats.txt
echo -e "m6A-hyper\t$chipHypervsM6aHyper\t$chipHypovsM6aHyper" >> hyper-hypo_stats.txt
echo -e "m6A-hypo\t$chipHypervsM6aHypo\t$chipHypovsM6aHypo" >> hyper-hypo_stats.txt





cd /data/zhoukr/hhl_setd2_m6a/analysis/joint-analysis/Hela
genomeSize=/data/zhoukr/reference/genome/genome_size/nature/human_hg19_chrsize.txt

chipPath=/data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/Hela/xls/
macs2PeaksFC.pl -cutoff 0.1 -nonOverlapA -nonOverlapB -x $chipPath/Hela_shCont_chip.xls $chipPath/Hela_shSetD2_chip.xls -o temp.txt
awk 'BEGIN{OFS="\t";FS="\t";cutoff=0.5;}
{
  if($7<cutoff){
    print $1,$2,$3,$4,$5,$6;
  }
}
' temp.txt > chip_hypo.bed

awk 'BEGIN{OFS="\t";FS="\t";cutoff=0.5;}
{
  if($7>1/cutoff){
    print $1,$2,$3,$4,$5,$6;
  }
}
' temp.txt > chip_hyper.bed

m6aPath=/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/Hela/xls/
exome2PeaksFC.pl -cutoff 0.1 -nonOverlapA -nonOverlapB -x $m6aPath/Hela_shCont_m6a.xls $m6aPath/Hela_shSetD2_m6a.xls -o temp.txt
awk 'BEGIN{OFS="\t";FS="\t";cutoff=0.5;}
{
  if($13<cutoff){
    print $1,$2,$3,$4"-"$5,$6;
  }
}
' temp.txt > m6a_hypo.bed

awk 'BEGIN{OFS="\t";FS="\t";cutoff=0.5;}
{
  if($13>1/cutoff){
    print $1,$2,$3,$4"-"$5,$6;
  }
}
' temp.txt > m6a_hyper.bed

rm -f temp.txt

chipHypovsM6aHypo=`bedtools fisher -a chip_hypo.bed -b m6a_hypo.bed -g $genomeSize | awk '{if(FNR==3){split($0,arr,": ");print arr[2];}}'`
chipHypovsM6aHyper=`bedtools fisher -a chip_hypo.bed -b m6a_hyper.bed -g $genomeSize | awk '{if(FNR==3){split($0,arr,": ");print arr[2];}}'`
chipHypervsM6aHypo=`bedtools fisher -a chip_hyper.bed -b m6a_hypo.bed -g $genomeSize | awk '{if(FNR==3){split($0,arr,": ");print arr[2];}}'`
chipHypervsM6aHyper=`bedtools fisher -a chip_hyper.bed -b m6a_hyper.bed -g $genomeSize | awk '{if(FNR==3){split($0,arr,": ");print arr[2];}}'`

echo -e "Type\tH3K36me3-hyper\tH3K36me3-hypo" > hyper-hypo_stats.txt
echo -e "m6A-hyper\t$chipHypervsM6aHyper\t$chipHypovsM6aHyper" >> hyper-hypo_stats.txt
echo -e "m6A-hypo\t$chipHypervsM6aHypo\t$chipHypovsM6aHypo" >> hyper-hypo_stats.txt

