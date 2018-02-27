#!/bin/sh
annotation=/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.gtf
bamPath=/data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/
rscript=/data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/Rscripts/Hela_DEGseq.r
ouput=/data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/Hela/DEGseq/

rm -rf $ouput
mkdir -p $ouput
cd $ouput

featureCounts -T 10 -a $annotation -g gene_id -F GTF -t gene -M -o temp.txt \
  $bamPath/shCont/IP_rep1/Hela_ChIP-seq_shCont_IP_rep1.fastq.sorted.bam \
  $bamPath/shCont/IP_rep2/Hela_ChIP-seq_shCont_IP_rep2.fastq.sorted.bam \
  $bamPath/shCont/input_rep1/Hela_ChIP-seq_shCont_input_rep1.fastq.sorted.bam \
  $bamPath/shCont/input_rep2/Hela_ChIP-seq_shCont_input_rep2.fastq.sorted.bam \
  $bamPath/shSetD2/IP_rep1/Hela_ChIP-seq_shSetD2_IP_rep1.fastq.sorted.bam \
  $bamPath/shSetD2/IP_rep2/Hela_ChIP-seq_shSetD2_IP_rep2.fastq.sorted.bam \
  $bamPath/shSetD2/input_rep1/Hela_ChIP-seq_shSetD2_input_rep1.fastq.sorted.bam \
  $bamPath/shSetD2/input_rep2/Hela_ChIP-seq_shSetD2_input_rep2.fastq.sorted.bam \
  $bamPath/shM14/IP_rep1/Hela_ChIP-seq_shM14_IP_rep1.fastq.sorted.bam \
  $bamPath/shM14/input_rep1/Hela_ChIP-seq_shM14_input_rep1.fastq.sorted.bam \
  $bamPath/shM3/IP_rep1/Hela_ChIP-seq_shM3_IP_rep1.fastq.sorted.bam \
  $bamPath/shM3/input_rep1/Hela_ChIP-seq_shM3_input_rep1.fastq.sorted.bam \
  $bamPath/shWTAP/IP_rep1/Hela_ChIP-seq_shWTAP_IP_rep1.fastq.sorted.bam \
  $bamPath/shWTAP/input_rep1/Hela_ChIP-seq_shWTAP_input_rep1.fastq.sorted.bam \
   > /dev/null 2>&1

awk 'BEGIN{OFS="\t";FS="\t";}
    {
        if(FNR>1){
            if(FNR==2){
                print $1,$2,$3,$4,$5,$6,"shContIPRep1","shContIPRep2","shContInputRep1","shContInputRep2",
                "shSetd2IPRep1","shSetd2IPRep2","shSetd2InputRep1","shSetd2InputRep2",
                "shM14IP","shM14Input","shM3IP","shM3Input","shWTAPIP","shWTAPInput";
            }else{
                print $0
            }
        }
    }
' temp.txt > Hela_ChIP-seq_featureCount.txt

rm -f temp.txt

$rscript > /dev/null 2>&1

awk 'BEGIN{OFS="\t";FS="\t";}
    {
        if(FNR>1){
            if($5>0 && $5!="NA"){
                print $1, $5;
            }
        }
    }
' shCont/output_score.txt > Hela_ChIP-seq_shCont_DEGseq_FC.txt

awk 'BEGIN{OFS="\t";FS="\t";}
    {
        if(FNR>1){
            if($5>0 && $5!="NA"){
                print $1, $5;
            }
        }
    }
' shSetD2/output_score.txt > Hela_ChIP-seq_shSetD2_DEGseq_FC.txt

awk 'BEGIN{OFS="\t";FS="\t";}
    {
        if(FNR>1){
            if($5>0 && $5!="NA"){
                print $1, $5;
            }
        }
    }
' shM14/output_score.txt > Hela_ChIP-seq_shM14_DEGseq_FC.txt

awk 'BEGIN{OFS="\t";FS="\t";}
    {
        if(FNR>1){
            if($5>0 && $5!="NA"){
                print $1, $5;
            }
        }
    }
' shM3/output_score.txt > Hela_ChIP-seq_shM3_DEGseq_FC.txt

awk 'BEGIN{OFS="\t";FS="\t";}
    {
        if(FNR>1){
            if($5>0 && $5!="NA"){
                print $1, $5;
            }
        }
    }
' shWTAP/output_score.txt > Hela_ChIP-seq_shWTAP_DEGseq_FC.txt

cat *_DEGseq_FC.txt | cut -f 1 | sort | uniq > Hela_ChIP-seq_union_gene.txt

awk 'BEGIN{OFS="\t";FS="\t";}
ARGIND==1{
  hashGene[$1]=1;
  hashGeneExp[$1,1]="NA";
  hashGeneExp[$1,2]="NA";
  hashGeneExp[$1,3]="NA";
  hashGeneExp[$1,4]="NA";
  hashGeneExp[$1,5]="NA";
}
ARGIND==2{
  hashGeneExp[$1,1]=$2;
}
ARGIND==3{
  hashGeneExp[$1,2]=$2;
}
ARGIND==4{
  hashGeneExp[$1,3]=$2;
}
ARGIND==5{
  hashGeneExp[$1,4]=$2;
}
ARGIND==6{
  hashGeneExp[$1,5]=$2;
}
END{
  for (gene in hashGene) {
    print gene, hashGeneExp[gene,1], hashGeneExp[gene,2],
      hashGeneExp[gene,3], hashGeneExp[gene,4], hashGeneExp[gene,5];
  }
}
' Hela_ChIP-seq_union_gene.txt Hela_ChIP-seq_shCont_DEGseq_FC.txt \
  Hela_ChIP-seq_shSetD2_DEGseq_FC.txt Hela_ChIP-seq_shM14_DEGseq_FC.txt \
  Hela_ChIP-seq_shM3_DEGseq_FC.txt Hela_ChIP-seq_shWTAP_DEGseq_FC.txt > Hela_ChIP-seq_DEGseq_union_FC.txt

cat *_DEGseq_FC.txt | cut -f 1 | sort | uniq -c | awk '$1 == 5 { print $2 }' > Hela_ChIP-seq_intersected_gene.txt

awk 'BEGIN{OFS="\t";FS="\t";}
ARGIND==1{
  hashGene[$1]=1;
  hashGeneExp[$1,1]="NA";
  hashGeneExp[$1,2]="NA";
  hashGeneExp[$1,3]="NA";
  hashGeneExp[$1,4]="NA";
  hashGeneExp[$1,5]="NA";
}
ARGIND==2{
  hashGeneExp[$1,1]=$2;
}
ARGIND==3{
  hashGeneExp[$1,2]=$2;
}
ARGIND==4{
  hashGeneExp[$1,3]=$2;
}
ARGIND==5{
  hashGeneExp[$1,4]=$2;
}
ARGIND==6{
  hashGeneExp[$1,5]=$2;
}
END{
  for (gene in hashGene) {
    print gene, hashGeneExp[gene,1], hashGeneExp[gene,2],
      hashGeneExp[gene,3], hashGeneExp[gene,4], hashGeneExp[gene,5];
  }
}
' Hela_ChIP-seq_intersected_gene.txt Hela_ChIP-seq_shCont_DEGseq_FC.txt \
  Hela_ChIP-seq_shSetD2_DEGseq_FC.txt Hela_ChIP-seq_shM14_DEGseq_FC.txt \
  Hela_ChIP-seq_shM3_DEGseq_FC.txt Hela_ChIP-seq_shWTAP_DEGseq_FC.txt > Hela_ChIP-seq_DEGseq_intersected_FC.txt


##Hela_ChIP-seq_DEGseq_FC.txt
##GeneID,shCont,shSetD2,shM14,shM3,shWTAP
