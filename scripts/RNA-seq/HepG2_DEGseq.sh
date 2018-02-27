#!/bin/sh
gtf=/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.gtf
bamPath=/data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/
rscript=/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/Rscripts/HepG2/HepG2_DEGseq.r
ouput=/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/HepG2/DEGseq/
bed6=/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.all.longest.exon.bed6

rm -rf $ouput
mkdir -p $ouput/DeGenes
cd $ouput

featureCounts -T 10 -a $gtf -g gene_id -F GTF -t gene -M -o temp.txt \
  $bamPath/shCont/HepG2_m6A-seq_shCont_IP_rep1.fastq.sorted.bam \
  $bamPath/shCont/HepG2_m6A-seq_shCont_IP_rep2.fastq.sorted.bam \
  $bamPath/shCont/HepG2_m6A-seq_shCont_IP_rep3.fastq.sorted.bam \
  $bamPath/shCont/HepG2_m6A-seq_shCont_input_rep1.fastq.sorted.bam \
  $bamPath/shCont/HepG2_m6A-seq_shCont_input_rep2.fastq.sorted.bam \
  $bamPath/shCont/HepG2_m6A-seq_shCont_input_rep3.fastq.sorted.bam \
  $bamPath/shSetD2/HepG2_m6A-seq_shSetD2_IP_rep1.fastq.sorted.bam \
  $bamPath/shSetD2/HepG2_m6A-seq_shSetD2_IP_rep2.fastq.sorted.bam \
  $bamPath/shSetD2/HepG2_m6A-seq_shSetD2_IP_rep3.fastq.sorted.bam \
  $bamPath/shSetD2/HepG2_m6A-seq_shSetD2_input_rep1.fastq.sorted.bam \
  $bamPath/shSetD2/HepG2_m6A-seq_shSetD2_input_rep2.fastq.sorted.bam \
  $bamPath/shSetD2/HepG2_m6A-seq_shSetD2_input_rep3.fastq.sorted.bam \
  $bamPath/shM14/HepG2_m6A-seq_shM14_IP_rep1.fastq.sorted.bam \
  $bamPath/shM14/HepG2_m6A-seq_shM14_IP_rep2.fastq.sorted.bam \
  $bamPath/shM14/HepG2_m6A-seq_shM14_IP_rep3.fastq.sorted.bam \
  $bamPath/shM14/HepG2_m6A-seq_shM14_input_rep1.fastq.sorted.bam \
  $bamPath/shM14/HepG2_m6A-seq_shM14_input_rep2.fastq.sorted.bam \
  $bamPath/shM14/HepG2_m6A-seq_shM14_input_rep3.fastq.sorted.bam \
  $bamPath/shM3/HepG2_m6A-seq_shM3_IP_rep1.fastq.sorted.bam \
  $bamPath/shM3/HepG2_m6A-seq_shM3_IP_rep2.fastq.sorted.bam \
  $bamPath/shM3/HepG2_m6A-seq_shM3_IP_rep3.fastq.sorted.bam \
  $bamPath/shM3/HepG2_m6A-seq_shM3_input_rep1.fastq.sorted.bam \
  $bamPath/shM3/HepG2_m6A-seq_shM3_input_rep2.fastq.sorted.bam \
  $bamPath/shM3/HepG2_m6A-seq_shM3_input_rep3.fastq.sorted.bam \
  $bamPath/shWTAP/HepG2_m6A-seq_shWTAP_IP_rep1.fastq.sorted.bam \
  $bamPath/shWTAP/HepG2_m6A-seq_shWTAP_IP_rep2.fastq.sorted.bam \
  $bamPath/shWTAP/HepG2_m6A-seq_shWTAP_IP_rep3.fastq.sorted.bam \
  $bamPath/shWTAP/HepG2_m6A-seq_shWTAP_input_rep1.fastq.sorted.bam \
  $bamPath/shWTAP/HepG2_m6A-seq_shWTAP_input_rep2.fastq.sorted.bam \
  $bamPath/shWTAP/HepG2_m6A-seq_shWTAP_input_rep3.fastq.sorted.bam \
   > /dev/null 2>&1

awk 'BEGIN{OFS="\t";FS="\t";}
    {
        if(FNR>1){
            if(FNR==2){
                print $1,$2,$3,$4,$5,$6,
                "shContIPRep1","shContIPRep2","shContIPRep3","shContInputRep1","shContInputRep2","shContInputRep3",
                "shSetD2IPRep1","shSetD2IPRep2","shSetD2IPRep3","shSetD2InputRep1","shSetD2InputRep2","shSetD2InputRep3",
                "shM14IPRep1","shM14IPRep2","shM14IPRep3","shM14InputRep1","shM14InputRep2","shM14InputRep3",
                "shM3IPRep1","shM3IPRep2","shM3IPRep3","shM3InputRep1","shM3InputRep2","shM3InputRep3",
                "shWTAPIPRep1","shWTAPIPRep2","shWTAPIPRep3","shWTAPInputRep1","shWTAPInputRep2","shWTAPInputRep3";
            }else{
                print $0
            }
        }
    }
' temp.txt > HepG2_RNA-seq_featureCount.txt

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
' shCont/output_score.txt > HepG2_RNA-seq_shCont_DEGseq_FC.txt

awk 'BEGIN{OFS="\t";FS="\t";}
    {
        if(FNR>1){
            if($5>0 && $5!="NA"){
                print $1, $5;
            }
        }
    }
' shSetD2/output_score.txt > HepG2_RNA-seq_shSetD2_DEGseq_FC.txt

awk 'BEGIN{OFS="\t";FS="\t";}
    {
        if(FNR>1){
            if($5>0 && $5!="NA"){
                print $1, $5;
            }
        }
    }
' shM14/output_score.txt > HepG2_RNA-seq_shM14_DEGseq_FC.txt

awk 'BEGIN{OFS="\t";FS="\t";}
    {
        if(FNR>1){
            if($5>0 && $5!="NA"){
                print $1, $5;
            }
        }
    }
' shM3/output_score.txt > HepG2_RNA-seq_shM3_DEGseq_FC.txt

awk 'BEGIN{OFS="\t";FS="\t";}
    {
        if(FNR>1){
            if($5>0 && $5!="NA"){
                print $1, $5;
            }
        }
    }
' shWTAP/output_score.txt > HepG2_RNA-seq_shWTAP_DEGseq_FC.txt

cat *_DEGseq_FC.txt | cut -f 1 | sort | uniq > HepG2_RNA-seq_union_gene.txt

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
' HepG2_RNA-seq_union_gene.txt HepG2_RNA-seq_shCont_DEGseq_FC.txt \
  HepG2_RNA-seq_shSetD2_DEGseq_FC.txt HepG2_RNA-seq_shM14_DEGseq_FC.txt \
  HepG2_RNA-seq_shM3_DEGseq_FC.txt HepG2_RNA-seq_shWTAP_DEGseq_FC.txt > HepG2_RNA-seq_DEGseq_union_FC.txt

cat *_DEGseq_FC.txt | cut -f 1 | sort | uniq -c | awk '$1 == 5 { print $2 }' > HepG2_RNA-seq_intersected_gene.txt

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
' HepG2_RNA-seq_intersected_gene.txt HepG2_RNA-seq_shCont_DEGseq_FC.txt \
  HepG2_RNA-seq_shSetD2_DEGseq_FC.txt HepG2_RNA-seq_shM14_DEGseq_FC.txt \
  HepG2_RNA-seq_shM3_DEGseq_FC.txt HepG2_RNA-seq_shWTAP_DEGseq_FC.txt > HepG2_RNA-seq_DEGseq_intersected_FC.txt


##HepG2_RNA-seq_DEGseq_FC.txt
##GeneID,shCont,shSetD2,shM14,shM3,shWTAP


cd DeGenes/

awk 'BEGIN{OFS="\t";FS="\t";}
ARGIND==1{
  split($4,geneIdI,":");
  geneIdi=geneIdI[1];
  geneNamei=geneIdI[2];
  hashGeneIdName[geneIdi]=geneNamei;
}
ARGIND==2{
  if(FNR==1){
    $1 = "\"GeneID\"\t\"GeneName\"";
    print $0;
  }else{
    $1 = $1"\t"hashGeneIdName[$1];
    print $0;
  }
}' $bed6 shSetD2/output_score.txt > HepG2_RNA-seq_shSetD2_DEGseq_DeGenes.txt

awk 'BEGIN{OFS="\t";}
ARGIND==1{
  split($4,geneIdI,":");
  geneIdi=geneIdI[1];
  geneNamei=geneIdI[2];
  hashGeneIdName[geneIdi]=geneNamei;
}
ARGIND==2{
  if(FNR==1){
    $1 = "\"GeneID\"\t\"GeneName\"";
    print $0;
  }else{
    $1 = $1"\t"hashGeneIdName[$1];
    print $0;
  }
}' $bed6 shM14/output_score.txt > HepG2_RNA-seq_shM14_DEGseq_DeGenes.txt

awk 'BEGIN{OFS="\t";}
ARGIND==1{
  split($4,geneIdI,":");
  geneIdi=geneIdI[1];
  geneNamei=geneIdI[2];
  hashGeneIdName[geneIdi]=geneNamei;
}
ARGIND==2{
  if(FNR==1){
    $1 = "\"GeneID\"\t\"GeneName\"";
    print $0;
  }else{
    $1 = $1"\t"hashGeneIdName[$1];
    print $0;
  }
}' $bed6 shM3/output_score.txt > HepG2_RNA-seq_shM3_DEGseq_DeGenes.txt

awk 'BEGIN{OFS="\t";}
ARGIND==1{
  split($4,geneIdI,":");
  geneIdi=geneIdI[1];
  geneNamei=geneIdI[2];
  hashGeneIdName[geneIdi]=geneNamei;
}
ARGIND==2{
  if(FNR==1){
    $1 = "\"GeneID\"\t\"GeneName\"";
    print $0;
  }else{
    $1 = $1"\t"hashGeneIdName[$1];
    print $0;
  }
}' $bed6 shWTAP/output_score.txt > HepG2_RNA-seq_shWTAP_DEGseq_DeGenes.txt


