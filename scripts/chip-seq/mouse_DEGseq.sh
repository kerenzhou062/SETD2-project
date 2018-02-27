#!/bin/sh
annotation=/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.vM9.annotation.gtf
mESCsbamPath=/data/zhoukr/hhl_setd2_m6a/mouse_ESCs_ChIP-seq
mEFsbamPath=/data/zhoukr/hhl_setd2_m6a/mouse_MEFs_ChIP-seq
rscript=/data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/Rscripts/mouse_DEGseq.r
ouput=/data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/mouse/DEGseq/

rm -rf $ouput
mkdir -p $ouput
cd $ouput

featureCounts -T 10 -a $annotation -g gene_id -F GTF -t gene -M -o temp.txt \
  $mESCsbamPath/Ctrl_D0/mESCs_ChIP-seq_Ctrl_D0_ChIP-rep1.sorted.bam \
  $mESCsbamPath/Ctrl_D0/mESCs_ChIP-seq_Ctrl_D0_ChIP-rep2.sorted.bam \
  $mESCsbamPath/Ctrl_D0/mESCs_ChIP-seq_Ctrl_D0_input-rep1.sorted.bam \
  $mESCsbamPath/Ctrl_D0/mESCs_ChIP-seq_Ctrl_D0_input-rep2.sorted.bam \
  $mESCsbamPath/SetD2-KD_D0/mESCs_ChIP-seq_SetD2-KD_D0_ChIP-rep1.sorted.bam \
  $mESCsbamPath/SetD2-KD_D0/mESCs_ChIP-seq_SetD2-KD_D0_ChIP-rep2.sorted.bam \
  $mESCsbamPath/SetD2-KD_D0/mESCs_ChIP-seq_SetD2-KD_D0_input-rep1.sorted.bam \
  $mESCsbamPath/SetD2-KD_D0/mESCs_ChIP-seq_SetD2-KD_D0_input-rep2.sorted.bam \
  $mEFsbamPath/SetD2-WT_NA/mMEFs_ChIP-seq_SetD2-WT_NA_ChIP-rep1.sorted.bam \
  $mEFsbamPath/SetD2-WT_NA/mMEFs_ChIP-seq_SetD2-WT_NA_ChIP-rep2.sorted.bam \
  $mEFsbamPath/SetD2-WT_NA/mMEFs_ChIP-seq_SetD2-WT_NA_input-rep1.sorted.bam \
  $mEFsbamPath/SetD2-WT_NA/mMEFs_ChIP-seq_SetD2-WT_NA_input-rep2.sorted.bam \
  $mEFsbamPath/SetD2-KO_NA/mMEFs_ChIP-seq_SetD2-KO_NA_ChIP-rep1.sorted.bam \
  $mEFsbamPath/SetD2-KO_NA/mMEFs_ChIP-seq_SetD2-KO_NA_ChIP-rep2.sorted.bam \
  $mEFsbamPath/SetD2-KO_NA/mMEFs_ChIP-seq_SetD2-KO_NA_input-rep1.sorted.bam \
  $mEFsbamPath/SetD2-KO_NA/mMEFs_ChIP-seq_SetD2-KO_NA_input-rep2.sorted.bam \
   > /dev/null 2>&1

awk 'BEGIN{OFS="\t";FS="\t";}
    {
        if(FNR>1){
            if(FNR==2){
                print $1,$2,$3,$4,$5,$6,
                "mESCsShContD0IPRep1","mESCsShContD0IPRep2","mESCsShContD0InputRep1","mESCsShContD0InputRep2",
                "mESCsShSetD2D0IPRep1","mESCsShSetD2D0IPRep2","mESCsShSetD2D0InputRep1","mESCsShSetD2D0InputRep2",
                "mEFsShContIPRep1","mEFsShContIPRep2","mEFsShContInputRep1","mEFsShContInputRep2",
                "mEFsShSetD2IPRep1","mEFsShSetD2IPRep2","mEFsShSetD2InputRep1","mEFsShSetD2InputRep2";
            }else{
                print $0
            }
        }
    }
' temp.txt > mouse_ChIP-seq_featureCount.txt

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
' mESCsShContD0/output_score.txt > mouse_ChIP-seq_mESCsShContD0_DEGseq_FC.txt

awk 'BEGIN{OFS="\t";FS="\t";}
    {
        if(FNR>1){
            if($5>0 && $5!="NA"){
                print $1, $5;
            }
        }
    }
' mESCsShSetD2D0/output_score.txt > mouse_ChIP-seq_mESCsShSetD2D0_DEGseq_FC.txt

awk 'BEGIN{OFS="\t";FS="\t";}
    {
        if(FNR>1){
            if($5>0 && $5!="NA"){
                print $1, $5;
            }
        }
    }
' mEFsShCont/output_score.txt > mouse_ChIP-seq_mEFsShCont_DEGseq_FC.txt

awk 'BEGIN{OFS="\t";FS="\t";}
    {
        if(FNR>1){
            if($5>0 && $5!="NA"){
                print $1, $5;
            }
        }
    }
' mEFsShSetD2/output_score.txt > mouse_ChIP-seq_mEFsShSetD2_DEGseq_FC.txt

cat *_DEGseq_FC.txt | cut -f 1 | sort | uniq > mouse_ChIP-seq_union_gene.txt

awk 'BEGIN{OFS="\t";FS="\t";}
ARGIND==1{
  hashGene[$1]=1;
  hashGeneExp[$1,1]="NA";
  hashGeneExp[$1,2]="NA";
  hashGeneExp[$1,3]="NA";
  hashGeneExp[$1,4]="NA";
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
END{
  for (gene in hashGene) {
    print gene, hashGeneExp[gene,1], hashGeneExp[gene,2],
      hashGeneExp[gene,3], hashGeneExp[gene,4];
  }
}
' mouse_ChIP-seq_union_gene.txt mouse_ChIP-seq_mESCsShContD0_DEGseq_FC.txt \
  mouse_ChIP-seq_mESCsShSetD2D0_DEGseq_FC.txt mouse_ChIP-seq_mEFsShCont_DEGseq_FC.txt \
  mouse_ChIP-seq_mEFsShSetD2_DEGseq_FC.txt > mouse_ChIP-seq_DEGseq_union_FC.txt

cat *_DEGseq_FC.txt | cut -f 1 | sort | uniq -c | awk '$1 == 4 { print $2 }' > mouse_ChIP-seq_intersected_gene.txt

awk 'BEGIN{OFS="\t";FS="\t";}
ARGIND==1{
  hashGene[$1]=1;
  hashGeneExp[$1,1]="NA";
  hashGeneExp[$1,2]="NA";
  hashGeneExp[$1,3]="NA";
  hashGeneExp[$1,4]="NA";
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
END{
  for (gene in hashGene) {
    print gene, hashGeneExp[gene,1], hashGeneExp[gene,2],
      hashGeneExp[gene,3], hashGeneExp[gene,4];
  }
}
' mouse_ChIP-seq_intersected_gene.txt mouse_ChIP-seq_mESCsShContD0_DEGseq_FC.txt \
  mouse_ChIP-seq_mESCsShSetD2D0_DEGseq_FC.txt mouse_ChIP-seq_mEFsShCont_DEGseq_FC.txt \
  mouse_ChIP-seq_mEFsShSetD2_DEGseq_FC.txt > mouse_ChIP-seq_DEGseq_intersected_FC.txt

##mouse_ChIP-seq_DEGseq_FC.txt
##GeneID,mESCsShContD0,mESCsShSetD2D0,mEFsShCont,mEFsShSetD2
