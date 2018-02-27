#!/bin/sh
gtf=/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.vM9.annotation.gtf
mESCsbamPath=/data/zhoukr/hhl_setd2_m6a/mouse_ESCs_m6A-seq
mEFsbamPath=/data/zhoukr/hhl_setd2_m6a/mouse_MEFs_m6A-seq
rscript=/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/Rscripts/mouse/mouse_DEGseq.r
ouput=/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/mouse/DEGseq/
bed6=/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.vM9.annotation.all.longest.exon.bed6

rm -rf $ouput
mkdir -p $ouput/DeGenes
cd $ouput

featureCounts -T 10 -a $gtf -g gene_id -F GTF -t gene -M -o temp.txt \
  $mESCsbamPath/Ctrl_D0/mESCs_m6A-seq_Ctrl_D0_IP.sorted.bam \
  $mESCsbamPath/Ctrl_D0/mESCs_m6A-seq_Ctrl_D0_input.sorted.bam \
  $mESCsbamPath/Ctrl_D6/mESCs_m6A-seq_Ctrl_D6_IP.sorted.bam \
  $mESCsbamPath/Ctrl_D6/mESCs_m6A-seq_Ctrl_D6_input.sorted.bam \
  $mESCsbamPath/SetD2-KD_D0/mESCs_m6A-seq_SetD2-KD_D0_IP.sorted.bam \
  $mESCsbamPath/SetD2-KD_D0/mESCs_m6A-seq_SetD2-KD_D0_input.sorted.bam \
  $mESCsbamPath/SetD2-KD_D6/mESCs_m6A-seq_SetD2-KD_D6_IP.sorted.bam \
  $mESCsbamPath/SetD2-KD_D6/mESCs_m6A-seq_SetD2-KD_D6_input.sorted.bam \
  $mEFsbamPath/SetD2-WT_NA/mMEFs_m6A-seq_SetD2-WT_NA_IP.sorted.bam \
  $mEFsbamPath/SetD2-WT_NA/mMEFs_m6A-seq_SetD2-WT_NA_input.sorted.bam \
  $mEFsbamPath/SetD2-KO_NA/mMEFs_m6A-seq_SetD2-KO_NA_IP.sorted.bam \
  $mEFsbamPath/SetD2-KO_NA/mMEFs_m6A-seq_SetD2-KO_NA_input.sorted.bam \
   > /dev/null 2>&1

awk 'BEGIN{OFS="\t";FS="\t";}
    {
        if(FNR>1){
            if(FNR==2){
                print $1,$2,$3,$4,$5,$6,
                "mESCsShContD0IP","mESCsShContD0Input","mESCsShContD6IP","mESCsShContD6Input",
                "mESCsShSetD2D0IP","mESCsShSetD2D0Input","mESCsShSetD2D6IP","mESCsShSetD2D6Input",
                "mEFsShContIP","mEFsShContInput", "mEFsShSetD2IP","mEFsShSetD2Input";
            }else{
                print $0
            }
        }
    }
' temp.txt > mouse_RNA-seq_featureCount.txt

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
' mESCsShContD0/output_score.txt > mouse_RNA-seq_mESCsShContD0_DEGseq_FC.txt

awk 'BEGIN{OFS="\t";FS="\t";}
    {
        if(FNR>1){
            if($5>0 && $5!="NA"){
                print $1, $5;
            }
        }
    }
' mESCsShContD6/output_score.txt > mouse_RNA-seq_mESCsShContD6_DEGseq_FC.txt

awk 'BEGIN{OFS="\t";FS="\t";}
    {
        if(FNR>1){
            if($5>0 && $5!="NA"){
                print $1, $5;
            }
        }
    }
' mESCsShSetD2D0/output_score.txt > mouse_RNA-seq_mESCsShSetD2D0_DEGseq_FC.txt

awk 'BEGIN{OFS="\t";FS="\t";}
    {
        if(FNR>1){
            if($5>0 && $5!="NA"){
                print $1, $5;
            }
        }
    }
' mESCsShSetD2D6/output_score.txt > mouse_RNA-seq_mESCsShSetD2D6_DEGseq_FC.txt

awk 'BEGIN{OFS="\t";FS="\t";}
    {
        if(FNR>1){
            if($5>0 && $5!="NA"){
                print $1, $5;
            }
        }
    }
' mEFsShCont/output_score.txt > mouse_RNA-seq_mEFsShCont_DEGseq_FC.txt

awk 'BEGIN{OFS="\t";FS="\t";}
    {
        if(FNR>1){
            if($5>0 && $5!="NA"){
                print $1, $5;
            }
        }
    }
' mEFsShSetD2/output_score.txt > mouse_RNA-seq_mEFsShSetD2_DEGseq_FC.txt

cat *_DEGseq_FC.txt | cut -f 1 | sort | uniq > mouse_RNA-seq_union_gene.txt

awk 'BEGIN{OFS="\t";FS="\t";}
ARGIND==1{
  hashGene[$1]=1;
  hashGeneExp[$1,1]="NA";
  hashGeneExp[$1,2]="NA";
  hashGeneExp[$1,3]="NA";
  hashGeneExp[$1,4]="NA";
  hashGeneExp[$1,5]="NA";
  hashGeneExp[$1,6]="NA";
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
ARGIND==7{
  hashGeneExp[$1,6]=$2;
}
END{
  for (gene in hashGene) {
    print gene, hashGeneExp[gene,1], hashGeneExp[gene,2],
      hashGeneExp[gene,3], hashGeneExp[gene,4], hashGeneExp[gene,5], hashGeneExp[gene,6];
  }
}
' mouse_RNA-seq_union_gene.txt mouse_RNA-seq_mESCsShContD0_DEGseq_FC.txt \
  mouse_RNA-seq_mESCsShContD6_DEGseq_FC.txt mouse_RNA-seq_mESCsShSetD2D0_DEGseq_FC.txt \
  mouse_RNA-seq_mESCsShSetD2D6_DEGseq_FC.txt mouse_RNA-seq_mEFsShCont_DEGseq_FC.txt \
  mouse_RNA-seq_mEFsShSetD2_DEGseq_FC.txt > mouse_RNA-seq_DEGseq_union_FC.txt

cat *_DEGseq_FC.txt | cut -f 1 | sort | uniq -c | awk '$1 == 6 { print $2 }' > mouse_RNA-seq_intersected_gene.txt

awk 'BEGIN{OFS="\t";FS="\t";}
ARGIND==1{
  hashGene[$1]=1;
  hashGeneExp[$1,1]="NA";
  hashGeneExp[$1,2]="NA";
  hashGeneExp[$1,3]="NA";
  hashGeneExp[$1,4]="NA";
  hashGeneExp[$1,5]="NA";
  hashGeneExp[$1,6]="NA";
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
ARGIND==7{
  hashGeneExp[$1,6]=$2;
}
END{
  for (gene in hashGene) {
    print gene, hashGeneExp[gene,1], hashGeneExp[gene,2],
      hashGeneExp[gene,3], hashGeneExp[gene,4], hashGeneExp[gene,5], hashGeneExp[gene,6];
  }
}
' mouse_RNA-seq_intersected_gene.txt mouse_RNA-seq_mESCsShContD0_DEGseq_FC.txt \
  mouse_RNA-seq_mESCsShContD6_DEGseq_FC.txt mouse_RNA-seq_mESCsShSetD2D0_DEGseq_FC.txt \
  mouse_RNA-seq_mESCsShSetD2D6_DEGseq_FC.txt mouse_RNA-seq_mEFsShCont_DEGseq_FC.txt \
  mouse_RNA-seq_mEFsShSetD2_DEGseq_FC.txt > mouse_RNA-seq_DEGseq_intersected_FC.txt

##mouse_RNA-seq_DEGseq_FC.txt
##GeneID,mESCsShContD0,mESCsShContD6,mESCsShSetD2D0,mESCsShSetD2D6,mEFsShCont,mEFsShSetD2

cd DeGenes/

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
}' $bed6 mESCsShContD6/output_score.txt > mouse_RNA-seq_mESCsShContD6_DEGseq_DeGenes.txt

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
}' $bed6 mESCsShSetD2D0/output_score.txt > mouse_RNA-seq_mESCsShSetD2D0_DEGseq_DeGenes.txt

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
}' $bed6 mESCsShSetD2D6/output_score.txt > mouse_RNA-seq_mESCsShSetD2D6_DEGseq_DeGenes.txt

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
}' $bed6 mEFsShSetD2/output_score.txt > mouse_RNA-seq_mEFsShSetD2_DEGseq_DeGenes.txt

