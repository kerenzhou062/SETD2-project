#!/bin/sh

### Hela, ChIP-seq
cd /data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/HepG2/input
annotation=/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.all.longest.exon.bed6

echo -e "GeneID\tgeneName\tshSETD2\tshMETTL14" > HepG2_RNA-seq_shSETD2vsshMETTL14_gene_FC.txt
awk 'BEGIN{OFS="\t";}
    ARGIND==1{
        split($4,geneIdI,":");
        geneIdi=geneIdI[1];
        geneNamei=geneIdI[2];
        hashGeneIdName[geneIdi]=geneNamei;
    }
    ARGIND==2{
        if(FNR>1){
            if($6<=0.05){
                hash[$2]= $3;
            }
        }
    }
    ARGIND==3{
        if(FNR>1){
            if($6<=0.05){
                if($2 in hash){

                    print $2,hashGeneIdName[$2],hash[$2],$3;
                }
            }
        }
    }' $annotation HepG2_RNA-seq_shCont-shSetD2_htseq-count-gene-edgeR.FC.txt HepG2_RNA-seq_shCont-shM14_htseq-count-gene-edgeR.FC.txt >> HepG2_RNA-seq_shSETD2vsshMETTL14_gene_FC.txt

echo -e "GeneID\tgeneName\tshSETD2\tshMETTL3" > HepG2_RNA-seq_shSETD2vsshMETTL3_gene_FC.txt
awk 'BEGIN{OFS="\t";}
    ARGIND==1{
        split($4,geneIdI,":");
        geneIdi=geneIdI[1];
        geneNamei=geneIdI[2];
        hashGeneIdName[geneIdi]=geneNamei;
    }
    ARGIND==2{
        if(FNR>1){
            if($6<=0.05){
                hash[$2]= $3;
            }
        }
    }
    ARGIND==3{
        if(FNR>1){
            if($6<=0.05){
                if($2 in hash){
                    print $2,hashGeneIdName[$2],hash[$2],$3;
                }
            }
        }
    }' $annotation HepG2_RNA-seq_shCont-shSetD2_htseq-count-gene-edgeR.FC.txt HepG2_RNA-seq_shCont-shM3_htseq-count-gene-edgeR.FC.txt >> HepG2_RNA-seq_shSETD2vsshMETTL3_gene_FC.txt

echo -e "GeneID\tgeneName\tshSETD2\tshWTAP" > HepG2_RNA-seq_shSETD2vsshWTAP_gene_FC.txt
awk 'BEGIN{OFS="\t";}
    ARGIND==1{
        split($4,geneIdI,":");
        geneIdi=geneIdI[1];
        geneNamei=geneIdI[2];
        hashGeneIdName[geneIdi]=geneNamei;
    }
    ARGIND==2{
        if(FNR>1){
            if($6<=0.05){
                hash[$2]= $3;
            }
        }
    }
    ARGIND==3{
        if(FNR>1){
            if($6<=0.05){
                if($2 in hash){
                    print $2,hashGeneIdName[$2],hash[$2],$3;
                }
            }
        }
    }' $annotation HepG2_RNA-seq_shCont-shSetD2_htseq-count-gene-edgeR.FC.txt HepG2_RNA-seq_shCont-shWTAP_htseq-count-gene-edgeR.FC.txt >> HepG2_RNA-seq_shSETD2vsshWTAP_gene_FC.txt
