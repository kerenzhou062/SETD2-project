#!/bin/sh

### Hela, ChIP-seq
cd /data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/mMEFs/IPvsInput
annotation=/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.vM9.annotation.all.longest.exon.bed6
geneExpDir=/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/mMEFs/RSEM

### Gene M6A level assessment
echo -e "GeneID\tgeneName-lowercase\tgeneName-uppercase\tInput(FPKM)\tIP(FPKM)\tGeneM6ALevel\tGeneM6ALevel(log2((IP+0.01)/(Input+0.01)))" > mMEFs_gene_m6A_level_shCont.txt
awk 'BEGIN{OFS="\t";}
    ARGIND==1{
        split($4,geneIdI,":");
        geneIdi=geneIdI[1];
        geneNamei=geneIdI[2];
        hashGeneIdName[geneIdi]=geneNamei;
    }
    ARGIND==2{
        if(FNR>1){
            EXPa[$1]= $7;
        }
    }
    ARGIND==3{
        if(FNR>1){
            EXPb[$1]= $7;
        }
    }
    END{
        for(geneID in EXPa){
            ExpA=EXPa[geneID];ExpB=EXPb[geneID];
            GeneM6ALevel=log((ExpB+0.01)/(ExpA+0.01))/log(2);
            geneName=hashGeneIdName[geneID];
            print geneID,geneName,toupper(geneName),ExpA,ExpB,2**GeneM6ALevel,GeneM6ALevel;
        }
    }' $annotation $geneExpDir/SetD2-WT_NA_input.genes.results $geneExpDir/SetD2-WT_NA_IP.genes.results >> mMEFs_gene_m6A_level_shCont.txt

echo -e "GeneID\tgeneName-lowercase\tgeneName-uppercase\tInput(FPKM)\tIP(FPKM)\tGeneM6ALevel\tGeneM6ALevel(log2((IP+0.01)/(Input+0.01)))" > mMEFs_gene_m6A_level_shSetD2.txt
awk 'BEGIN{OFS="\t";}
    ARGIND==1{
        split($4,geneIdI,":");
        geneIdi=geneIdI[1];
        geneNamei=geneIdI[2];
        hashGeneIdName[geneIdi]=geneNamei;
    }
    ARGIND==2{
        if(FNR>1){
            EXPa[$1]= $7;
        }
    }
    ARGIND==3{
        if(FNR>1){
            EXPb[$1]= $7;
        }
    }
    END{
        for(geneID in EXPa){
            ExpA=EXPa[geneID];ExpB=EXPb[geneID];
            GeneM6ALevel=log((ExpB+0.01)/(ExpA+0.01))/log(2);
            geneName=hashGeneIdName[geneID];
            print geneID,geneName,toupper(geneName),ExpA,ExpB,2**GeneM6ALevel,GeneM6ALevel;
        }
    }' $annotation $geneExpDir/SetD2-KO_NA_input.genes.results $geneExpDir/SetD2-KO_NA_IP.genes.results >> mMEFs_gene_m6A_level_shSetD2.txt
