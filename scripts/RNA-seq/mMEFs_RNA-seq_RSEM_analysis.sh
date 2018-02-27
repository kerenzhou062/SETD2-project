#!/bin/sh

### mMEFs, RNA-seq
cd /data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/mMEFs/RSEM
annotation=/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.vM9.annotation.all.longest.exon.bed6

echo -e "GeneID\tgeneName\tshCont(FPKM)\tshSetD2(FPKM)" > mMEFs_RNA-seq_RSEM_input.txt
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
        slen=asorti(EXPa,EXPaS);
        for (i=1;i<=slen;i++){
            geneID=EXPaS[i];
            if(geneID in hashGeneIdName){
                ExpA=EXPa[EXPaS[i]];ExpB=EXPb[geneID];
                geneName=hashGeneIdName[geneID];
                print geneID,geneName,ExpA,ExpB;
            }
        }
    }' $annotation SetD2-WT_NA_input.genes.results SetD2-KO_NA_input.genes.results >> mMEFs_RNA-seq_RSEM_input.txt

echo -e "GeneID\tgeneName\tshCont(FPKM)\tshSetD2(FPKM)" > mMEFs_RNA-seq_RSEM_IP.txt
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
        slen=asorti(EXPa,EXPaS);
        for (i=1;i<=slen;i++){
            geneID=EXPaS[i];
            if(geneID in hashGeneIdName){
                ExpA=EXPa[EXPaS[i]];ExpB=EXPb[geneID];
                geneName=hashGeneIdName[geneID];
                print geneID,geneName,ExpA,ExpB;
            }
        }
    }' $annotation SetD2-WT_NA_IP.genes.results SetD2-KO_NA_IP.genes.results >> mMEFs_RNA-seq_RSEM_IP.txt
