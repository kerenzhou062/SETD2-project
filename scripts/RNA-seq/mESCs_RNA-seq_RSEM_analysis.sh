#!/bin/sh

### mESCs, RNA-seq
cd /data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/mESCs/RSEM/
annotation=/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.vM9.annotation.all.longest.exon.bed6

echo -e "GeneID\tgeneName\tshCont-D0(FPKM)\tshSetD2-D0(FPKM)" > mESCs_RNA-seq_RSEM_D0_input.txt
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
    }' $annotation Ctrl_D0_input.genes.results SetD2-KD_D0_input.genes.results >> mESCs_RNA-seq_RSEM_D0_input.txt

echo -e "GeneID\tgeneName\tshCont-D6(FPKM)\tshSetD2-D6(FPKM)" > mESCs_RNA-seq_RSEM_D6_input.txt
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
    }' $annotation Ctrl_D6_input.genes.results SetD2-KD_D6_input.genes.results >> mESCs_RNA-seq_RSEM_D6_input.txt

echo -e "GeneID\tgeneName\tshCont-D0(FPKM)\tshSetD2-D0(FPKM)" > mESCs_RNA-seq_RSEM_D0_IP.txt
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
    }' $annotation Ctrl_D0_IP.genes.results SetD2-KD_D0_IP.genes.results >> mESCs_RNA-seq_RSEM_D0_IP.txt

echo -e "GeneID\tgeneName\tshCont-D6(FPKM)\tshSetD2-D6(FPKM)" > mESCs_RNA-seq_RSEM_D6_IP.txt
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
    }' $annotation Ctrl_D6_IP.genes.results SetD2-KD_D6_IP.genes.results >> mESCs_RNA-seq_RSEM_D6_IP.txt

