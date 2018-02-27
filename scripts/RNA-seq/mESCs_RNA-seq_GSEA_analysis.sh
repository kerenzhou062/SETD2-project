#!/bin/sh

### mESCs, RNA-seq
cd /data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/mESCs/GSEA/
geneExpDir=/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/mESCs/RSEM
annotation=/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.vM9.annotation.all.longest.exon.bed6

## .cls file
echo -e "4\t4\t1\n#shContD0\tshContD6\tshSetD2D0\tshSetD2D6\nshContD0\tshContD6\tshSetD2D0\tshSetD2D6" > mESCs_Phenotype.cls

## TXT: Text file format for expression dataset (*.txt)
echo -e "NAME\tDESCRIPTION\tshContD0\tshContD6\tshSetD2D0\tshSetD2D6" > mESCs_gene_level_GSEA_geneUppercase.txt
awk 'BEGIN{OFS="\t";}
    ARGIND==1{
        split($4,geneIdI,":");
        geneIdi=geneIdI[1];
        geneNamei=geneIdI[2];
        geneType=geneIdI[3];
        if(geneType=="protein_coding"){
            hashGeneIdName[geneIdi]=geneNamei;
        }
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
    ARGIND==4{
        if(FNR>1){
            EXPc[$1]= $7;
        }
    }
    ARGIND==5{
        if(FNR>1){
            EXPd[$1]= $7;
        }
    }
    END{
        slen=asorti(EXPa,EXPaS);
        for (i=1;i<=slen;i++){
            geneID=EXPaS[i];
            if(geneID in hashGeneIdName){
                ExpA=EXPa[EXPaS[i]];ExpB=EXPb[geneID];ExpC=EXPc[geneID];ExpD=EXPd[geneID];
                geneName=hashGeneIdName[geneID];
                if(ExpA<1&&ExpB<1&&ExpC<1&&ExpD<1){}else{
                    print toupper(geneName),"na",ExpA,ExpB,ExpC,ExpD;
                }
            }
        }
    }' $annotation $geneExpDir/Ctrl_D0_input.genes.results $geneExpDir/Ctrl_D6_input.genes.results \
    $geneExpDir/SetD2-KD_D0_input.genes.results $geneExpDir/SetD2-KD_D6_input.genes.results >> mESCs_gene_level_GSEA_geneUppercase.txt

echo -e "NAME\tDESCRIPTION\tshContD0\tshContD6\tshSetD2D0\tshSetD2D6" > mESCs_gene_level_GSEA_geneLowercase.txt
awk 'BEGIN{OFS="\t";}
    ARGIND==1{
        split($4,geneIdI,":");
        geneIdi=geneIdI[1];
        geneNamei=geneIdI[2];
        geneType=geneIdI[3];
        if(geneType=="protein_coding"){
            hashGeneIdName[geneIdi]=geneNamei;
        }
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
    ARGIND==4{
        if(FNR>1){
            EXPc[$1]= $7;
        }
    }
    ARGIND==5{
        if(FNR>1){
            EXPd[$1]= $7;
        }
    }
    END{
        slen=asorti(EXPa,EXPaS);
        for (i=1;i<=slen;i++){
            geneID=EXPaS[i];
            if(geneID in hashGeneIdName){
                ExpA=EXPa[EXPaS[i]];ExpB=EXPb[geneID];ExpC=EXPc[geneID];ExpD=EXPd[geneID];
                geneName=hashGeneIdName[geneID];
                if(ExpA<1&&ExpB<1&&ExpC<1&&ExpD<1){}else{
                    print tolower(geneName),"na",ExpA,ExpB,ExpC,ExpD;
                }
            }
        }
    }' $annotation $geneExpDir/Ctrl_D0_input.genes.results $geneExpDir/Ctrl_D6_input.genes.results \
    $geneExpDir/SetD2-KD_D0_input.genes.results $geneExpDir/SetD2-KD_D6_input.genes.results >> mESCs_gene_level_GSEA_geneLowercase.txt

echo -e "NAME\tDESCRIPTION\tshContD0\tshContD6\tshSetD2D0\tshSetD2D6" > mESCs_gene_m6a_tagsEnrichment_GSEA_geneUppercase.txt
awk 'BEGIN{OFS="\t";}
    ARGIND==1{
        split($4,geneIdI,":");
        geneIdi=geneIdI[1];
        geneNamei=geneIdI[2];
        geneType=geneIdI[3];
        if(geneType=="protein_coding"){
            hashGeneIdName[geneIdi]=geneNamei;
        }
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
    ARGIND==4{
        if(FNR>1){
            EXPc[$1]= $7;
        }
    }
    ARGIND==5{
        if(FNR>1){
            EXPd[$1]= $7;
        }
    }
    END{
        slen=asorti(EXPa,EXPaS);
        for (i=1;i<=slen;i++){
            geneID=EXPaS[i];
            if(geneID in hashGeneIdName){
                ExpA=EXPa[EXPaS[i]];ExpB=EXPb[geneID];ExpC=EXPc[geneID];ExpD=EXPd[geneID];
                geneName=hashGeneIdName[geneID];
                if(ExpA<1&&ExpB<1&&ExpC<1&&ExpD<1){}else{
                    print toupper(geneName),"na",ExpA,ExpB,ExpC,ExpD;
                }
            }
        }
    }' $annotation $geneExpDir/Ctrl_D0_IP.genes.results $geneExpDir/Ctrl_D6_IP.genes.results \
    $geneExpDir/SetD2-KD_D0_IP.genes.results $geneExpDir/SetD2-KD_D6_IP.genes.results >> mESCs_gene_m6a_tagsEnrichment_GSEA_geneUppercase.txt

echo -e "NAME\tDESCRIPTION\tshContD0\tshContD6\tshSetD2D0\tshSetD2D6" > mESCs_gene_m6a_tagsEnrichment_GSEA_geneLowercase.txt
awk 'BEGIN{OFS="\t";}
    ARGIND==1{
        split($4,geneIdI,":");
        geneIdi=geneIdI[1];
        geneNamei=geneIdI[2];
        geneType=geneIdI[3];
        if(geneType=="protein_coding"){
            hashGeneIdName[geneIdi]=geneNamei;
        }
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
    ARGIND==4{
        if(FNR>1){
            EXPc[$1]= $7;
        }
    }
    ARGIND==5{
        if(FNR>1){
            EXPd[$1]= $7;
        }
    }
    END{
        slen=asorti(EXPa,EXPaS);
        for (i=1;i<=slen;i++){
            geneID=EXPaS[i];
            if(geneID in hashGeneIdName){
                ExpA=EXPa[EXPaS[i]];ExpB=EXPb[geneID];ExpC=EXPc[geneID];ExpD=EXPd[geneID];
                geneName=hashGeneIdName[geneID];
                if(ExpA<1&&ExpB<1&&ExpC<1&&ExpD<1){}else{
                    print tolower(geneName),"na",ExpA,ExpB,ExpC,ExpD;
                }
            }
        }
    }' $annotation $geneExpDir/Ctrl_D0_IP.genes.results $geneExpDir/Ctrl_D6_IP.genes.results \
    $geneExpDir/SetD2-KD_D0_IP.genes.results $geneExpDir/SetD2-KD_D6_IP.genes.results >> mESCs_gene_m6a_tagsEnrichment_GSEA_geneLowercase.txt
