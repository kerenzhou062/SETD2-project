#!/bin/sh

### Hela, RNA-seq
cd /data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/Hela/RSEM/
annotation=/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.all.longest.exon.bed6

echo -e "GeneID\tgeneName\trep1\trep2\trep3\tMeanFPKM\tlog2MeanFPKM\trep1_log2\trep2_log2\trep3_log2" > Hela_RNA-seq_RSEM_shCont_input.txt
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
    ARGIND==4{
        if(FNR>1){
            EXPc[$1]= $7;
        }
    }
    END{
        slen=asorti(EXPa,EXPaS);
        for (i=1;i<=slen;i++){
            geneID=EXPaS[i];
            if(geneID in hashGeneIdName){
                ExpA=EXPa[EXPaS[i]];ExpB=EXPb[geneID];ExpC=EXPc[geneID];
                ExpMean=(ExpA+ExpB+ExpC)/3;
                geneName=hashGeneIdName[geneID];
                if (ExpMean>=0){
                    print geneID,geneName,ExpA,ExpB,ExpC,ExpMean,log(ExpMean)/log(2),log(ExpA)/log(2),log(ExpB)/log(2),log(ExpC)/log(2);
                }
            }
        }
    }' $annotation shCont_input_rep1.genes.results shCont_input_rep2.genes.results shCont_input_rep3.genes.results >> Hela_RNA-seq_RSEM_shCont_input.txt

echo -e "GeneID\tgeneName\trep1\trep2\trep3\tMeanFPKM\tlog2MeanFPKM\trep1_log2\trep2_log2\trep3_log2" > Hela_RNA-seq_RSEM_shCont_IP.txt
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
    ARGIND==4{
        if(FNR>1){
            EXPc[$1]= $7;
        }
    }
    END{
        slen=asorti(EXPa,EXPaS);
        for (i=1;i<=slen;i++){
            geneID=EXPaS[i];
            if(geneID in hashGeneIdName){
                ExpA=EXPa[EXPaS[i]];ExpB=EXPb[geneID];ExpC=EXPc[geneID];
                ExpMean=(ExpA+ExpB+ExpC)/3;
                geneName=hashGeneIdName[geneID];
                if (ExpMean>=0){
                    print geneID,geneName,ExpA,ExpB,ExpC,ExpMean,log(ExpMean)/log(2),log(ExpA)/log(2),log(ExpB)/log(2),log(ExpC)/log(2);
                }
            }
        }
    }' $annotation shCont_IP_rep1.genes.results shCont_IP_rep2.genes.results shCont_IP_rep3.genes.results >> Hela_RNA-seq_RSEM_shCont_IP.txt

echo -e "GeneID\tgeneName\trep1\trep2\trep3\tMeanFPKM\tlog2MeanFPKM\trep1_log2\trep2_log2\trep3_log2" > Hela_RNA-seq_RSEM_shSetD2_input.txt
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
    ARGIND==4{
        if(FNR>1){
            EXPc[$1]= $7;
        }
    }
    END{
        slen=asorti(EXPa,EXPaS);
        for (i=1;i<=slen;i++){
            geneID=EXPaS[i];
            if(geneID in hashGeneIdName){
                ExpA=EXPa[EXPaS[i]];ExpB=EXPb[geneID];ExpC=EXPc[geneID];
                ExpMean=(ExpA+ExpB+ExpC)/3;
                geneName=hashGeneIdName[geneID];
                if (ExpMean>=0){
                    print geneID,geneName,ExpA,ExpB,ExpC,ExpMean,log(ExpMean)/log(2),log(ExpA)/log(2),log(ExpB)/log(2),log(ExpC)/log(2);
                }
            }
        }
    }' $annotation shSetD2_input_rep1.genes.results shSetD2_input_rep2.genes.results shSetD2_input_rep3.genes.results >> Hela_RNA-seq_RSEM_shSetD2_input.txt

echo -e "GeneID\tgeneName\trep1\trep2\trep3\tMeanFPKM\tlog2MeanFPKM\trep1_log2\trep2_log2\trep3_log2" > Hela_RNA-seq_RSEM_shSetD2_IP.txt
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
    ARGIND==4{
        if(FNR>1){
            EXPc[$1]= $7;
        }
    }
    END{
        slen=asorti(EXPa,EXPaS);
        for (i=1;i<=slen;i++){
            geneID=EXPaS[i];
            if(geneID in hashGeneIdName){
                ExpA=EXPa[EXPaS[i]];ExpB=EXPb[geneID];ExpC=EXPc[geneID];
                ExpMean=(ExpA+ExpB+ExpC)/3;
                geneName=hashGeneIdName[geneID];
                if (ExpMean>=0){
                    print geneID,geneName,ExpA,ExpB,ExpC,ExpMean,log(ExpMean)/log(2),log(ExpA)/log(2),log(ExpB)/log(2),log(ExpC)/log(2);
                }
            }
        }
    }' $annotation shSetD2_IP_rep1.genes.results shSetD2_IP_rep2.genes.results shSetD2_IP_rep3.genes.results >> Hela_RNA-seq_RSEM_shSetD2_IP.txt

echo -e "GeneID\tgeneName\trep1\trep2\trep3\tMeanFPKM\tlog2MeanFPKM\trep1_log2\trep2_log2\trep3_log2" > Hela_RNA-seq_RSEM_shM14_input.txt
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
    ARGIND==4{
        if(FNR>1){
            EXPc[$1]= $7;
        }
    }
    END{
        slen=asorti(EXPa,EXPaS);
        for (i=1;i<=slen;i++){
            geneID=EXPaS[i];
            if(geneID in hashGeneIdName){
                ExpA=EXPa[EXPaS[i]];ExpB=EXPb[geneID];ExpC=EXPc[geneID];
                ExpMean=(ExpA+ExpB+ExpC)/3;
                geneName=hashGeneIdName[geneID];
                if (ExpMean>=0){
                    print geneID,geneName,ExpA,ExpB,ExpC,ExpMean,log(ExpMean)/log(2),log(ExpA)/log(2),log(ExpB)/log(2),log(ExpC)/log(2);
                }
            }
        }
    }' $annotation shM14_input_rep1.genes.results shM14_input_rep2.genes.results shM14_input_rep3.genes.results >> Hela_RNA-seq_RSEM_shM14_input.txt

echo -e "GeneID\tgeneName\trep1\trep2\trep3\tMeanFPKM\tlog2MeanFPKM\trep1_log2\trep2_log2\trep3_log2" > Hela_RNA-seq_RSEM_shM14_IP.txt
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
    ARGIND==4{
        if(FNR>1){
            EXPc[$1]= $7;
        }
    }
    END{
        slen=asorti(EXPa,EXPaS);
        for (i=1;i<=slen;i++){
            geneID=EXPaS[i];
            if(geneID in hashGeneIdName){
                ExpA=EXPa[EXPaS[i]];ExpB=EXPb[geneID];ExpC=EXPc[geneID];
                ExpMean=(ExpA+ExpB+ExpC)/3;
                geneName=hashGeneIdName[geneID];
                if (ExpMean>=0){
                    print geneID,geneName,ExpA,ExpB,ExpC,ExpMean,log(ExpMean)/log(2),log(ExpA)/log(2),log(ExpB)/log(2),log(ExpC)/log(2);
                }
            }
        }
    }' $annotation shM14_IP_rep1.genes.results shM14_IP_rep2.genes.results shM14_IP_rep3.genes.results >> Hela_RNA-seq_RSEM_shM14_IP.txt

echo -e "GeneID\tgeneName\trep1\trep2\trep3\tMeanFPKM\tlog2MeanFPKM\trep1_log2\trep2_log2\trep3_log2" > Hela_RNA-seq_RSEM_shM3_input.txt
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
    ARGIND==4{
        if(FNR>1){
            EXPc[$1]= $7;
        }
    }
    END{
        slen=asorti(EXPa,EXPaS);
        for (i=1;i<=slen;i++){
            geneID=EXPaS[i];
            if(geneID in hashGeneIdName){
                ExpA=EXPa[EXPaS[i]];ExpB=EXPb[geneID];ExpC=EXPc[geneID];
                ExpMean=(ExpA+ExpB+ExpC)/3;
                geneName=hashGeneIdName[geneID];
                if (ExpMean>=0){
                    print geneID,geneName,ExpA,ExpB,ExpC,ExpMean,log(ExpMean)/log(2),log(ExpA)/log(2),log(ExpB)/log(2),log(ExpC)/log(2);
                }
            }
        }
    }' $annotation shM3_input_rep1.genes.results shM3_input_rep2.genes.results shM3_input_rep3.genes.results >> Hela_RNA-seq_RSEM_shM3_input.txt

echo -e "GeneID\tgeneName\trep1\trep2\trep3\tMeanFPKM\tlog2MeanFPKM\trep1_log2\trep2_log2\trep3_log2" > Hela_RNA-seq_RSEM_shM3_IP.txt
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
    ARGIND==4{
        if(FNR>1){
            EXPc[$1]= $7;
        }
    }
    END{
        slen=asorti(EXPa,EXPaS);
        for (i=1;i<=slen;i++){
            geneID=EXPaS[i];
            if(geneID in hashGeneIdName){
                ExpA=EXPa[EXPaS[i]];ExpB=EXPb[geneID];ExpC=EXPc[geneID];
                ExpMean=(ExpA+ExpB+ExpC)/3;
                geneName=hashGeneIdName[geneID];
                if (ExpMean>=0){
                    print geneID,geneName,ExpA,ExpB,ExpC,ExpMean,log(ExpMean)/log(2),log(ExpA)/log(2),log(ExpB)/log(2),log(ExpC)/log(2);
                }
            }
        }
    }' $annotation shM3_IP_rep1.genes.results shM3_IP_rep2.genes.results shM3_IP_rep3.genes.results >> Hela_RNA-seq_RSEM_shM3_IP.txt

echo -e "GeneID\tgeneName\trep1\trep2\trep3\tMeanFPKM\tlog2MeanFPKM\trep1_log2\trep2_log2\trep3_log2" > Hela_RNA-seq_RSEM_shWTAP_input.txt
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
    ARGIND==4{
        if(FNR>1){
            EXPc[$1]= $7;
        }
    }
    END{
        slen=asorti(EXPa,EXPaS);
        for (i=1;i<=slen;i++){
            geneID=EXPaS[i];
            if(geneID in hashGeneIdName){
                ExpA=EXPa[EXPaS[i]];ExpB=EXPb[geneID];ExpC=EXPc[geneID];
                ExpMean=(ExpA+ExpB+ExpC)/3;
                geneName=hashGeneIdName[geneID];
                if (ExpMean>=0){
                    print geneID,geneName,ExpA,ExpB,ExpC,ExpMean,log(ExpMean)/log(2),log(ExpA)/log(2),log(ExpB)/log(2),log(ExpC)/log(2);
                }
            }
        }
    }' $annotation shWTAP_input_rep1.genes.results shWTAP_input_rep2.genes.results shWTAP_input_rep3.genes.results >> Hela_RNA-seq_RSEM_shWTAP_input.txt

echo -e "GeneID\tgeneName\trep1\trep2\trep3\tMeanFPKM\tlog2MeanFPKM\trep1_log2\trep2_log2\trep3_log2" > Hela_RNA-seq_RSEM_shWTAP_IP.txt
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
    ARGIND==4{
        if(FNR>1){
            EXPc[$1]= $7;
        }
    }
    END{
        slen=asorti(EXPa,EXPaS);
        for (i=1;i<=slen;i++){
            geneID=EXPaS[i];
            if(geneID in hashGeneIdName){
                ExpA=EXPa[EXPaS[i]];ExpB=EXPb[geneID];ExpC=EXPc[geneID];
                ExpMean=(ExpA+ExpB+ExpC)/3;
                geneName=hashGeneIdName[geneID];
                if (ExpMean>=0){
                    print geneID,geneName,ExpA,ExpB,ExpC,ExpMean,log(ExpMean)/log(2),log(ExpA)/log(2),log(ExpB)/log(2),log(ExpC)/log(2);
                }
            }
        }
    }' $annotation shWTAP_IP_rep1.genes.results shWTAP_IP_rep2.genes.results shWTAP_IP_rep3.genes.results >> Hela_RNA-seq_RSEM_shWTAP_IP.txt

