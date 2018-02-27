#!/bin/sh

### HepG2, RNA-seq
cd /data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/HepG2/RSEM/
targetFolder=/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/HepG2/RSEM-correlation

##### correlation log2(FPKM)

echo -e "GeneID\tshSETD2\tshCont" > $targetFolder/HepG2_RNA-seq_RSEM_shSetD2vsshCont_log2_input.txt
awk 'BEGIN{OFS="\t";}
    ARGIND==1{
        if(FNR>1){
            EXPa[$1]= $7;
        }
    }
    ARGIND==2{
        if(FNR>1){
            EXPb[$1]= $7;
        }
    }
    END{
        for (geneID in EXPa){
            if(geneID in EXPb){
                ExpA=EXPa[geneID];ExpB=EXPb[geneID];
                print geneID,ExpA,ExpB;
            }
        }
    }' HepG2_RNA-seq_RSEM_shSetD2_input.txt HepG2_RNA-seq_RSEM_shCont_input.txt >> $targetFolder/HepG2_RNA-seq_RSEM_shSetD2vsshCont_log2_input.txt

echo -e "GeneID\tshSETD2\tshMETTL14" > $targetFolder/HepG2_RNA-seq_RSEM_shSetD2vsshM14_log2_input.txt
awk 'BEGIN{OFS="\t";}
    ARGIND==1{
        if(FNR>1){
            EXPa[$1]= $7;
        }
    }
    ARGIND==2{
        if(FNR>1){
            EXPb[$1]= $7;
        }
    }
    END{
        for (geneID in EXPa){
            if(geneID in EXPb){
                ExpA=EXPa[geneID];ExpB=EXPb[geneID];
                print geneID,ExpA,ExpB;
            }
        }
    }' HepG2_RNA-seq_RSEM_shSetD2_input.txt HepG2_RNA-seq_RSEM_shM14_input.txt >> $targetFolder/HepG2_RNA-seq_RSEM_shSetD2vsshM14_log2_input.txt

echo -e "GeneID\tshSETD2\tshMETTL3" > $targetFolder/HepG2_RNA-seq_RSEM_shSetD2vsshM3_log2_input.txt
awk 'BEGIN{OFS="\t";}
    ARGIND==1{
        if(FNR>1){
            EXPa[$1]= $7;
        }
    }
    ARGIND==2{
        if(FNR>1){
            EXPb[$1]= $7;
        }
    }
    END{
        for (geneID in EXPa){
            if(geneID in EXPb){
                ExpA=EXPa[geneID];ExpB=EXPb[geneID];
                print geneID,ExpA,ExpB;
            }
        }
    }' HepG2_RNA-seq_RSEM_shSetD2_input.txt HepG2_RNA-seq_RSEM_shM3_input.txt >> $targetFolder/HepG2_RNA-seq_RSEM_shSetD2vsshM3_log2_input.txt

echo -e "GeneID\tshSETD2\tshWTAP" > $targetFolder/HepG2_RNA-seq_RSEM_shSetD2vsshWTAP_log2_input.txt
awk 'BEGIN{OFS="\t";}
    ARGIND==1{
        if(FNR>1){
            EXPa[$1]= $7;
        }
    }
    ARGIND==2{
        if(FNR>1){
            EXPb[$1]= $7;
        }
    }
    END{
        for (geneID in EXPa){
            if(geneID in EXPb){
                ExpA=EXPa[geneID];ExpB=EXPb[geneID];
                print geneID,ExpA,ExpB;
            }
        }
    }' HepG2_RNA-seq_RSEM_shSetD2_input.txt HepG2_RNA-seq_RSEM_shWTAP_input.txt >> $targetFolder/HepG2_RNA-seq_RSEM_shSetD2vsshWTAP_log2_input.txt

##### correlation log2(FPKM)

echo -e "GeneID\tshSETD2\tshCont" > $targetFolder/HepG2_RNA-seq_RSEM_shSetD2vsshCont_FPKM_input.txt
awk 'BEGIN{OFS="\t";}
    ARGIND==1{
        if(FNR>1){
            EXPa[$1]= $6;
        }
    }
    ARGIND==2{
        if(FNR>1){
            EXPb[$1]= $6;
        }
    }
    END{
        for (geneID in EXPa){
            if(geneID in EXPb){
                ExpA=EXPa[geneID];ExpB=EXPb[geneID];
                print geneID,ExpA,ExpB;
            }
        }
    }' HepG2_RNA-seq_RSEM_shSetD2_input.txt HepG2_RNA-seq_RSEM_shCont_input.txt >> $targetFolder/HepG2_RNA-seq_RSEM_shSetD2vsshCont_FPKM_input.txt

echo -e "GeneID\tshSETD2\tshMETTL14" > $targetFolder/HepG2_RNA-seq_RSEM_shSetD2vsshM14_FPKM_input.txt
awk 'BEGIN{OFS="\t";}
    ARGIND==1{
        if(FNR>1){
            EXPa[$1]= $6;
        }
    }
    ARGIND==2{
        if(FNR>1){
            EXPb[$1]= $6;
        }
    }
    END{
        for (geneID in EXPa){
            if(geneID in EXPb){
                ExpA=EXPa[geneID];ExpB=EXPb[geneID];
                print geneID,ExpA,ExpB;
            }
        }
    }' HepG2_RNA-seq_RSEM_shSetD2_input.txt HepG2_RNA-seq_RSEM_shM14_input.txt >> $targetFolder/HepG2_RNA-seq_RSEM_shSetD2vsshM14_FPKM_input.txt

echo -e "GeneID\tshSETD2\tshMETTL3" > $targetFolder/HepG2_RNA-seq_RSEM_shSetD2vsshM3_FPKM_input.txt
awk 'BEGIN{OFS="\t";}
    ARGIND==1{
        if(FNR>1){
            EXPa[$1]= $6;
        }
    }
    ARGIND==2{
        if(FNR>1){
            EXPb[$1]= $6;
        }
    }
    END{
        for (geneID in EXPa){
            if(geneID in EXPb){
                ExpA=EXPa[geneID];ExpB=EXPb[geneID];
                print geneID,ExpA,ExpB;
            }
        }
    }' HepG2_RNA-seq_RSEM_shSetD2_input.txt HepG2_RNA-seq_RSEM_shM3_input.txt >> $targetFolder/HepG2_RNA-seq_RSEM_shSetD2vsshM3_FPKM_input.txt

echo -e "GeneID\tshSETD2\tshWTAP" > $targetFolder/HepG2_RNA-seq_RSEM_shSetD2vsshWTAP_FPKM_input.txt
awk 'BEGIN{OFS="\t";}
    ARGIND==1{
        if(FNR>1){
            EXPa[$1]= $6;
        }
    }
    ARGIND==2{
        if(FNR>1){
            EXPb[$1]= $6;
        }
    }
    END{
        for (geneID in EXPa){
            if(geneID in EXPb){
                ExpA=EXPa[geneID];ExpB=EXPb[geneID];
                print geneID,ExpA,ExpB;
            }
        }
    }' HepG2_RNA-seq_RSEM_shSetD2_input.txt HepG2_RNA-seq_RSEM_shWTAP_input.txt >> $targetFolder/HepG2_RNA-seq_RSEM_shSetD2vsshWTAP_FPKM_input.txt

