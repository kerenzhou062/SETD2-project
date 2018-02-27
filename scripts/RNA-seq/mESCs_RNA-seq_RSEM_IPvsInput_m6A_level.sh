#!/bin/sh

### Hela, ChIP-seq
cd /data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/mESCs/IPvsInput
annotation=/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.vM9.annotation.all.longest.exon.bed6
geneExpDir=/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/mESCs/RSEM

### Gene M6A level assessment
echo -e "GeneID\tgeneName-lowercase\tgeneName-uppercase\tInput(FPKM)\tIP(FPKM)\tGeneM6ALevel\tGeneM6ALevel(log2((IP+0.01)/(Input+0.01)))" > mESCs_gene_m6A_level_shCont-D0.txt
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
    }' $annotation $geneExpDir/Ctrl_D0_input.genes.results $geneExpDir/Ctrl_D0_IP.genes.results >> mESCs_gene_m6A_level_shCont-D0.txt

echo -e "GeneID\tgeneName-lowercase\tgeneName-uppercase\tInput(FPKM)\tIP(FPKM)\tGeneM6ALevel\tGeneM6ALevel(log2((IP+0.01)/(Input+0.01)))" > mESCs_gene_m6A_level_shSetD2-D0.txt
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
    }' $annotation $geneExpDir/SetD2-KD_D0_input.genes.results $geneExpDir/SetD2-KD_D0_IP.genes.results >> mESCs_gene_m6A_level_shSetD2-D0.txt

echo -e "GeneID\tgeneName-lowercase\tgeneName-uppercase\tInput(FPKM)\tIP(FPKM)\tGeneM6ALevel\tGeneM6ALevel(log2((IP+0.01)/(Input+0.01)))" > mESCs_gene_m6A_level_shCont-D6.txt
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
    }' $annotation $geneExpDir/Ctrl_D6_input.genes.results $geneExpDir/Ctrl_D6_IP.genes.results >> mESCs_gene_m6A_level_shCont-D6.txt

echo -e "GeneID\tgeneName-lowercase\tgeneName-uppercase\tInput(FPKM)\tIP(FPKM)\tGeneM6ALevel\tGeneM6ALevel(log2((IP+0.01)/(Input+0.01)))" > mESCs_gene_m6A_level_shSetD2-D6.txt
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
    }' $annotation $geneExpDir/SetD2-KD_D6_input.genes.results $geneExpDir/SetD2-KD_D6_IP.genes.results >> mESCs_gene_m6A_level_shSetD2-D6.txt

### For GSEA analysis

echo -e "4\t4\t1\n#shContD0\tshContD6\tshSetD2D0\tshSetD2D6\nshContD0\tshContD6\tshSetD2D0\tshSetD2D6" > mESCs_m6a_Phenotype.cls

echo -e "NAME\tDESCRIPTION\tshContD0\tshContD6\tshSetD2D0\tshSetD2D6" > mESCs_gene_m6A_level_GSEA_geneUppercase.txt
awk 'BEGIN{OFS="\t";}
    ARGIND==1{
        if(FNR>1){
            if($4>=1 && $7>=1){
                if($3 in LevelA){
                    if($7>LevelA[$3]){
                        LevelA[$3]=$7;
                    }
                }else{
                    LevelA[$3]= $7;
                }
            }
        }
    }
    ARGIND==2{
        if(FNR>1){
            if($3 in LevelB){
                if($7>LevelB[$3]){
                    LevelB[$3]= $7;
                }
            }else{
                LevelB[$3]= $7;
            }
        }
    }
    ARGIND==3{
        if(FNR>1){
            if($3 in LevelC){
                if($7>LevelC[$3]){
                    LevelC[$3]= $7;
                }
            }else{
                LevelC[$3]= $7;
            }
        }
    }
    ARGIND==4{
        if(FNR>1){
            if($3 in LevelD){
                if($7>LevelD[$3]){
                    LevelD[$3]= $7;
                }
            }else{
                LevelD[$3]= $7;
            }
        }
    }
    END{
        for(geneName in LevelA){
            Levela=LevelA[geneName];Levelb=LevelB[geneName];Levelc=LevelC[geneName];Leveld=LevelD[geneName];
            print geneName,"na",Levela,Levelb,Levelc,Leveld;
        }
    }' mESCs_gene_m6A_level_shCont-D0.txt mESCs_gene_m6A_level_shCont-D6.txt \
       mESCs_gene_m6A_level_shSetD2-D0.txt mESCs_gene_m6A_level_shSetD2-D6.txt >> mESCs_gene_m6A_level_GSEA_geneUppercase.txt

echo -e "NAME\tDESCRIPTION\tshContD0\tshContD6\tshSetD2D0\tshSetD2D6" > mESCs_gene_m6A_level_GSEA_geneLowercase.txt
awk 'BEGIN{OFS="\t";}
    ARGIND==1{
        if(FNR>1){
            if($4>=1 && $7>=1){
                if($2 in LevelA){
                    if($7>LevelA[$2]){
                        LevelA[$2]=$7;
                    }
                }else{
                    LevelA[$2]= $7;
                }
            }
        }
    }
    ARGIND==2{
        if(FNR>1){
            if($2 in LevelB){
                if($7>LevelB[$2]){
                    LevelB[$2]= $7;
                }
            }else{
                LevelB[$2]= $7;
            }
        }
    }
    ARGIND==3{
        if(FNR>1){
            if($2 in LevelC){
                if($7>LevelC[$2]){
                    LevelC[$2]= $7;
                }
            }else{
                LevelC[$2]= $7;
            }
        }
    }
    ARGIND==4{
        if(FNR>1){
            if($2 in LevelD){
                if($7>LevelD[$2]){
                    LevelD[$2]= $7;
                }
            }else{
                LevelD[$2]= $7;
            }
        }
    }
    END{
        for(geneName in LevelA){
            Levela=LevelA[geneName];Levelb=LevelB[geneName];Levelc=LevelC[geneName];Leveld=LevelD[geneName];
            print geneName,"na",Levela,Levelb,Levelc,Leveld;
        }
    }' mESCs_gene_m6A_level_shCont-D0.txt mESCs_gene_m6A_level_shCont-D6.txt \
       mESCs_gene_m6A_level_shSetD2-D0.txt mESCs_gene_m6A_level_shSetD2-D6.txt >> mESCs_gene_m6A_level_GSEA_geneLowercase.txt
