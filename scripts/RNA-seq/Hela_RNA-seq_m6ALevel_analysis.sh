#!/bin/sh

inputFile=/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/Hela/RSEM/Hela_RNA-seq_RSEM_shCont_input.txt
ipFile=/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/Hela/RSEM/Hela_RNA-seq_RSEM_shCont_IP.txt
output=/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/Hela/m6A_Level/

rm -rf $output
mkdir $output
cd $output

awk 'BEGIN{OFS="\t";}
    ARGIND==1{
        if(FNR > 1){
            if ($5 > 0) {
                hashGeneInput[$1]=$5;
            }
        }
    }
    ARGIND==2{
        if(FNR > 1){
            if ($5 > 0) {
                hashGeneIp[$1]=$5;
            }
        }
    }
    END{
        slen=asorti(hashGeneInput,hashGeneInputS);
        for (i=1;i<=slen;i++){
            geneID=hashGeneInputS[i];
            if(geneID in hashGeneIp){
                FC = hashGeneIp[geneID]/hashGeneInput[geneID]
                if (FC > 1) {
                    print geneID, FC;
                }
            }
        }
    }' $inputFile $ipFile > Hela_m6a_RSEM_shCont_FC.txt
