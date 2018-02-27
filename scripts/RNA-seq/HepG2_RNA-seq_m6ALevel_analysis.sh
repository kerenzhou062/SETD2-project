#!/bin/sh

inputFile=/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/HepG2/RSEM/HepG2_RNA-seq_RSEM_shCont_input.txt
ipFile=/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/HepG2/RSEM/HepG2_RNA-seq_RSEM_shCont_IP.txt
output=/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/HepG2/m6A_Level/

rm -rf $output
mkdir $output
cd $output

awk 'BEGIN{OFS="\t";}
    ARGIND==1{
        if(FNR > 1){
            if ($6 > 0) {
                hashGeneInput[$1]=$6;
            }
        }
    }
    ARGIND==2{
        if(FNR > 1){
            if ($6 > 0) {
                hashGeneIp[$1]=$6;
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
    }' $inputFile $ipFile > HepG2_m6a_RSEM_shCont_FC.txt
