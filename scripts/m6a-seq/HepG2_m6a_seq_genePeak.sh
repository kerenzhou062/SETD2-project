#!/bin/sh
xlsFolder=/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/xls
output=/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/xls/genePeak

rm -rf $output
mkdir $output
cd $output

awk 'BEGIN{OFS="\t";FS="\t";}
    {
        geneID=$4;folEnri=$15;
        if(geneID in hashGene){
            if(folEnri > hashGene[geneID]){
                hashGene[geneID] = folEnri;
            }
        }else{
            hashGene[geneID] = folEnri;
        }
    }
    END{
        for (i in hashGene) {
            print i, hashGene[i];
        }
    }
' $xlsFolder/HepG2_shCont_m6a.xls | sort -k 1,1 > $output/HepG2_shCont_gene_m6a.txt
