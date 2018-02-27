#!/bin/sh

### HepG2, join-analysis
### gene expression in shCont: boxplot stats(0.00000 1.68182 3.05976 4.45516 8.58436), log2 transform(MeanFPKM>=1)
### ref gene m6A Level (log2) in shCont: boxplot stats (1.000148 1.477855 2.107638 2.966892 5.195494)
### ref gene H3k36me3 Level (log2) in shCont: boxplot stats (1.000006 1.358726 1.781892 2.311330 3.736260)
cd /data/zhoukr/hhl_setd2_m6a/analysis/joint-analysis/HepG2
annotation=/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.all.longest.exon.bed6
HistoneDir=/data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/htseq-count
m6ADir=/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/HepG2/IPvsInput
geneExpDir=/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/HepG2/RSEM

## all
echo -e "GeneID\tGeneName\tlog2FPKM\tlog2M6A" > HepG2_RNA-seq_gene_Exp_m6A_shCont_stats.txt
awk 'BEGIN{OFS="\t";
        cutExp=1;
        cutM6A=1;
        cutHistone=1;
    }
    ARGIND==1{
        split($4,geneIdI,":");
        geneIdi=geneIdI[1];
        geneNamei=geneIdI[2];
        hashGeneIdName[geneIdi]=geneNamei;
    }ARGIND==2{
        if(FNR>1){
            if($6>=cutExp){
                GeneExp[$1]=$7;
            }
        }
    }
    ARGIND==3{
        if(FNR>1){
            if($3>=cutM6A){
                if($2 in GeneExp){
                    print $2,hashGeneIdName[$2],GeneExp[$2],$3;
                }
            }
        }
    }' $annotation $geneExpDir/HepG2_RNA-seq_RSEM_shCont_input.txt $m6ADir/HepG2_m6A-seq_shCont_htseq-count_gene-edgeR.FC.txt >> HepG2_RNA-seq_gene_Exp_m6A_shCont_stats.txt

echo -e "GeneID\tGeneName\tlog2FPKM\tlog2H3K4me3" > HepG2_RNA-seq_gene_Exp_histone_shCont_stats.txt
awk 'BEGIN{OFS="\t";
        cutExp=1;
        cutM6A=1;
        cutHistone=1;
    }
    ARGIND==1{
        split($4,geneIdI,":");
        geneIdi=geneIdI[1];
        geneNamei=geneIdI[2];
        hashGeneIdName[geneIdi]=geneNamei;
    }ARGIND==2{
        if(FNR>1){
            if($6>=cutExp){
                GeneExp[$1]=$7;
            }
        }
    }
    ARGIND==3{
        if(FNR>1){
            if($3>=cutHistone){
                if($2 in GeneExp){
                    print $2,hashGeneIdName[$2],GeneExp[$2],$3;
                }
            }
        }
    }' $annotation $geneExpDir/HepG2_RNA-seq_RSEM_shCont_input.txt $HistoneDir/HepG2_chip-seq_shCont_htseq-count-edgeR.FC.txt >> HepG2_RNA-seq_gene_Exp_histone_shCont_stats.txt

echo -e "GeneID\tGeneName\tlog2M6A\tlog2H3K4me3" > HepG2_RNA-seq_gene_m6A_histone_shCont_stats.txt
awk 'BEGIN{OFS="\t";
        cutExp=1;
        cutM6A=1;
        cutHistone=1;
    }
    ARGIND==1{
        split($4,geneIdI,":");
        geneIdi=geneIdI[1];
        geneNamei=geneIdI[2];
        hashGeneIdName[geneIdi]=geneNamei;
    }ARGIND==2{
        if(FNR>1){
            if($3>=cutM6A){
                GeneExp[$2]=$3;
            }
        }
    }
    ARGIND==3{
        if(FNR>1){
            if($3>=cutHistone){
                if($2 in GeneExp){
                    print $2,hashGeneIdName[$2],GeneExp[$2],$3;
                }
            }
        }
    }' $annotation $m6ADir/HepG2_m6A-seq_shCont_htseq-count_gene-edgeR.FC.txt $HistoneDir/HepG2_chip-seq_shCont_htseq-count-edgeR.FC.txt >> HepG2_RNA-seq_gene_m6A_histone_shCont_stats.txt


### high
echo -e "GeneID\tGeneName\tlog2FPKM\tlog2M6A" > high/HepG2_RNA-seq_gene_Exp_m6A_shCont_stats.txt
awk 'BEGIN{OFS="\t";
        highExp=4.45516;lowExp=1.477855;
        highM6A=2.966892;lowM6A=1.477855;
        highHistone=2.311330;lowHistone=1.358726;
    }
    ARGIND==1{
        split($4,geneIdI,":");
        geneIdi=geneIdI[1];
        geneNamei=geneIdI[2];
        hashGeneIdName[geneIdi]=geneNamei;
    }ARGIND==2{
        if(FNR>1){
            if($6>=highExp){
                GeneExp[$1]=$7;
            }
        }
    }
    ARGIND==3{
        if(FNR>1){
            if($3>=highM6A){
                if($2 in GeneExp){
                    print $2,hashGeneIdName[$2],GeneExp[$2],$3;
                }
            }
        }
    }' $annotation $geneExpDir/HepG2_RNA-seq_RSEM_shCont_input.txt $m6ADir/HepG2_m6A-seq_shCont_htseq-count_gene-edgeR.FC.txt >> high/HepG2_RNA-seq_gene_Exp_m6A_shCont_stats.txt

echo -e "GeneID\tGeneName\tlog2FPKM\tlog2H3K4me3" > high/HepG2_RNA-seq_gene_Exp_histone_shCont_stats.txt
awk 'BEGIN{OFS="\t";
        highExp=4.45516;lowExp=1.477855;
        highM6A=2.966892;lowM6A=1.477855;
        highHistone=2.311330;lowHistone=1.358726;
    }
    ARGIND==1{
        split($4,geneIdI,":");
        geneIdi=geneIdI[1];
        geneNamei=geneIdI[2];
        hashGeneIdName[geneIdi]=geneNamei;
    }ARGIND==2{
        if(FNR>1){
            if($6>=highExp){
                GeneExp[$1]=$7;
            }
        }
    }
    ARGIND==3{
        if(FNR>1){
            if($3>=highHistone){
                if($2 in GeneExp){
                    print $2,hashGeneIdName[$2],GeneExp[$2],$3;
                }
            }
        }
    }' $annotation $geneExpDir/HepG2_RNA-seq_RSEM_shCont_input.txt $HistoneDir/HepG2_chip-seq_shCont_htseq-count-edgeR.FC.txt >> high/HepG2_RNA-seq_gene_Exp_histone_shCont_stats.txt

echo -e "GeneID\tGeneName\tlog2M6A\tlog2H3K4me3" > high/HepG2_RNA-seq_gene_m6A_histone_shCont_stats.txt
awk 'BEGIN{OFS="\t";
        highExp=4.45516;lowExp=1.477855;
        highM6A=2.966892;lowM6A=1.477855;
        highHistone=2.311330;lowHistone=1.358726;
    }
    ARGIND==1{
        split($4,geneIdI,":");
        geneIdi=geneIdI[1];
        geneNamei=geneIdI[2];
        hashGeneIdName[geneIdi]=geneNamei;
    }ARGIND==2{
        if(FNR>1){
            if($3>=highM6A){
                GeneExp[$2]=$3;
            }
        }
    }
    ARGIND==3{
        if(FNR>1){
            if($3>=highHistone){
                if($2 in GeneExp){
                    print $2,hashGeneIdName[$2],GeneExp[$2],$3;
                }
            }
        }
    }' $annotation $m6ADir/HepG2_m6A-seq_shCont_htseq-count_gene-edgeR.FC.txt $HistoneDir/HepG2_chip-seq_shCont_htseq-count-edgeR.FC.txt >> high/HepG2_RNA-seq_gene_m6A_histone_shCont_stats.txt

### modest
echo -e "GeneID\tGeneName\tlog2FPKM\tlog2M6A" > modest/HepG2_RNA-seq_gene_Exp_m6A_shCont_stats.txt
awk 'BEGIN{OFS="\t";
        highExp=4.45516;lowExp=1.477855;
        highM6A=2.966892;lowM6A=1.477855;
        highHistone=2.311330;lowHistone=1.358726;
    }
    ARGIND==1{
        split($4,geneIdI,":");
        geneIdi=geneIdI[1];
        geneNamei=geneIdI[2];
        hashGeneIdName[geneIdi]=geneNamei;
    }ARGIND==2{
        if(FNR>1){
            if(lowExp<$6<highExp){
                GeneExp[$1]=$7;
            }
        }
    }
    ARGIND==3{
        if(FNR>1){
            if(lowM6A<$3<highM6A){
                if($2 in GeneExp){
                    print $2,hashGeneIdName[$2],GeneExp[$2],$3;
                }
            }
        }
    }' $annotation $geneExpDir/HepG2_RNA-seq_RSEM_shCont_input.txt $m6ADir/HepG2_m6A-seq_shCont_htseq-count_gene-edgeR.FC.txt >> modest/HepG2_RNA-seq_gene_Exp_m6A_shCont_stats.txt

echo -e "GeneID\tGeneName\tlog2FPKM\tlog2H3K4me3" > modest/HepG2_RNA-seq_gene_Exp_histone_shCont_stats.txt
awk 'BEGIN{OFS="\t";
        highExp=4.45516;lowExp=1.477855;
        highM6A=2.966892;lowM6A=1.477855;
        highHistone=2.311330;lowHistone=1.358726;
    }
    ARGIND==1{
        split($4,geneIdI,":");
        geneIdi=geneIdI[1];
        geneNamei=geneIdI[2];
        hashGeneIdName[geneIdi]=geneNamei;
    }ARGIND==2{
        if(FNR>1){
            if(lowExp<$6<highExp){
                GeneExp[$1]=$7;
            }
        }
    }
    ARGIND==3{
        if(FNR>1){
            if(lowHistone<$3<highHistone){
                if($2 in GeneExp){
                    print $2,hashGeneIdName[$2],GeneExp[$2],$3;
                }
            }
        }
    }' $annotation $geneExpDir/HepG2_RNA-seq_RSEM_shCont_input.txt $HistoneDir/HepG2_chip-seq_shCont_htseq-count-edgeR.FC.txt >> modest/HepG2_RNA-seq_gene_Exp_histone_shCont_stats.txt

echo -e "GeneID\tGeneName\tlog2M6A\tlog2H3K4me3" > modest/HepG2_RNA-seq_gene_m6A_histone_shCont_stats.txt
awk 'BEGIN{OFS="\t";
        highExp=4.45516;lowExp=1.477855;
        highM6A=2.966892;lowM6A=1.477855;
        highHistone=2.311330;lowHistone=1.358726;
    }
    ARGIND==1{
        split($4,geneIdI,":");
        geneIdi=geneIdI[1];
        geneNamei=geneIdI[2];
        hashGeneIdName[geneIdi]=geneNamei;
    }ARGIND==2{
        if(FNR>1){
            if(lowM6A<$3<highM6A){
                GeneExp[$2]=$3;
            }
        }
    }
    ARGIND==3{
        if(FNR>1){
            if(lowHistone<$3<highHistone){
                if($2 in GeneExp){
                    print $2,hashGeneIdName[$2],GeneExp[$2],$3;
                }
            }
        }
    }' $annotation $m6ADir/HepG2_m6A-seq_shCont_htseq-count_gene-edgeR.FC.txt $HistoneDir/HepG2_chip-seq_shCont_htseq-count-edgeR.FC.txt >> modest/HepG2_RNA-seq_gene_m6A_histone_shCont_stats.txt


##low
echo -e "GeneID\tGeneName\tlog2FPKM\tlog2M6A" > low/HepG2_RNA-seq_gene_Exp_m6A_shCont_stats.txt
awk 'BEGIN{OFS="\t";
        highExp=4.45516;lowExp=1.477855;
        highM6A=2.966892;lowM6A=1.477855;
        highHistone=2.311330;lowHistone=1.358726;
    }
    ARGIND==1{
        split($4,geneIdI,":");
        geneIdi=geneIdI[1];
        geneNamei=geneIdI[2];
        hashGeneIdName[geneIdi]=geneNamei;
    }ARGIND==2{
        if(FNR>1){
            if($6<=lowExp){
                GeneExp[$1]=$7;
            }
        }
    }
    ARGIND==3{
        if(FNR>1){
            if($3<=lowM6A){
                if($2 in GeneExp){
                    print $2,hashGeneIdName[$2],GeneExp[$2],$3;
                }
            }
        }
    }' $annotation $geneExpDir/HepG2_RNA-seq_RSEM_shCont_input.txt $m6ADir/HepG2_m6A-seq_shCont_htseq-count_gene-edgeR.FC.txt >> low/HepG2_RNA-seq_gene_Exp_m6A_shCont_stats.txt

echo -e "GeneID\tGeneName\tlog2FPKM\tlog2H3K4me3" > low/HepG2_RNA-seq_gene_Exp_histone_shCont_stats.txt
awk 'BEGIN{OFS="\t";
        highExp=4.45516;lowExp=1.477855;
        highM6A=2.966892;lowM6A=1.477855;
        highHistone=2.311330;lowHistone=1.358726;
    }
    ARGIND==1{
        split($4,geneIdI,":");
        geneIdi=geneIdI[1];
        geneNamei=geneIdI[2];
        hashGeneIdName[geneIdi]=geneNamei;
    }ARGIND==2{
        if(FNR>1){
            if($6<=lowExp){
                GeneExp[$1]=$7;
            }
        }
    }
    ARGIND==3{
        if(FNR>1){
            if($3<=lowHistone){
                if($2 in GeneExp){
                    print $2,hashGeneIdName[$2],GeneExp[$2],$3;
                }
            }
        }
    }' $annotation $geneExpDir/HepG2_RNA-seq_RSEM_shCont_input.txt $HistoneDir/HepG2_chip-seq_shCont_htseq-count-edgeR.FC.txt >> low/HepG2_RNA-seq_gene_Exp_histone_shCont_stats.txt

echo -e "GeneID\tGeneName\tlog2M6A\tlog2H3K4me3" > low/HepG2_RNA-seq_gene_m6A_histone_shCont_stats.txt
awk 'BEGIN{OFS="\t";
        highExp=4.45516;lowExp=1.477855;
        highM6A=2.966892;lowM6A=1.477855;
        highHistone=2.311330;lowHistone=1.358726;
    }
    ARGIND==1{
        split($4,geneIdI,":");
        geneIdi=geneIdI[1];
        geneNamei=geneIdI[2];
        hashGeneIdName[geneIdi]=geneNamei;
    }ARGIND==2{
        if(FNR>1){
            if($3<=lowM6A){
                GeneExp[$2]=$3;
            }
        }
    }
    ARGIND==3{
        if(FNR>1){
            if($3<=lowHistone){
                if($2 in GeneExp){
                    print $2,hashGeneIdName[$2],GeneExp[$2],$3;
                }
            }
        }
    }' $annotation $m6ADir/HepG2_m6A-seq_shCont_htseq-count_gene-edgeR.FC.txt $HistoneDir/HepG2_chip-seq_shCont_htseq-count-edgeR.FC.txt >> low/HepG2_RNA-seq_gene_m6A_histone_shCont_stats.txt
