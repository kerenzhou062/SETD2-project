#!/bin/sh

### Hela, ChIP-seq
cd /data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/HepG2/m6A_Exp_stats
m6ADir=/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/HepG2/IPvsInput
geneExpDir=/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/HepG2/RSEM
### gene expression in shCont: boxplot stats(0.00000 1.68182 3.05976 4.45516 8.58436), log2 transform(MeanFPKM>=1)
### highExpCut=4.45516,lowExpCut=1.68182
### gene m6A FC Determined by Huilin Huang: non-targets = lower than 1.2 fold reduction, hypo-low = 1.2-1.5 fold reduction, hypo-modest = 1.5-2 fold reduction, hypo-high = more than 2 fold reduction

echo -e "## gene expression and m6a FC statistics in HepG2" > HepG2_RNA-seq_gene_m6A_Exp_shSetD2_stats.txt
awk 'BEGIN{OFS="\t";
        m6AHypoLow=-log(1.2)/log(2);m6AHypoModest=-log(1.5)/log(2);m6AHypoHigh=-1;
        m6AHyperLow=log(1.2)/log(2);m6AHyperModest=log(1.5)/log(2);m6AHyperHigh=1;
        highExpCut=4.45516;lowExpCut=1.68182;
        dowmM6ACount=0;upM6ACount=0;totalM6ACount=0;
    }
    ARGIND==1{
        if(FNR>1){
            if($3<=m6AHypoHigh){
                hashHypoHighM6A[$1]=$2"\t"$3"\t"$4;
            }else if(m6AHypoHigh<$3 && $3<m6AHypoModest){
                hashHypoModestM6A[$1]=$2"\t"$3"\t"$4;
            }else if(m6AHypoModest<=$3 && $3<=m6AHypoLow){
                hashHypoLowM6A[$1]=$2"\t"$3"\t"$4;
            }else if(m6AHyperLow<=$3 && $3<=m6AHyperModest){
                hashHyperLowM6A[$1]=$2"\t"$3"\t"$4;
            }else if(m6AHyperModest<$3 && $3<m6AHyperHigh){
                hashHyperModestM6A[$1]=$2"\t"$3"\t"$4;
            }else if($3>=m6AHyperHigh){
                hashHyperHighM6A[$1]=$2"\t"$3"\t"$4;
            }
        }
        totalM6ACount++;
        if($3<0){
            dowmM6ACount++;
        }else if($3>=0){
            upM6ACount++;
        }
    }
    ARGIND==2{
        if(FNR>1){
            hashGeneExp[$1]=$6;
            if($7<=lowExpCut){
                hashLowExp[$1]=$6;
            }else if(lowExpCut<$7 && $7<highExpCut){
                hashModestExp[$1]=$6;
            }else if($7>=highExpCut){
                hashHighExp[$1]=$6;
            }
        }
    }
    END{
        print "##gene m6A FC in shSetD2 vs. shCont:";
        print "##Determined by Huilin: non-targets = lower than 1.2 fold reduction, hypo-low = 1.2-1.5 fold reduction, hypo-modest = 1.5-2 fold reduction, hypo-high = more than 2 fold reduction";
        print "##ref gene m6A Level (log2) in shCont: boxplot stats (1.000148 1.481705 2.108830 2.962674 5.163127)";
        print "HypoHigh\tHypoModest\tHypoLow\tHypoTotal\tHyperLow\tHyperModest\tHyperHigh\tHyperTotal";
        print length(hashHypoHighM6A),length(hashHypoModestM6A),length(hashHypoLowM6A),dowmM6ACount,length(hashHyperLowM6A),length(hashHyperModestM6A),length(hashHyperHighM6A),upM6ACount;
        print "\n##gene Expression(FPKM) in shCont:";
        print "##Determined by boxplot stats(0.00000 1.68182 3.05976 4.45516 8.58436), log2 transform(MeanFPKM>=1)";
        print "##highExpCut=4.709380,lowExpCut=2.194300";
        print "HighExp\tModestExp\tLowExp";
        print length(hashHighExp),length(hashModestExp),length(hashLowExp),"\n";

        HyperHighM6ACount=0;
        HyperModestM6ACount=0;
        HyperLowM6ACount=0;
        HypoHighM6ACount=0;
        HypoModestM6ACount=0;
        HypoLowM6ACount=0;
        for(gene in hashHighExp){
            if(gene in hashHyperHighM6A){
                HyperHighM6ACount++;
            }
            if(gene in hashHyperModestM6A){
                HyperModestM6ACount++;
            }
            if(gene in hashHyperLowM6A){
                HyperLowM6ACount++;
            }
            if(gene in hashHypoHighM6A){
                HypoHighM6ACount++;
                hashHypoHighM6AInExpHigh[gene]=1;
            }
            if(gene in hashHypoModestM6A){
                HypoModestM6ACount++;
                hashHypoModestM6AInExpHigh[gene]=1;
            }
            if(gene in hashHypoLowM6A){
                HypoLowM6ACount++;
                hashHypoLowM6AInExpHigh[gene]=1;
            }
        }
        print "##high Expresion Group:";
        print "Type\tCount\tPercentage";
        print "HypoLowM6ACount",HypoLowM6ACount,HypoLowM6ACount/totalM6ACount*100;
        print "HypoModestM6ACount",HypoModestM6ACount,HypoModestM6ACount/totalM6ACount*100;
        print "HypoHighM6ACount",HypoHighM6ACount,HypoHighM6ACount/totalM6ACount*100;
        print "HyperLowM6ACount",HyperLowM6ACount,HyperLowM6ACount/totalM6ACount*100;
        print "HyperModestM6ACount",HyperModestM6ACount,HyperModestM6ACount/totalM6ACount*100;
        print "HyperHighM6ACount",HyperHighM6ACount,HyperHighM6ACount/totalM6ACount*100, "\n";

        HyperHighM6ACount=0;
        HyperModestM6ACount=0;
        HyperLowM6ACount=0;
        HypoHighM6ACount=0;
        HypoModestM6ACount=0;
        HypoLowM6ACount=0;
        for(gene in hashModestExp){
            if(gene in hashHyperHighM6A){
                HyperHighM6ACount++;
            }
            if(gene in hashHyperModestM6A){
                HyperModestM6ACount++;
            }
            if(gene in hashHyperLowM6A){
                HyperLowM6ACount++;
            }
            if(gene in hashHypoHighM6A){
                HypoHighM6ACount++;
            }
            if(gene in hashHypoModestM6A){
                HypoModestM6ACount++;
            }
            if(gene in hashHypoLowM6A){
                HypoLowM6ACount++;
            }
        }
        print "##modest Expresion Group:";
        print "Type\tCount\tPercentage";
        print "HypoLowM6ACount",HypoLowM6ACount,HypoLowM6ACount/totalM6ACount*100;
        print "HypoModestM6ACount",HypoModestM6ACount,HypoModestM6ACount/totalM6ACount*100;
        print "HypoHighM6ACount",HypoHighM6ACount,HypoHighM6ACount/totalM6ACount*100;
        print "HyperLowM6ACount",HyperLowM6ACount,HyperLowM6ACount/totalM6ACount*100;
        print "HyperModestM6ACount",HyperModestM6ACount,HyperModestM6ACount/totalM6ACount*100;
        print "HyperHighM6ACount",HyperHighM6ACount,HyperHighM6ACount/totalM6ACount*100, "\n";

        HyperHighM6ACount=0;
        HyperModestM6ACount=0;
        HyperLowM6ACount=0;
        HypoHighM6ACount=0;
        HypoModestM6ACount=0;
        HypoLowM6ACount=0;
        for(gene in hashLowExp){
            if(gene in hashHyperHighM6A){
                HyperHighM6ACount++;
            }
            if(gene in hashHyperModestM6A){
                HyperModestM6ACount++;
            }
            if(gene in hashHyperLowM6A){
                HyperLowM6ACount++;
            }
            if(gene in hashHypoHighM6A){
                HypoHighM6ACount++;
            }
            if(gene in hashHypoModestM6A){
                HypoModestM6ACount++;
            }
            if(gene in hashHypoLowM6A){
                HypoLowM6ACount++;
            }
        }
        print "##low Expresion Group:";
        print "Type\tCount\tPercentage";
        print "HypoLowM6ACount",HypoLowM6ACount,HypoLowM6ACount/totalM6ACount*100;
        print "HypoModestM6ACount",HypoModestM6ACount,HypoModestM6ACount/totalM6ACount*100;
        print "HypoHighM6ACount",HypoHighM6ACount,HypoHighM6ACount/totalM6ACount*100;
        print "HyperLowM6ACount",HyperLowM6ACount,HyperLowM6ACount/totalM6ACount*100;
        print "HyperModestM6ACount",HyperModestM6ACount,HyperModestM6ACount/totalM6ACount*100;
        print "HyperHighM6ACount",HyperHighM6ACount,HyperHighM6ACount/totalM6ACount*100, "\n";
        print "";
        print "High hypo m6A in high gene expression:";
        print "GeneID\tGeneName\tGeneM6AFC(log2)\tGeneM6ALevelInshCont(log2)\tGeneMeanFPKM";
        for(gene in hashHypoHighM6AInExpHigh){
            print gene,hashHypoHighM6A[gene],hashGeneExp[gene];
        }
        print "";
        print "Modest hypo m6A in high gene expression:";
        print "GeneID\tGeneName\tGeneM6AFC(log2)\tGeneM6ALevelInshCont(log2)\tGeneMeanFPKM";
        for(gene in hashHypoModestM6AInExpHigh){
            print gene,hashHypoModestM6A[gene],hashGeneExp[gene];
        }
        print "";
        print "Low hypo m6A in high gene expression:";
        print "GeneID\tGeneName\tGeneM6AFC(log2)\tGeneM6ALevelInshCont(log2)\tGeneMeanFPKM";
        for(gene in hashHypoLowM6AInExpHigh){
            print gene,hashHypoLowM6A[gene],hashGeneExp[gene];
        }
    }' $m6ADir/HepG2_RNA-seq_gene_m6A_level_FC_shSetD2vsshCont.txt $geneExpDir/HepG2_RNA-seq_RSEM_shCont_input.txt >> HepG2_RNA-seq_gene_m6A_Exp_shSetD2_stats.txt

echo -e "## gene expression and m6a FC statistics in HepG2" > HepG2_RNA-seq_gene_m6A_Exp_shM14_stats.txt
awk 'BEGIN{OFS="\t";
        m6AHypoLow=-log(1.2)/log(2);m6AHypoModest=-log(1.5)/log(2);m6AHypoHigh=-1;
        m6AHyperLow=log(1.2)/log(2);m6AHyperModest=log(1.5)/log(2);m6AHyperHigh=1;
        highExpCut=4.45516;lowExpCut=1.68182;
        dowmM6ACount=0;upM6ACount=0;totalM6ACount=0;
    }
    ARGIND==1{
        if(FNR>1){
            if($3<=m6AHypoHigh){
                hashHypoHighM6A[$1]=$2"\t"$3"\t"$4;
            }else if(m6AHypoHigh<$3 && $3<m6AHypoModest){
                hashHypoModestM6A[$1]=$2"\t"$3"\t"$4;
            }else if(m6AHypoModest<=$3 && $3<=m6AHypoLow){
                hashHypoLowM6A[$1]=$2"\t"$3"\t"$4;
            }else if(m6AHyperLow<=$3 && $3<=m6AHyperModest){
                hashHyperLowM6A[$1]=$2"\t"$3"\t"$4;
            }else if(m6AHyperModest<$3 && $3<m6AHyperHigh){
                hashHyperModestM6A[$1]=$2"\t"$3"\t"$4;
            }else if($3>=m6AHyperHigh){
                hashHyperHighM6A[$1]=$2"\t"$3"\t"$4;
            }
        }
        totalM6ACount++;
        if($3<0){
            dowmM6ACount++;
        }else if($3>=0){
            upM6ACount++
        }
    }
    ARGIND==2{
        if(FNR>1){
            hashGeneExp[$1]=$6;
            if($7<=lowExpCut){
                hashLowExp[$1]=$6;
            }else if(lowExpCut<$7 && $7<highExpCut){
                hashModestExp[$1]=$6;
            }else if($7>=highExpCut){
                hashHighExp[$1]=$6;
            }
        }
    }
    END{
        print "##gene m6A FC in shM14 vs. shCont:";
        print "##Determined by Huilin: non-targets = lower than 1.2 fold reduction, hypo-low = 1.2-1.5 fold reduction, hypo-modest = 1.5-2 fold reduction, hypo-high = more than 2 fold reduction";
        print "##ref gene m6A Level (log2) in shCont: boxplot stats (1.000148 1.481705 2.108830 2.962674 5.163127)";
        print "HypoHigh\tHypoModest\tHypoLow\tHypoTotal\tHyperLow\tHyperModest\tHyperHigh\tHyperTotal";
        print length(hashHypoHighM6A),length(hashHypoModestM6A),length(hashHypoLowM6A),dowmM6ACount,length(hashHyperLowM6A),length(hashHyperModestM6A),length(hashHyperHighM6A),upM6ACount;
        print "\n##gene Expression(FPKM) in shCont:";
        print "##Determined by boxplot stats(0.00000 1.68182 3.05976 4.45516 8.58436), log2 transform(MeanFPKM>=1)";
        print "##highExpCut=4.709380,lowExpCut=2.194300";
        print "HighExp\tModestExp\tLowExp";
        print length(hashHighExp),length(hashModestExp),length(hashLowExp),"\n";

        HyperHighM6ACount=0;
        HyperModestM6ACount=0;
        HyperLowM6ACount=0;
        HypoHighM6ACount=0;
        HypoModestM6ACount=0;
        HypoLowM6ACount=0;
        for(gene in hashHighExp){
            if(gene in hashHyperHighM6A){
                HyperHighM6ACount++;
            }
            if(gene in hashHyperModestM6A){
                HyperModestM6ACount++;
            }
            if(gene in hashHyperLowM6A){
                HyperLowM6ACount++;
            }
            if(gene in hashHypoHighM6A){
                HypoHighM6ACount++;
                hashHypoHighM6AInExpHigh[gene]=1;
            }
            if(gene in hashHypoModestM6A){
                HypoModestM6ACount++;
                hashHypoModestM6AInExpHigh[gene]=1;
            }
            if(gene in hashHypoLowM6A){
                HypoLowM6ACount++;
                hashHypoLowM6AInExpHigh[gene]=1;
            }
        }
        print "##high Expresion Group:";
        print "Type\tCount\tPercentage";
        print "HypoLowM6ACount",HypoLowM6ACount,HypoLowM6ACount/totalM6ACount*100;
        print "HypoModestM6ACount",HypoModestM6ACount,HypoModestM6ACount/totalM6ACount*100;
        print "HypoHighM6ACount",HypoHighM6ACount,HypoHighM6ACount/totalM6ACount*100;
        print "HyperLowM6ACount",HyperLowM6ACount,HyperLowM6ACount/totalM6ACount*100;
        print "HyperModestM6ACount",HyperModestM6ACount,HyperModestM6ACount/totalM6ACount*100;
        print "HyperHighM6ACount",HyperHighM6ACount,HyperHighM6ACount/totalM6ACount*100, "\n";

        HyperHighM6ACount=0;
        HyperModestM6ACount=0;
        HyperLowM6ACount=0;
        HypoHighM6ACount=0;
        HypoModestM6ACount=0;
        HypoLowM6ACount=0;
        for(gene in hashModestExp){
            if(gene in hashHyperHighM6A){
                HyperHighM6ACount++;
            }
            if(gene in hashHyperModestM6A){
                HyperModestM6ACount++;
            }
            if(gene in hashHyperLowM6A){
                HyperLowM6ACount++;
            }
            if(gene in hashHypoHighM6A){
                HypoHighM6ACount++;
            }
            if(gene in hashHypoModestM6A){
                HypoModestM6ACount++;
            }
            if(gene in hashHypoLowM6A){
                HypoLowM6ACount++;
            }
        }
        print "##modest Expresion Group:";
        print "Type\tCount\tPercentage";
        print "HypoLowM6ACount",HypoLowM6ACount,HypoLowM6ACount/totalM6ACount*100;
        print "HypoModestM6ACount",HypoModestM6ACount,HypoModestM6ACount/totalM6ACount*100;
        print "HypoHighM6ACount",HypoHighM6ACount,HypoHighM6ACount/totalM6ACount*100;
        print "HyperLowM6ACount",HyperLowM6ACount,HyperLowM6ACount/totalM6ACount*100;
        print "HyperModestM6ACount",HyperModestM6ACount,HyperModestM6ACount/totalM6ACount*100;
        print "HyperHighM6ACount",HyperHighM6ACount,HyperHighM6ACount/totalM6ACount*100, "\n";

        HyperHighM6ACount=0;
        HyperModestM6ACount=0;
        HyperLowM6ACount=0;
        HypoHighM6ACount=0;
        HypoModestM6ACount=0;
        HypoLowM6ACount=0;
        for(gene in hashLowExp){
            if(gene in hashHyperHighM6A){
                HyperHighM6ACount++;
            }
            if(gene in hashHyperModestM6A){
                HyperModestM6ACount++;
            }
            if(gene in hashHyperLowM6A){
                HyperLowM6ACount++;
            }
            if(gene in hashHypoHighM6A){
                HypoHighM6ACount++;
            }
            if(gene in hashHypoModestM6A){
                HypoModestM6ACount++;
            }
            if(gene in hashHypoLowM6A){
                HypoLowM6ACount++;
            }
        }
        print "##low Expresion Group:";
        print "Type\tCount\tPercentage";
        print "HypoLowM6ACount",HypoLowM6ACount,HypoLowM6ACount/totalM6ACount*100;
        print "HypoModestM6ACount",HypoModestM6ACount,HypoModestM6ACount/totalM6ACount*100;
        print "HypoHighM6ACount",HypoHighM6ACount,HypoHighM6ACount/totalM6ACount*100;
        print "HyperLowM6ACount",HyperLowM6ACount,HyperLowM6ACount/totalM6ACount*100;
        print "HyperModestM6ACount",HyperModestM6ACount,HyperModestM6ACount/totalM6ACount*100;
        print "HyperHighM6ACount",HyperHighM6ACount,HyperHighM6ACount/totalM6ACount*100, "\n";
        print "";
        print "High hypo m6A in high gene expression:";
        print "GeneID\tGeneName\tGeneM6AFC(log2)\tGeneM6ALevelInshCont(log2)\tGeneMeanFPKM";
        for(gene in hashHypoHighM6AInExpHigh){
            print gene,hashHypoHighM6A[gene],hashGeneExp[gene];
        }
        print "";
        print "Modest hypo m6A in high gene expression:";
        print "GeneID\tGeneName\tGeneM6AFC(log2)\tGeneM6ALevelInshCont(log2)\tGeneMeanFPKM";
        for(gene in hashHypoModestM6AInExpHigh){
            print gene,hashHypoModestM6A[gene],hashGeneExp[gene];
        }
        print "";
        print "Low hypo m6A in high gene expression:";
        print "GeneID\tGeneName\tGeneM6AFC(log2)\tGeneM6ALevelInshCont(log2)\tGeneMeanFPKM";
        for(gene in hashHypoLowM6AInExpHigh){
            print gene,hashHypoLowM6A[gene],hashGeneExp[gene];
        }
    }' $m6ADir/HepG2_RNA-seq_gene_m6A_level_FC_shM14vsshCont.txt $geneExpDir/HepG2_RNA-seq_RSEM_shCont_input.txt >> HepG2_RNA-seq_gene_m6A_Exp_shM14_stats.txt

echo -e "## gene expression and m6a FC statistics in HepG2" > HepG2_RNA-seq_gene_m6A_Exp_shM3_stats.txt
awk 'BEGIN{OFS="\t";
        m6AHypoLow=-log(1.2)/log(2);m6AHypoModest=-log(1.5)/log(2);m6AHypoHigh=-1;
        m6AHyperLow=log(1.2)/log(2);m6AHyperModest=log(1.5)/log(2);m6AHyperHigh=1;
        highExpCut=4.45516;lowExpCut=1.68182;
        dowmM6ACount=0;upM6ACount=0;totalM6ACount=0;
    }
    ARGIND==1{
        if(FNR>1){
            if($3<=m6AHypoHigh){
                hashHypoHighM6A[$1]=$2"\t"$3"\t"$4;
            }else if(m6AHypoHigh<$3 && $3<m6AHypoModest){
                hashHypoModestM6A[$1]=$2"\t"$3"\t"$4;
            }else if(m6AHypoModest<=$3 && $3<=m6AHypoLow){
                hashHypoLowM6A[$1]=$2"\t"$3"\t"$4;
            }else if(m6AHyperLow<=$3 && $3<=m6AHyperModest){
                hashHyperLowM6A[$1]=$2"\t"$3"\t"$4;
            }else if(m6AHyperModest<$3 && $3<m6AHyperHigh){
                hashHyperModestM6A[$1]=$2"\t"$3"\t"$4;
            }else if($3>=m6AHyperHigh){
                hashHyperHighM6A[$1]=$2"\t"$3"\t"$4;
            }
        }
        totalM6ACount++;
        if($3<0){
            dowmM6ACount++;
        }else if($3>=0){
            upM6ACount++
        }
    }
    ARGIND==2{
        if(FNR>1){
            hashGeneExp[$1]=$6;
            if($7<=lowExpCut){
                hashLowExp[$1]=$6;
            }else if(lowExpCut<$7 && $7<highExpCut){
                hashModestExp[$1]=$6;
            }else if($7>=highExpCut){
                hashHighExp[$1]=$6;
            }
        }
    }
    END{
        print "##gene m6A FC in shM3 vs. shCont:";
        print "##Determined by Huilin: non-targets = lower than 1.2 fold reduction, hypo-low = 1.2-1.5 fold reduction, hypo-modest = 1.5-2 fold reduction, hypo-high = more than 2 fold reduction";
        print "##ref gene m6A Level (log2) in shCont: boxplot stats (1.000148 1.481705 2.108830 2.962674 5.163127)";
        print "HypoHigh\tHypoModest\tHypoLow\tHypoTotal\tHyperLow\tHyperModest\tHyperHigh\tHyperTotal";
        print length(hashHypoHighM6A),length(hashHypoModestM6A),length(hashHypoLowM6A),dowmM6ACount,length(hashHyperLowM6A),length(hashHyperModestM6A),length(hashHyperHighM6A),upM6ACount;
        print "\n##gene Expression(FPKM) in shCont:";
        print "##Determined by boxplot stats(0.00000 1.68182 3.05976 4.45516 8.58436), log2 transform(MeanFPKM>=1)";
        print "##highExpCut=4.709380,lowExpCut=2.194300";
        print "HighExp\tModestExp\tLowExp";
        print length(hashHighExp),length(hashModestExp),length(hashLowExp),"\n";

        HyperHighM6ACount=0;
        HyperModestM6ACount=0;
        HyperLowM6ACount=0;
        HypoHighM6ACount=0;
        HypoModestM6ACount=0;
        HypoLowM6ACount=0;
        for(gene in hashHighExp){
            if(gene in hashHyperHighM6A){
                HyperHighM6ACount++;
            }
            if(gene in hashHyperModestM6A){
                HyperModestM6ACount++;
            }
            if(gene in hashHyperLowM6A){
                HyperLowM6ACount++;
            }
            if(gene in hashHypoHighM6A){
                HypoHighM6ACount++;
                hashHypoHighM6AInExpHigh[gene]=1;
            }
            if(gene in hashHypoModestM6A){
                HypoModestM6ACount++;
                hashHypoModestM6AInExpHigh[gene]=1;
            }
            if(gene in hashHypoLowM6A){
                HypoLowM6ACount++;
                hashHypoLowM6AInExpHigh[gene]=1;
            }
        }
        print "##high Expresion Group:";
        print "Type\tCount\tPercentage";
        print "HypoLowM6ACount",HypoLowM6ACount,HypoLowM6ACount/totalM6ACount*100;
        print "HypoModestM6ACount",HypoModestM6ACount,HypoModestM6ACount/totalM6ACount*100;
        print "HypoHighM6ACount",HypoHighM6ACount,HypoHighM6ACount/totalM6ACount*100;
        print "HyperLowM6ACount",HyperLowM6ACount,HyperLowM6ACount/totalM6ACount*100;
        print "HyperModestM6ACount",HyperModestM6ACount,HyperModestM6ACount/totalM6ACount*100;
        print "HyperHighM6ACount",HyperHighM6ACount,HyperHighM6ACount/totalM6ACount*100, "\n";

        HyperHighM6ACount=0;
        HyperModestM6ACount=0;
        HyperLowM6ACount=0;
        HypoHighM6ACount=0;
        HypoModestM6ACount=0;
        HypoLowM6ACount=0;
        for(gene in hashModestExp){
            if(gene in hashHyperHighM6A){
                HyperHighM6ACount++;
            }
            if(gene in hashHyperModestM6A){
                HyperModestM6ACount++;
            }
            if(gene in hashHyperLowM6A){
                HyperLowM6ACount++;
            }
            if(gene in hashHypoHighM6A){
                HypoHighM6ACount++;
            }
            if(gene in hashHypoModestM6A){
                HypoModestM6ACount++;
            }
            if(gene in hashHypoLowM6A){
                HypoLowM6ACount++;
            }
        }
        print "##modest Expresion Group:";
        print "Type\tCount\tPercentage";
        print "HypoLowM6ACount",HypoLowM6ACount,HypoLowM6ACount/totalM6ACount*100;
        print "HypoModestM6ACount",HypoModestM6ACount,HypoModestM6ACount/totalM6ACount*100;
        print "HypoHighM6ACount",HypoHighM6ACount,HypoHighM6ACount/totalM6ACount*100;
        print "HyperLowM6ACount",HyperLowM6ACount,HyperLowM6ACount/totalM6ACount*100;
        print "HyperModestM6ACount",HyperModestM6ACount,HyperModestM6ACount/totalM6ACount*100;
        print "HyperHighM6ACount",HyperHighM6ACount,HyperHighM6ACount/totalM6ACount*100, "\n";

        HyperHighM6ACount=0;
        HyperModestM6ACount=0;
        HyperLowM6ACount=0;
        HypoHighM6ACount=0;
        HypoModestM6ACount=0;
        HypoLowM6ACount=0;
        for(gene in hashLowExp){
            if(gene in hashHyperHighM6A){
                HyperHighM6ACount++;
            }
            if(gene in hashHyperModestM6A){
                HyperModestM6ACount++;
            }
            if(gene in hashHyperLowM6A){
                HyperLowM6ACount++;
            }
            if(gene in hashHypoHighM6A){
                HypoHighM6ACount++;
            }
            if(gene in hashHypoModestM6A){
                HypoModestM6ACount++;
            }
            if(gene in hashHypoLowM6A){
                HypoLowM6ACount++;
            }
        }
        print "##low Expresion Group:";
        print "Type\tCount\tPercentage";
        print "HypoLowM6ACount",HypoLowM6ACount,HypoLowM6ACount/totalM6ACount*100;
        print "HypoModestM6ACount",HypoModestM6ACount,HypoModestM6ACount/totalM6ACount*100;
        print "HypoHighM6ACount",HypoHighM6ACount,HypoHighM6ACount/totalM6ACount*100;
        print "HyperLowM6ACount",HyperLowM6ACount,HyperLowM6ACount/totalM6ACount*100;
        print "HyperModestM6ACount",HyperModestM6ACount,HyperModestM6ACount/totalM6ACount*100;
        print "HyperHighM6ACount",HyperHighM6ACount,HyperHighM6ACount/totalM6ACount*100, "\n";
        print "";
        print "High hypo m6A in high gene expression:";
        print "GeneID\tGeneName\tGeneM6AFC(log2)\tGeneM6ALevelInshCont(log2)\tGeneMeanFPKM";
        for(gene in hashHypoHighM6AInExpHigh){
            print gene,hashHypoHighM6A[gene],hashGeneExp[gene];
        }
        print "";
        print "Modest hypo m6A in high gene expression:";
        print "GeneID\tGeneName\tGeneM6AFC(log2)\tGeneM6ALevelInshCont(log2)\tGeneMeanFPKM";
        for(gene in hashHypoModestM6AInExpHigh){
            print gene,hashHypoModestM6A[gene],hashGeneExp[gene];
        }
        print "";
        print "Low hypo m6A in high gene expression:";
        print "GeneID\tGeneName\tGeneM6AFC(log2)\tGeneM6ALevelInshCont(log2)\tGeneMeanFPKM";
        for(gene in hashHypoLowM6AInExpHigh){
            print gene,hashHypoLowM6A[gene],hashGeneExp[gene];
        }
    }' $m6ADir/HepG2_RNA-seq_gene_m6A_level_FC_shM3vsshCont.txt $geneExpDir/HepG2_RNA-seq_RSEM_shCont_input.txt >> HepG2_RNA-seq_gene_m6A_Exp_shM3_stats.txt

echo -e "## gene expression and m6a FC statistics in HepG2" > HepG2_RNA-seq_gene_m6A_Exp_shWTAP_stats.txt
awk 'BEGIN{OFS="\t";
        m6AHypoLow=-log(1.2)/log(2);m6AHypoModest=-log(1.5)/log(2);m6AHypoHigh=-1;
        m6AHyperLow=log(1.2)/log(2);m6AHyperModest=log(1.5)/log(2);m6AHyperHigh=1;
        highExpCut=4.45516;lowExpCut=1.68182;
        dowmM6ACount=0;upM6ACount=0;totalM6ACount=0;
    }
    ARGIND==1{
        if(FNR>1){
            if($3<=m6AHypoHigh){
                hashHypoHighM6A[$1]=$2"\t"$3"\t"$4;
            }else if(m6AHypoHigh<$3 && $3<m6AHypoModest){
                hashHypoModestM6A[$1]=$2"\t"$3"\t"$4;
            }else if(m6AHypoModest<=$3 && $3<=m6AHypoLow){
                hashHypoLowM6A[$1]=$2"\t"$3"\t"$4;
            }else if(m6AHyperLow<=$3 && $3<=m6AHyperModest){
                hashHyperLowM6A[$1]=$2"\t"$3"\t"$4;
            }else if(m6AHyperModest<$3 && $3<m6AHyperHigh){
                hashHyperModestM6A[$1]=$2"\t"$3"\t"$4;
            }else if($3>=m6AHyperHigh){
                hashHyperHighM6A[$1]=$2"\t"$3"\t"$4;
            }
        }
        totalM6ACount++;
        if($3<0){
            dowmM6ACount++;
        }else if($3>=0){
            upM6ACount++
        }
    }
    ARGIND==2{
        if(FNR>1){
            hashGeneExp[$1]=$6;
            if($7<=lowExpCut){
                hashLowExp[$1]=$6;
            }else if(lowExpCut<$7 && $7<highExpCut){
                hashModestExp[$1]=$6;
            }else if($7>=highExpCut){
                hashHighExp[$1]=$6;
            }
        }
    }
    END{
        print "##gene m6A FC in shWTAP vs. shCont:";
        print "##Determined by Huilin: non-targets = lower than 1.2 fold reduction, hypo-low = 1.2-1.5 fold reduction, hypo-modest = 1.5-2 fold reduction, hypo-high = more than 2 fold reduction";
        print "##ref gene m6A Level (log2) in shCont: boxplot stats (1.000148 1.481705 2.108830 2.962674 5.163127)";
        print "HypoHigh\tHypoModest\tHypoLow\tHypoTotal\tHyperLow\tHyperModest\tHyperHigh\tHyperTotal";
        print length(hashHypoHighM6A),length(hashHypoModestM6A),length(hashHypoLowM6A),dowmM6ACount,length(hashHyperLowM6A),length(hashHyperModestM6A),length(hashHyperHighM6A),upM6ACount;
        print "\n##gene Expression(FPKM) in shCont:";
        print "##Determined by boxplot stats(0.00000 1.68182 3.05976 4.45516 8.58436), log2 transform(MeanFPKM>=1)";
        print "##highExpCut=4.709380,lowExpCut=2.194300";
        print "HighExp\tModestExp\tLowExp";
        print length(hashHighExp),length(hashModestExp),length(hashLowExp),"\n";

        HyperHighM6ACount=0;
        HyperModestM6ACount=0;
        HyperLowM6ACount=0;
        HypoHighM6ACount=0;
        HypoModestM6ACount=0;
        HypoLowM6ACount=0;
        for(gene in hashHighExp){
            if(gene in hashHyperHighM6A){
                HyperHighM6ACount++;
            }
            if(gene in hashHyperModestM6A){
                HyperModestM6ACount++;
            }
            if(gene in hashHyperLowM6A){
                HyperLowM6ACount++;
            }
            if(gene in hashHypoHighM6A){
                HypoHighM6ACount++;
                hashHypoHighM6AInExpHigh[gene]=1;
            }
            if(gene in hashHypoModestM6A){
                HypoModestM6ACount++;
                hashHypoModestM6AInExpHigh[gene]=1;
            }
            if(gene in hashHypoLowM6A){
                HypoLowM6ACount++;
                hashHypoLowM6AInExpHigh[gene]=1;
            }
        }
        print "##high Expresion Group:";
        print "Type\tCount\tPercentage";
        print "HypoLowM6ACount",HypoLowM6ACount,HypoLowM6ACount/totalM6ACount*100;
        print "HypoModestM6ACount",HypoModestM6ACount,HypoModestM6ACount/totalM6ACount*100;
        print "HypoHighM6ACount",HypoHighM6ACount,HypoHighM6ACount/totalM6ACount*100;
        print "HyperLowM6ACount",HyperLowM6ACount,HyperLowM6ACount/totalM6ACount*100;
        print "HyperModestM6ACount",HyperModestM6ACount,HyperModestM6ACount/totalM6ACount*100;
        print "HyperHighM6ACount",HyperHighM6ACount,HyperHighM6ACount/totalM6ACount*100, "\n";

        HyperHighM6ACount=0;
        HyperModestM6ACount=0;
        HyperLowM6ACount=0;
        HypoHighM6ACount=0;
        HypoModestM6ACount=0;
        HypoLowM6ACount=0;
        for(gene in hashModestExp){
            if(gene in hashHyperHighM6A){
                HyperHighM6ACount++;
            }
            if(gene in hashHyperModestM6A){
                HyperModestM6ACount++;
            }
            if(gene in hashHyperLowM6A){
                HyperLowM6ACount++;
            }
            if(gene in hashHypoHighM6A){
                HypoHighM6ACount++;
            }
            if(gene in hashHypoModestM6A){
                HypoModestM6ACount++;
            }
            if(gene in hashHypoLowM6A){
                HypoLowM6ACount++;
            }
        }
        print "##modest Expresion Group:";
        print "Type\tCount\tPercentage";
        print "HypoLowM6ACount",HypoLowM6ACount,HypoLowM6ACount/totalM6ACount*100;
        print "HypoModestM6ACount",HypoModestM6ACount,HypoModestM6ACount/totalM6ACount*100;
        print "HypoHighM6ACount",HypoHighM6ACount,HypoHighM6ACount/totalM6ACount*100;
        print "HyperLowM6ACount",HyperLowM6ACount,HyperLowM6ACount/totalM6ACount*100;
        print "HyperModestM6ACount",HyperModestM6ACount,HyperModestM6ACount/totalM6ACount*100;
        print "HyperHighM6ACount",HyperHighM6ACount,HyperHighM6ACount/totalM6ACount*100, "\n";

        HyperHighM6ACount=0;
        HyperModestM6ACount=0;
        HyperLowM6ACount=0;
        HypoHighM6ACount=0;
        HypoModestM6ACount=0;
        HypoLowM6ACount=0;
        for(gene in hashLowExp){
            if(gene in hashHyperHighM6A){
                HyperHighM6ACount++;
            }
            if(gene in hashHyperModestM6A){
                HyperModestM6ACount++;
            }
            if(gene in hashHyperLowM6A){
                HyperLowM6ACount++;
            }
            if(gene in hashHypoHighM6A){
                HypoHighM6ACount++;
            }
            if(gene in hashHypoModestM6A){
                HypoModestM6ACount++;
            }
            if(gene in hashHypoLowM6A){
                HypoLowM6ACount++;
            }
        }
        print "##low Expresion Group:";
        print "Type\tCount\tPercentage";
        print "HypoLowM6ACount",HypoLowM6ACount,HypoLowM6ACount/totalM6ACount*100;
        print "HypoModestM6ACount",HypoModestM6ACount,HypoModestM6ACount/totalM6ACount*100;
        print "HypoHighM6ACount",HypoHighM6ACount,HypoHighM6ACount/totalM6ACount*100;
        print "HyperLowM6ACount",HyperLowM6ACount,HyperLowM6ACount/totalM6ACount*100;
        print "HyperModestM6ACount",HyperModestM6ACount,HyperModestM6ACount/totalM6ACount*100;
        print "HyperHighM6ACount",HyperHighM6ACount,HyperHighM6ACount/totalM6ACount*100, "\n";
        print "";
        print "High hypo m6A in high gene expression:";
        print "GeneID\tGeneName\tGeneM6AFC(log2)\tGeneM6ALevelInshCont(log2)\tGeneMeanFPKM";
        for(gene in hashHypoHighM6AInExpHigh){
            print gene,hashHypoHighM6A[gene],hashGeneExp[gene];
        }
        print "";
        print "Modest hypo m6A in high gene expression:";
        print "GeneID\tGeneName\tGeneM6AFC(log2)\tGeneM6ALevelInshCont(log2)\tGeneMeanFPKM";
        for(gene in hashHypoModestM6AInExpHigh){
            print gene,hashHypoModestM6A[gene],hashGeneExp[gene];
        }
        print "";
        print "Low hypo m6A in high gene expression:";
        print "GeneID\tGeneName\tGeneM6AFC(log2)\tGeneM6ALevelInshCont(log2)\tGeneMeanFPKM";
        for(gene in hashHypoLowM6AInExpHigh){
            print gene,hashHypoLowM6A[gene],hashGeneExp[gene];
        }
    }' $m6ADir/HepG2_RNA-seq_gene_m6A_level_FC_shWTAPvsshCont.txt $geneExpDir/HepG2_RNA-seq_RSEM_shCont_input.txt >> HepG2_RNA-seq_gene_m6A_Exp_shWTAP_stats.txt
