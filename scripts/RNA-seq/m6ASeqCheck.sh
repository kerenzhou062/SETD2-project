#!/bin/sh
RNASeqDir=/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/
statsFile=/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/m6A_gene_correlation_stats.txt
## Hela
perl m6ASeqCorrelationCheck.pl input $RNASeqDir/Hela/RSEM $RNASeqDir/Hela/correlationStats
perl m6ASeqCorrelationCheck.pl IP $RNASeqDir/Hela/RSEM $RNASeqDir/Hela/correlationStats
perl ./m6ASeqSampleRepCorrelationCheck.pl $RNASeqDir/Hela/RSEM $RNASeqDir/Hela/correlationStats

## hepG2
perl m6ASeqCorrelationCheck.pl input $RNASeqDir/HepG2/RSEM $RNASeqDir/HepG2/correlationStats
perl m6ASeqCorrelationCheck.pl IP $RNASeqDir/HepG2/RSEM $RNASeqDir/HepG2/correlationStats
perl ./m6ASeqSampleRepCorrelationCheck.pl $RNASeqDir/HepG2/RSEM $RNASeqDir/HepG2/correlationStats

echo "#stats of gene level(FPKM) correlations between replicates" > $statsFile
echo "#Hela control rep4 and rep5 were retrive from Hechuan, (pubmed:23453015, Niu Y, et al. Genomics Proteomics Bioinformatics 2013)" >> $statsFile
echo "HepG2-IP:" >> $statsFile
cat $RNASeqDir/HepG2/correlationStats/m6A-seq_IP_correlationStats.txt >> $statsFile
echo "" >> $statsFile
echo "Hela-IP:" >> $statsFile
cat $RNASeqDir/Hela/correlationStats/m6A-seq_IP_correlationStats.txt >> $statsFile
echo -e "\n\n" >> $statsFile
echo "HepG2-input:" >> $statsFile
cat $RNASeqDir/HepG2/correlationStats/m6A-seq_input_correlationStats.txt >> $statsFile
echo "" >> $statsFile
echo "Hela-input:" >> $statsFile
cat $RNASeqDir/Hela/correlationStats/m6A-seq_input_correlationStats.txt >> $statsFile

HelaRsem=/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/Hela/RSEM
#### rsem summary
cd $HelaRsem
echo -e "GeneID\tgeneName\tshCont-rep1\tshCont-rep2\tshCont-rep3\tshSetD2-rep1\tshSetD2-rep2\tshSetD2-rep3\tshM14-rep1\tshM14-rep2\tshM14-rep3\tshM3-rep1\tshM3-rep2\tshM3-rep3\tshWTAP-rep1\tshWTAP-rep2\tshWTAP-rep3" \
  > $RNASeqDir/Hela/correlationStats/Hela_input_rsem_summary.txt

awk 'BEGIN{OFS="\t";}
    ARGIND==1{
        if(FNR>1){
            EXPa[$1]= $2"\t"$3"\t"$4"\t"$5;
        }
    }
    ARGIND==2{
        if(FNR>1){
            EXPb[$1]= $3"\t"$4"\t"$5;
        }
    }
    ARGIND==3{
        if(FNR>1){
            EXPc[$1]= $3"\t"$4"\t"$5;
        }
    }
    ARGIND==4{
        if(FNR>1){
            EXPd[$1]= $3"\t"$4"\t"$5;
        }
    }
    ARGIND==5{
        if(FNR>1){
            EXPe[$1]= $3"\t"$4"\t"$5;
        }
    }
    END{
        slen=asorti(EXPa,EXPaS);
        for (i=1;i<=slen;i++){
            id=EXPaS[i];
            if(id in EXPa){
                print id,EXPa[id],EXPb[id],EXPc[id],EXPd[id],EXPe[id];
            }
        }
    }' Hela_RNA-seq_RSEM_shCont_input.txt Hela_RNA-seq_RSEM_shSetD2_input.txt Hela_RNA-seq_RSEM_shM14_input.txt \
      Hela_RNA-seq_RSEM_shM3_input.txt Hela_RNA-seq_RSEM_shWTAP_input.txt >> $RNASeqDir/Hela/correlationStats/Hela_input_rsem_summary.txt

echo -e "GeneID\tgeneName\tshCont-rep1\tshCont-rep2\tshCont-rep3\tshSetD2-rep1\tshSetD2-rep2\tshSetD2-rep3\tshM14-rep1\tshM14-rep2\tshM14-rep3\tshM3-rep1\tshM3-rep2\tshM3-rep3\tshWTAP-rep1\tshWTAP-rep2\tshWTAP-rep3" \
  > $RNASeqDir/Hela/correlationStats/Hela_IP_rsem_summary.txt
awk 'BEGIN{OFS="\t";}
    ARGIND==1{
        if(FNR>1){
            EXPa[$1]= $2"\t"$3"\t"$4"\t"$5;
        }
    }
    ARGIND==2{
        if(FNR>1){
            EXPb[$1]= $3"\t"$4"\t"$5;
        }
    }
    ARGIND==3{
        if(FNR>1){
            EXPc[$1]= $3"\t"$4"\t"$5;
        }
    }
    ARGIND==4{
        if(FNR>1){
            EXPd[$1]= $3"\t"$4"\t"$5;
        }
    }
    ARGIND==5{
        if(FNR>1){
            EXPe[$1]= $3"\t"$4"\t"$5;
        }
    }
    END{
        slen=asorti(EXPa,EXPaS);
        for (i=1;i<=slen;i++){
            id=EXPaS[i];
            if(id in EXPa){
                print id,EXPa[id],EXPb[id],EXPc[id],EXPd[id],EXPe[id];
            }
        }
    }' Hela_RNA-seq_RSEM_shCont_IP.txt Hela_RNA-seq_RSEM_shSetD2_IP.txt Hela_RNA-seq_RSEM_shM14_IP.txt \
      Hela_RNA-seq_RSEM_shM3_IP.txt Hela_RNA-seq_RSEM_shWTAP_IP.txt >> $RNASeqDir/Hela/correlationStats/Hela_IP_rsem_summary.txt

HepG2Rsem=/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/HepG2/RSEM
#### rsem summary
cd $HepG2Rsem
echo -e "GeneID\tgeneName\tshCont-rep1\tshCont-rep2\tshCont-rep3\tshSetD2-rep1\tshSetD2-rep2\tshSetD2-rep3\tshM14-rep1\tshM14-rep2\tshM14-rep3\tshM3-rep1\tshM3-rep2\tshM3-rep3\tshWTAP-rep1\tshWTAP-rep2\tshWTAP-rep3" \
  > $RNASeqDir/HepG2/correlationStats/HepG2_input_rsem_summary.txt
awk 'BEGIN{OFS="\t";}
    ARGIND==1{
        if(FNR>1){
            EXPa[$1]= $2"\t"$3"\t"$4"\t"$5;
        }
    }
    ARGIND==2{
        if(FNR>1){
            EXPb[$1]= $3"\t"$4"\t"$5;
        }
    }
    ARGIND==3{
        if(FNR>1){
            EXPc[$1]= $3"\t"$4"\t"$5;
        }
    }
    ARGIND==4{
        if(FNR>1){
            EXPd[$1]= $3"\t"$4"\t"$5;
        }
    }
    ARGIND==5{
        if(FNR>1){
            EXPe[$1]= $3"\t"$4"\t"$5;
        }
    }
    END{
        slen=asorti(EXPa,EXPaS);
        for (i=1;i<=slen;i++){
            id=EXPaS[i];
            if(id in EXPa){
                print id,EXPa[id],EXPb[id],EXPc[id],EXPd[id],EXPe[id];
            }
        }
    }' HepG2_RNA-seq_RSEM_shCont_input.txt HepG2_RNA-seq_RSEM_shSetD2_input.txt HepG2_RNA-seq_RSEM_shM14_input.txt \
      HepG2_RNA-seq_RSEM_shM3_input.txt HepG2_RNA-seq_RSEM_shWTAP_input.txt >> $RNASeqDir/HepG2/correlationStats/HepG2_input_rsem_summary.txt

echo -e "GeneID\tgeneName\tshCont-rep1\tshCont-rep2\tshCont-rep3\tshSetD2-rep1\tshSetD2-rep2\tshSetD2-rep3\tshM14-rep1\tshM14-rep2\tshM14-rep3\tshM3-rep1\tshM3-rep2\tshM3-rep3\tshWTAP-rep1\tshWTAP-rep2\tshWTAP-rep3" \
  > $RNASeqDir/HepG2/correlationStats/HepG2_IP_rsem_summary.txt
awk 'BEGIN{OFS="\t";}
    ARGIND==1{
        if(FNR>1){
            EXPa[$1]= $2"\t"$3"\t"$4"\t"$5;
        }
    }
    ARGIND==2{
        if(FNR>1){
            EXPb[$1]= $3"\t"$4"\t"$5;
        }
    }
    ARGIND==3{
        if(FNR>1){
            EXPc[$1]= $3"\t"$4"\t"$5;
        }
    }
    ARGIND==4{
        if(FNR>1){
            EXPd[$1]= $3"\t"$4"\t"$5;
        }
    }
    ARGIND==5{
        if(FNR>1){
            EXPe[$1]= $3"\t"$4"\t"$5;
        }
    }
    END{
        slen=asorti(EXPa,EXPaS);
        for (i=1;i<=slen;i++){
            id=EXPaS[i];
            if(id in EXPa){
                print id,EXPa[id],EXPb[id],EXPc[id],EXPd[id],EXPe[id];
            }
        }
    }' HepG2_RNA-seq_RSEM_shCont_IP.txt HepG2_RNA-seq_RSEM_shSetD2_IP.txt HepG2_RNA-seq_RSEM_shM14_IP.txt \
      HepG2_RNA-seq_RSEM_shM3_IP.txt HepG2_RNA-seq_RSEM_shWTAP_IP.txt >> $RNASeqDir/HepG2/correlationStats/HepG2_IP_rsem_summary.txt

