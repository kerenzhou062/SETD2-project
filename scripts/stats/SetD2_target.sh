#!/bin/sh
echo -e "##Total target gene:$targetGeneNum"

### HepG2
cd /data/zhoukr/hhl_setd2_m6a/HepG2_ChIP-seq
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x shCont/macs1.4/rep1/shCont_rep1_peaks.xls -o /public/zhoukr/test_data/shCont_peaks.bed
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x shSetD2/macs1.4/rep1/shSetD2_rep1_peaks.xls -o /public/zhoukr/test_data/shSetD2_peaks.bed

cd /public/zhoukr/test_data
targetGeneNum=`cat SETD2_target_gene.annotation.bed6 | wc -l`
shContTargetGeneNum=`bedtools intersect -a SETD2_target_gene.annotation.bed6 -b shCont_peaks.bed -wa | cut -f 1,2,3,4,5,6 | sort | uniq | wc -l`
shContPeakNum=`bedtools intersect -a shCont_peaks.bed -b SETD2_target_gene.annotation.bed6 -wa | cut -f 1,2,3,4,5,6 | sort | uniq | wc -l`

shSetD2TargetGeneNum=`bedtools intersect -a SETD2_target_gene.annotation.bed6 -b shSetD2_peaks.bed -wa | cut -f 1,2,3,4,5,6 | sort | uniq | wc -l`
shSetD2PeakNum=`bedtools intersect -a shSetD2_peaks.bed -b SETD2_target_gene.annotation.bed6 -wa | cut -f 1,2,3,4,5,6 | sort | uniq | wc -l`

echo -e "#sample\tshContGeneNum\tshContPeakNum\tshSetD2GeneNum\tshSetD2PeakNum"
echo -e "HepG2_rep1\t${shContTargetGeneNum}\t${shContPeakNum}\t${shSetD2TargetGeneNum}\t${shSetD2PeakNum}"
shContPeakInter=`bedtools intersect -a shCont_peaks.bed -b SETD2_target_gene.annotation.bed6 -wa -wb | cut -f 1,2,3,5,10 | awk 'OFS="\t" {print $1,$2,$3,$5,$4,"HepG2_rep1","shCont"}'`
shSetD2PeakInter=`bedtools intersect -a shSetD2_peaks.bed -b SETD2_target_gene.annotation.bed6 -wa -wb | cut -f 1,2,3,5,10 | awk 'OFS="\t" {print $1,$2,$3,$5,$4,"HepG2_rep1","shSetD2"}'`

cd /data/zhoukr/hhl_setd2_m6a/HepG2_ChIP-seq
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x shCont/macs1.4/rep2/shCont_rep2_peaks.xls -o /public/zhoukr/test_data/shCont_peaks.bed
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x shSetD2/macs1.4/rep2/shSetD2_rep2_peaks.xls -o /public/zhoukr/test_data/shSetD2_peaks.bed

cd /public/zhoukr/test_data
shContTargetGeneNum=`bedtools intersect -a SETD2_target_gene.annotation.bed6 -b shCont_peaks.bed -wa | cut -f 1,2,3,4,5,6 | sort | uniq | wc -l`
shContPeakNum=`bedtools intersect -a shCont_peaks.bed -b SETD2_target_gene.annotation.bed6 -wa | cut -f 1,2,3,4,5,6 | sort | uniq | wc -l`

shSetD2TargetGeneNum=`bedtools intersect -a SETD2_target_gene.annotation.bed6 -b shSetD2_peaks.bed -wa | cut -f 1,2,3,4,5,6 | sort | uniq | wc -l`
shSetD2PeakNum=`bedtools intersect -a shSetD2_peaks.bed -b SETD2_target_gene.annotation.bed6 -wa | cut -f 1,2,3,4,5,6 | sort | uniq | wc -l`

echo -e "HepG2_rep2\t${shContTargetGeneNum}\t${shContPeakNum}\t${shSetD2TargetGeneNum}\t${shSetD2PeakNum}\n"
echo -e "#peakChromssome\tpeakStart\tpeakEnd\tgeneName\tfoldEnrichment\tsample\ttype"
echo -e "${shContPeakInter}\n${shSetD2PeakInter}"
shContPeakInter=`bedtools intersect -a shCont_peaks.bed -b SETD2_target_gene.annotation.bed6 -wa -wb | cut -f 1,2,3,5,10 | awk 'OFS="\t" {print $1,$2,$3,$5,$4,"HepG2_rep2","shCont"}'`
shSetD2PeakInter=`bedtools intersect -a shSetD2_peaks.bed -b SETD2_target_gene.annotation.bed6 -wa -wb | cut -f 1,2,3,5,10 | awk 'OFS="\t" {print $1,$2,$3,$5,$4,"HepG2_rep2","shSetD2"}'`
echo -e "${shContPeakInter}${shSetD2PeakInter}"
echo -e "\n"

### Hela
cd /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x shCont/macs1.4/rep1/shCont_rep1_peaks.xls -o /public/zhoukr/test_data/shCont_peaks.bed
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x shSetD2/macs1.4/rep1/shSetD2_rep1_peaks.xls -o /public/zhoukr/test_data/shSetD2_peaks.bed

cd /public/zhoukr/test_data
targetGeneNum=`cat SETD2_target_gene.annotation.bed6 | wc -l`
shContTargetGeneNum=`bedtools intersect -a SETD2_target_gene.annotation.bed6 -b shCont_peaks.bed -wa | cut -f 1,2,3,4,5,6 | sort | uniq | wc -l`
shContPeakNum=`bedtools intersect -a shCont_peaks.bed -b SETD2_target_gene.annotation.bed6 -wa | cut -f 1,2,3,4,5,6 | sort | uniq | wc -l`

shSetD2TargetGeneNum=`bedtools intersect -a SETD2_target_gene.annotation.bed6 -b shSetD2_peaks.bed -wa | cut -f 1,2,3,4,5,6 | sort | uniq | wc -l`
shSetD2PeakNum=`bedtools intersect -a shSetD2_peaks.bed -b SETD2_target_gene.annotation.bed6 -wa | cut -f 1,2,3,4,5,6 | sort | uniq | wc -l`

echo -e "#sample\tshContGeneNum\tshContPeakNum\tshSetD2GeneNum\tshSetD2PeakNum"
echo -e "Hela_rep1\t${shContTargetGeneNum}\t${shContPeakNum}\t${shSetD2TargetGeneNum}\t${shSetD2PeakNum}"
shContPeakInter=`bedtools intersect -a shCont_peaks.bed -b SETD2_target_gene.annotation.bed6 -wa -wb | cut -f 1,2,3,5,10 | awk 'OFS="\t" {print $1,$2,$3,$5,$4,"Hela_rep1","shCont"}'`
shSetD2PeakInter=`bedtools intersect -a shSetD2_peaks.bed -b SETD2_target_gene.annotation.bed6 -wa -wb | cut -f 1,2,3,5,10 | awk 'OFS="\t" {print $1,$2,$3,$5,$4,"Hela_rep1","shSetD2"}'`

cd /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x shCont/macs1.4/rep2/shCont_rep2_peaks.xls -o /public/zhoukr/test_data/shCont_peaks.bed
macsXlsPeakFilter.pl -pval 1e-5 -fdr 0.05 -x shSetD2/macs1.4/rep2/shSetD2_rep2_peaks.xls -o /public/zhoukr/test_data/shSetD2_peaks.bed

cd /public/zhoukr/test_data
shContTargetGeneNum=`bedtools intersect -a SETD2_target_gene.annotation.bed6 -b shCont_peaks.bed -wa | cut -f 1,2,3,4,5,6 | sort | uniq | wc -l`
shContPeakNum=`bedtools intersect -a shCont_peaks.bed -b SETD2_target_gene.annotation.bed6 -wa | cut -f 1,2,3,4,5,6 | sort | uniq | wc -l`

shSetD2TargetGeneNum=`bedtools intersect -a SETD2_target_gene.annotation.bed6 -b shSetD2_peaks.bed -wa | cut -f 1,2,3,4,5,6 | sort | uniq | wc -l`
shSetD2PeakNum=`bedtools intersect -a shSetD2_peaks.bed -b SETD2_target_gene.annotation.bed6 -wa | cut -f 1,2,3,4,5,6 | sort | uniq | wc -l`

echo -e "Hela_rep2\t${shContTargetGeneNum}\t${shContPeakNum}\t${shSetD2TargetGeneNum}\t${shSetD2PeakNum}\n"
echo -e "#peakChromssome\tpeakStart\tpeakEnd\tgeneName\tfoldEnrichment\tsample\ttype"
echo -e "${shContPeakInter}\n${shSetD2PeakInter}"
shContPeakInter=`bedtools intersect -a shCont_peaks.bed -b SETD2_target_gene.annotation.bed6 -wa -wb | cut -f 1,2,3,5,10 | awk 'OFS="\t" {print $1,$2,$3,$5,$4,"Hela_rep2","shCont"}'`
shSetD2PeakInter=`bedtools intersect -a shSetD2_peaks.bed -b SETD2_target_gene.annotation.bed6 -wa -wb | cut -f 1,2,3,5,10 | awk 'OFS="\t" {print $1,$2,$3,$5,$4,"Hela_rep2","shSetD2"}'`
echo -e "${shContPeakInter}\n${shSetD2PeakInter}"
echo -e "\n"
