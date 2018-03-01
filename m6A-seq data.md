The data processing procedure for m6A-seq data were listed as  below:

---
## Check sequencing reads qualities ##
The qualities of sequencing reads were checked by fastQC software:
```bash
fastqc fastq1 fastq2 ...

```

## Map sequenceing reads to the genome ##
All sequencing reads were aligned to the human (hg19) and mouse (mm10) genomes by using TopHat (version v2.1.1) software with following commands:

```bash
#using tophat to map sequencing
tophat -p 2 -g 1 --library-type=fr-firststrand -G gencode.v24lift37.annotation.gtf -o ./ bowtie2Index_UCSC/hg19 IP.fastq
rename accepted_hits.bam IP.fastq.bam accepted_hits.bam
mv IP.fastq.bam 

#sorted bam (samtools v1.3.1)
samtools sort --threads 16 -m 2G -O bam -o IP.fastq.sorted.bam IP.fastq.bam
rm IP.fastq.bam
samtools index -b IP.fastq.sorted.bam

```

## Gernerate tiled data file(TDF) ##
We used IGVtools (version 2.3.98) to generate TDF of each sample:
```bash
igvtools count -z 7 -w 25 -e 250 shCont_input_rep1.fastq.sorted.bam shCont_input_rep1.cov.tdf hg19.chrom.sizes > /dev/null 2>&1

```

## Peak calling ##
We used exomePeak (version 2.8.0) to detect enriched m6A sites with suggested parameters:
```R
#!/usr/bin/env Rscript
#exomepeak Script 2
#R script
#Define parameters and load library
library("exomePeak")
setwd("./shCont/")

INPUT1_BAM = "HepG2_m6A-seq_shCont_input_rep1.fastq.sorted.bam"
INPUT2_BAM = "HepG2_m6A-seq_shCont_input_rep2.fastq.sorted.bam"
INPUT3_BAM = "HepG2_m6A-seq_shCont_input_rep3.fastq.sorted.bam"
IP1_BAM = "HepG2_m6A-seq_shCont_IP_rep1.fastq.sorted.bam"
IP2_BAM = "HepG2_m6A-seq_shCont_IP_rep2.fastq.sorted.bam"
IP3_BAM = "HepG2_m6A-seq_shCont_IP_rep3.fastq.sorted.bam"

GTF_ANNO = "/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.gtf"
#peak calling
exomepeak(GENE_ANNO_GTF=GTF_ANNO, IP_BAM=c(IP1_BAM, IP2_BAM), INPUT_BAM=c(INPUT1_BAM, INPUT2_BAM), EXPERIMENT_NAME="shCont_m6A_rep1-2")
exomepeak(GENE_ANNO_GTF=GTF_ANNO, IP_BAM=c(IP1_BAM, IP3_BAM), INPUT_BAM=c(INPUT1_BAM, INPUT3_BAM), EXPERIMENT_NAME="shCont_m6A_rep1-3")
exomepeak(GENE_ANNO_GTF=GTF_ANNO, IP_BAM=c(IP2_BAM, IP3_BAM), INPUT_BAM=c(INPUT2_BAM, INPUT3_BAM), EXPERIMENT_NAME="shCont_m6A_rep2-3")
exomepeak(GENE_ANNO_GTF=GTF_ANNO, IP_BAM=c(IP1_BAM, IP2_BAM, IP3_BAM), INPUT_BAM=c(INPUT1_BAM, INPUT2_BAM, INPUT3_BAM), EXPERIMENT_NAME="shCont_m6A_rep1-2-3")
exomepeak(GENE_ANNO_GTF=GTF_ANNO, IP_BAM=c(IP1_BAM), INPUT_BAM=c(INPUT1_BAM), EXPERIMENT_NAME="shCont_m6A_rep1")
exomepeak(GENE_ANNO_GTF=GTF_ANNO, IP_BAM=c(IP2_BAM), INPUT_BAM=c(INPUT2_BAM), EXPERIMENT_NAME="shCont_m6A_rep2")
exomepeak(GENE_ANNO_GTF=GTF_ANNO, IP_BAM=c(IP3_BAM), INPUT_BAM=c(INPUT3_BAM), EXPERIMENT_NAME="shCont_m6A_rep3")

```

## Venn plots test of m6A peaks across samples ##
We gernerate peaks without strand and duplicate peaks
```bash
exomePeakToSummit.pl -p 1e-5 -fdr 1e-3 -fold 2 -original -x  shCont_m6A_rep1/con_peak.xls -o HepG2_shCont_m6a.tmp
awk '{key=$1""$2""$3; hash[key]=$1"\t"$2"\t"$3"\t"$6;}END{for(i in hash){print hash[i];}}' OFS="\t" HepG2_shCont_m6a.tmp | sort | uniq > HepG2_shCont_m6a_rep1.bed

```

We use ChIPpeakAnno R package to draw Venn plots of across samples
```R
bedFile <- read.table('HepG2_shCont_m6a_rep1.bed', header=FALSE)
shCont_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

bedFile <- read.table('HepG2_shSetD2_m6a_rep1.bed', header=FALSE)
shSetD2_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

bedFile <- read.table('HepG2_shM14_m6a_rep1.bed', header=FALSE)
shM14_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

bedFile <- read.table('HepG2_shM3_m6a_rep1.bed', header=FALSE)
shM3_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

bedFile <- read.table('HepG2_shWTAP_m6a_rep1.bed', header=FALSE)
shWTAP_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3), strand = bedFile$V4)

pdf("HepG2_m6a_shTargets_rep1_venn.pdf", width=8)
makeVennDiagram(Peaks=list(shCont_granges, shSetD2_granges, shM14_granges, shM3_granges, shWTAP_granges),
  NameOfPeaks=c("shNS", "shSETD2", "shMETTL14", "shMETTL3", "shWTAP"), margin=0.05,
  fill=c("orange","brown","cyan","green","magenta"), fontface = "bold", alpha = 0.5, cex=1, cat.cex=1)
dev.off()


```

## Definition of writer- and Setd2-responsive m6A sites ##
```bash
#fold change of intersected peaks
log2(FC) = log2((shSetd2+0.1) / (shCont+0.1))
#cutoff for Setd2-responsive sites
log2(FC) < log2(0.5)

awk '{if(FNR>1){if($14<(1/2)) {for(i=1;i<12;i++) {printf("%s\t",$i);} printf("%s", $12); printf "\n";}}}' HepG2_m6a_peak_FC.txt > HepG2_m6a_shCont_SetD2_responsive.bed12
awk '{if(FNR>1){if($15<(1/2)) {for(i=1;i<12;i++) {printf("%s\t",$i);} printf("%s", $12); printf "\n";}}}' HepG2_m6a_peak_FC.txt > HepG2_m6a_shCont_M14_responsive.bed12

#responsive sites in shSetd2 samples
bedtools intersect -a HepG2_shSetD2_m6a.bed12 -b HepG2_m6a_shCont_M14_responsive.bed12 -s -split -u | sort  > ./temp.bed12
bedtools intersect -a temp.bed12 -b HepG2_m6a_shCont_SetD2_responsive.bed12 -s -split -v | sort | uniq > HepG2_m6a_shSetD2_M14_responsive.bed12

#identify share-writer responsive sites
awk '{if(FNR>1){if($15<(1/2) && $16<(1/2) && $17<(1/2)) {for(i=1;i<12;i++) {printf("%s\t",$i);} printf("%s", $12); printf "\n";}}}' ./../xls/allPeak/HepG2_m6a_peak_FC.txt > HepG2_m6a_shCont_writerShare_responsive.bed12
bedtools intersect -a HepG2_shSetD2_m6a.bed12 -b HepG2_m6a_shCont_writerShare_responsive.bed12 -s -split -u | sort  > ./temp.bed12
bedtools intersect -a temp.bed12 -b HepG2_m6a_shCont_SetD2_responsive.bed12 -s -split -v | sort | uniq > HepG2_m6a_shSetD2_writerShare_responsive.bed12

#identify all-writer responsive sites
awk '{if(FNR>1){if($15<(1/2) || $16<(1/2) || $17<(1/2)) {for(i=1;i<12;i++) {printf("%s\t",$i);} printf("%s", $12); printf "\n";}}}' ./../xls/allPeak/HepG2_m6a_peak_FC.txt > HepG2_m6a_shCont_writerAll_responsive.bed12
bedtools intersect -a HepG2_shSetD2_m6a.bed12 -b HepG2_m6a_shCont_writerAll_responsive.bed12 -s -split -u | sort  > ./temp.bed12
bedtools intersect -a temp.bed12 -b HepG2_m6a_shCont_SetD2_responsive.bed12 -s -split -v | sort | uniq > HepG2_m6a_shSetD2_writerAll_responsive.bed12

#identify non-writer responsive sites
cat HepG2_m6a_shCont_M14_responsive.bed12 HepG2_m6a_shCont_M3_responsive.bed12 HepG2_m6a_shCont_WTAP_responsive.bed12 | sort -t $'\t' -k1,1V -k2,2n | uniq > HepG2_shCont_allWriter_responsive.bed12
awk '{if(FNR>1){if($15>=(1/2)&&$16>=(1/2)&&$17>=(1/2)) {for(i=1;i<12;i++) {printf("%s\t",$i);} printf("%s", $12); printf "\n";}}}' ./../xls/allPeak/HepG2_m6a_peak_FC.txt > HepG2_m6a_shCont_nonWriter_responsive.bed12
bedtools intersect -a HepG2_shSetD2_m6a.bed12 -b HepG2_m6a_shCont_nonWriter_responsive.bed12 -s -split -u | sort  > ./temp.bed12
bedtools intersect -a temp.bed12 -b HepG2_m6a_shCont_SetD2_responsive.bed12 -s -split -v | sort | uniq > HepG2_m6a_shSetD2_nonWriter_responsive.bed12

#intersect with H3K36me3 sites
bedtools intersect -a HepG2_m6a_shCont_M14_responsive.bed12 -b HepG2_macs_shCont_peaks.bed -wa | sort -t $'\t' -k 1,1 -k 2,2n | uniq > HepG2_m6a_shCont_M14_responsive_H3K36me3.bed12
bedtools intersect -a HepG2_m6a_shSetD2_M14_responsive.bed12 -b HepG2_macs_shCont_peaks.bed -wa | sort -t $'\t' -k 1,1 -k 2,2n | uniq > HepG2_m6a_shSetD2_M14_responsive_H3K36me3.bed12

#related bash scripts
*m6A-seq_FC.sh
*m6A-seq_bin_FC.sh
```

## draw the distribution of m6A on mRNA genes ##
```bash
#draw bin distribution based on count and bed12, overlapped with H3K36me3 peaks
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shCont_M14_responsive_H3K36me3.bed12 -bed6 gencode.v24lift37.annotation.mRNA.longest.exon.bed6 -o ./HepG2_m6a_shCont_M14_responsive_H3K36me3.bed12bin
bedBinDistribution.pl -strand -smooth move -span 10 -t count --input HepG2_m6a_shSetD2_M14_responsive_H3K36me3.bed12 -bed6 gencode.v24lift37.annotation.mRNA.longest.exon.bed6 -o ./HepG2_m6a_shSetD2_M14_responsive_H3K36me3.bed12bin
paste HepG2_m6a_shCont_M14_responsive_H3K36me3.bed12bin HepG2_m6a_shSetD2_M14_responsive_H3K36me3.bed12bin | cut -f 1,2,3,6 > HepG2_m6a_M14_responsive_H3K36me3.bed12bin
sed -i '1i region\tbin\tshCont\tshSetD2' HepG2_m6a_M14_responsive_H3K36me3.bed12bin
rm -f HepG2_m6a_shCont_M14_responsive_H3K36me3.bed12bin HepG2_m6a_shSetD2_M14_responsive_H3K36me3.bed12bin

#gene type counts and region counts of SetD2 responsive m6a sites
#using regionDistribution.pl and geneDistribution.pl
geneDistribution.pl -strand --input HepG2_m6a_shCont_SetD2_responsive.bed12 -bed6 gencode.v24lift37.annotation.all.longest.exon.bed6 -o HepG2_m6a_shCont_SetD2_responsive.gene
sed -i '1i geneType\tpeakNumber' HepG2_m6a_shCont_SetD2_responsive.gene
regionDistribution.pl -strand -size 200 -f '5utr,cds,stopCodon,3utr' --input HepG2_m6a_shCont_SetD2_responsive.bed12 -bed6 gencode.v24lift37.annotation.mRNA.longest.exon.bed6 -o HepG2_m6a_shCont_SetD2_responsive.region
sed -i '1i region\tpeakNumber\tenrichment' HepG2_m6a_shCont_SetD2_responsive.region

#related bash scripts
*m6A-seq_bin_FC.sh

```

## Venn plots test of m6A peaks across writer- and Setd2-responsive sites ##
We gernerate responsive sites based on peak_FC.txt.
```bash
awk 'BEGIN{OFS="\t";};{ if(NR>1){ cutoff=1/2;if($14<cutoff){ gsub(/\.[0-9]+/,"" ,$4);print $4"-"NR; } } }' HepG2_m6a_peak_FC.txt | sort | uniq > ./HepG2_shSetD2_FC_gene.txt
awk 'BEGIN{OFS="\t";};{ if(NR>1){ cutoff=1/2;if($15<cutoff){ gsub(/\.[0-9]+/,"" ,$4);print $4"-"NR; } } }' HepG2_m6a_peak_FC.txt | sort | uniq > ./HepG2_shM14_FC_gene.txt
awk 'BEGIN{OFS="\t";};{ if(NR>1){ cutoff=1/2;if($16<cutoff){ gsub(/\.[0-9]+/,"" ,$4);print $4"-"NR; } } }' HepG2_m6a_peak_FC.txt | sort | uniq > ./HepG2_shM3_FC_gene.txt
awk 'BEGIN{OFS="\t";};{ if(NR>1){ cutoff=1/2;if($17<cutoff){ gsub(/\.[0-9]+/,"" ,$4);print $4"-"NR; } } }' HepG2_m6a_peak_FC.txt | sort | uniq > ./HepG2_shWTAP_FC_gene.txt

#related bash scripts
*m6A-seq_FC_Venn.sh
```

Then we used VennDiagram R package (codes in Rscripts/HepG2/draw_venn.r) to draw Venn plots.

## Cumulative plots of writer-repsonsive sites ##
We gernerate responsive sites based on peak_FC.txt.
```bash
#using cumulativePlot.pl
awk '{if(NR>1){cutoff=log(0.5)/log(2); if ($15 < cutoff) {print $14;}}}' HepG2_m6a_peak_FC_log2.txt > HepG2_shM14_FC_log2.txt
awk '{if(NR>1){cutoff=log(0.5)/log(2); if ($16 < cutoff) {print $14;}}}' HepG2_m6a_peak_FC_log2.txt > HepG2_shM3_FC_log2.txt
awk '{if(NR>1){cutoff=log(0.5)/log(2); if ($17 < cutoff) {print $14;}}}' HepG2_m6a_peak_FC_log2.txt > HepG2_shWTAP_FC_log2.txt
awk '{if(NR>1){cutoff=log(0.5)/log(2); if (($15 < cutoff)&&($16 < cutoff)&&($17 < cutoff)) {print $14;}}}' HepG2_m6a_peak_FC_log2.txt > HepG2_shareTargets_FC_log2.txt
awk '{if(NR>1){cutoff=log(0.5)/log(2); if (($15 >= cutoff)&&($16 >= cutoff)&&($17 >= cutoff)) {print $14;}}}' HepG2_m6a_peak_FC_log2.txt > HepG2_nonTargets_FC_log2.txt

cumulativePlot.pl -header false -interval 0.01 --input HepG2_shM14_FC_log2.txt HepG2_shM3_FC_log2.txt HepG2_shWTAP_FC_log2.txt HepG2_shareTargets_FC_log2.txt HepG2_nonTargets_FC_log2.txt -o HepG2_cumulativePlot.txt
sed -i '1i FC\tshM14\tshM3\tshWTAP\tshare\tnonTarget' HepG2_cumulativePlot.txt

#related bash scripts
*m6A-seq_FC_cumulative.sh
```

## Circos plot ##
We used circos.pl to draw circos plot. The coordinates of exomePeak results (in bed12 format) and their corresonding values of fold Encichment were adopted as the genomic intervals and corresponding values.
```bash
exome2PeaksFC.pl -cutoff 0.1 -nonOverlapA -nonOverlapB -x HepG2_shSetD2_m6a.xls HepG2_shCont_m6a.xls -o temp.xls
awk '{gsub(/chr/,"hs",$1);$13=log(1/$13)/log(2);print $1,$2,$3,$13;}' OFS="\t" temp.xls > HepG2_m6a_shSetD2_FC_circos.txt

partial_ideogram_conf="
<ideogram>\n
\n
<spacing>\n
default = 0.002r\n
break   = 0.2r\n
</spacing>\n
\n
<<include ideogram.position.conf>>\n
<<include ideogram.label.conf>>\n
<<include bands.conf>>\n
\n
radius*       = 0.92r\n
\n
</ideogram>\n

"
partial_ideogram_label_conf="
show_label       = yes\n
label_font       = bold\n
label_radius     = dims(image,radius)-95p\n
label_size       = 80\n
label_parallel   = yes\n
label_case       = lower\n
label_format     = eval(sprintf(\"chr%s\",var(label)))\n
"
all_ideogram_conf="
<ideogram>\n
\n
<spacing>\n
default = 0.002r\n
break   = 0.2r\n
</spacing>\n
\n
<<include ideogram.position.conf>>\n
<<include ideogram.label.conf>>\n
<<include bands.conf>>\n
\n
radius*       = 0.92r\n
\n
</ideogram>\n

"
all_ideogram_label_conf="
show_label       = yes\n
label_font       = bold\n
label_radius     = dims(image,radius)-95p\n
label_size       = 80\n
label_parallel   = yes\n
label_case       = upper\n
label_format     = eval(sprintf(\"%s\",var(label)))\n
"


cd /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/circos/conf/histogram/
echo -e $partial_ideogram_conf > ideogram.conf
echo -e $partial_ideogram_label_conf > ideogram.label.conf
circos -conf histogram.partial.FC.conf > /dev/null 2>&1
echo -e $all_ideogram_conf > ideogram.conf
echo -e $all_ideogram_label_conf > ideogram.label.conf
circos -conf histogram.FC.conf > /dev/null 2>&1
rm -f histogram.partial.FC.svg histogram.FC.svg

#related bash scrits
*m6A-seq_circos_bin.sh

```

## Correlation between H3K36me3 sites, heterochromatin regions and m6A sites ##
Teterochromatin regions were usually modified by H3K9me3 and H3K27me3. The H3K9me3 and H3K27me3 peak data were downloaded from ENCODE project (Accessions:ENCFF533JQH (H3K9me3) and ENCFF042EDV (H3K27me3)).

The regioneR R package was used to perform permutation test between regions of H3K36me3 and m6A or H3K9me3 and m6A.

```bash
HepG2_H3K36me3=/data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/bed6/HepG2_macs_shCont_peaks.bed
chipseq=/data/zhoukr/hhl_setd2_m6a/others/encode/
bed12Path=/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/bed12
num=`bedtools intersect -a $bed12Path/HepG2_shCont_m6a.bed12 -b $chipseq/narrowPeak_H3K9me3.bed -wa |sort|uniq| wc -l`
num=`bedtools intersect -a $bed12Path/HepG2_shCont_m6a.bed12 -b $HepG2_H3K36me3 -wa |sort|uniq| wc -l`

awk '{print $1, $2, $3}' $HepG2_H3K36me3 > HepG2_macs_shCont_peaks.bed
awk '{print $1, $2, $3}' $bed12Path/HepG2_shCont_m6a.bed12 > HepG2_shCont_m6a.bed
awk '{print $1, $2, $3}' $chipseq/narrowPeak_H3K9me3.bed > narrowPeak_H3K9me3.bed

#related bash scripts
HepG2_m6a_seq_heterochromatin.sh
```

Permutation test
```R
#!/usr/bin/env Rscript
#exomepeak Script 2
#R script
#Define parameters and load library
setwd("/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/heterochomatin")
library("rtracklayer")
library("GenomicRanges")
library("regioneR")
genome <- filterChromosomes(getGenome("hg19"))
m6a <- import("HepG2_shCont_m6a.bed")
H3K36me3 <- import("HepG2_macs_shCont_peaks.bed")
H3K9me3 <- import("narrowPeak_H3K9me3.bed")
m6a_H3K36me3 <- overlapPermTest(A=m6a, B=H3K36me3, ntimes=1000, genome=genome, verbose=FALSE, alternative='greater')

sink("m6a_H3K36me3_permutation_test.txt")
print(summary(m6a_H3K36me3))
sink()
pdf("m6a_H3K36me3_permutation_test.pdf")
plot(m6a_H3K36me3)
dev.off()

m6a_H3K9me3 <- overlapPermTest(A=m6a, B=H3K9me3, ntimes=1000, genome=genome, verbose=FALSE, alternative='greater')
sink("m6a_H3K9me3_permutation_test.txt")
print(summary(m6a_H3K9me3))
sink()
pdf("m6a_H3K9me3_permutation_test.pdf")
plot(m6a_H3K9me3)
dev.off()

```
