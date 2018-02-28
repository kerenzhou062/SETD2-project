The data processing procedure for ChIP-seq data were as listed below:

---

## Map sequenceing reads to the genome ##
All sequencing reads were aligned to the human (hg19) and mouse (mm10) genomes by using bowtie (version 1.1.2) software with following commands:

```bash
# using bowtie to map sequencing
bowtie -t -p 16 -v 2 -m 1 --best --strata --sam bowtieIndex_UCSC/hg19 fastq.map.sam
samtools view -h -bS -F 4 --threads 32 fastq.map.sam fastq.map.bam

# convert sam to bam (samtools v1.3.1)
samtools sort --threads 16 -m 2G -O bam -o fastqMap.sorted.bam fastq.map.bam
samtools index -b fastq.map.sorted.bam

# delete temp files
rm -f fastq.map.sam fastq.map.bam
```

## Peak calling ##
We used MACS (version 1.4.2) to detect enriched H3K36me3 sites withwith the suggested parameters for peak calling of H3K36me3:
```bash
# using samtools to merge replicates
samtools merge fastq.merge.bam fastq.rep1.bam fastq.rep2.bam

# identifing peaks of H3K36me3 from mapped bam
macs14 -t IP.fastq.merge.bam -c input.fastq.merge.bam -g hs -n shM3 --nomodel --shiftsize 147 -B -S --call-subpeaks > macs14.log 2>&1 &

# Cutoff for significantly enriched peak
p-value < 1e-5 
fdr < 0.05

# macs2 is an alternative option 
# but not used in this project
macs2 callpeak --broad -B -g hs -t IP_rep1.fastq.sorted.bam IP_rep2.fastq.rep2.bam -c input_rep1.fastq.sorted.bam input_rep2.fastq.sorted.bam  -n shCont_macs2 > shCont.macs2.log 2>&1 &

```

> refference:
> 1. Identifying ChIP-seq enrichment using MACS, Nat Protoc., 2012
> 2. shiftsize was set as recommended by "Use Model-Based Analysis of ChIP-Seq (MACS) to Analyze Short Reads Generated by Sequencing Protein–DNA Interactions in Embryonic Stem Cells")


## Definition of Setd2-responsive H3K36me3 sites ##
```bash
# fold change of intersected peaks
log2(FC) = log2((shSetd2+0.1) / (shCont+0.1))

# related bash scripts
*chip-seq_FC.sh

```

## draw the distribution of H3K36me3 on mRNA genes ##
```bash
# using bedBinDistribution.pl
bedBinDistribution.pl --feature 'promoter,5utr,cds,3utr' -span 5 -smooth move -t count --input H3K36me3_macs_peaks.bed6 -bed6 gencode.v24lift37.annotation.mRNA.longest.exon+promoter.bed6 -o H3K36me3_macs_peaks.bin 

# using regionDistribution.pl
regionDistribution.pl --feature 'promoter,5utr,cds,stopCodon,3utr' --input H3K36me3_macs_peaks.bed6 -bed6 gencode.v24lift37.annotation.mRNA.longest.exon+promoter.bed6 -o H3K36me3_macs_peaks.region

#related bash scripts
*chip-seq_bin.sh

```

## Calculate read counts of genes ##
We used featureCounts to calculate the read counts of genes in IP and input samples.
```bash
# using featureCounts
featureCounts -T 10 -a gencode.v24lift37.annotation.gtf -g gene_id -F GTF -t gene -M -o counts.txt bam1 bam2...

#related bash scrits
*DEGseq.sh

```

## Calculate Fold Changes of genes ##
We used DEGSeq to calculate the Fold Changes of genes between IP and input samples. The FC will be defined as modification level of H3K36me3 of genes.
```R
#!/usr/bin/env Rscript
# exomepeak Script 2
# R script
# Define parameters and load library
library("DEGseq")
directory="/data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/DEGseq"
setwd(directory)
featureCount="HepG2_ChIP-seq_featureCount.txt"
colCount <- ncol(read.table(file=featureCount, nrows=1))
countMatrix <- readGeneExp(file=featureCount, geneCol=1, valCol=c(7:colCount))

DEGexp(geneExpMatrix1=countMatrix, geneCol1=1, expCol1=c(2,3), groupLabel1="ip",
    geneExpMatrix2=countMatrix, geneCol2=1, expCol2=c(4,5), groupLabel2="input",
    method="MARS", outputDir=paste(directory, "shCont", sep="/"))

DEGexp(geneExpMatrix1=countMatrix, geneCol1=1, expCol1=c(6,7), groupLabel1="ip",
    geneExpMatrix2=countMatrix, geneCol2=1, expCol2=c(8,9), groupLabel2="input",
    method="MARS", outputDir=paste(directory, "shSetD2", sep="/"))

DEGexp(geneExpMatrix1=countMatrix, geneCol1=1, expCol1=c(10), groupLabel1="ip",
    geneExpMatrix2=countMatrix, geneCol2=1, expCol2=c(11), groupLabel2="input",
    method="MARS", outputDir=paste(directory, "shM14", sep="/"))

DEGexp(geneExpMatrix1=countMatrix, geneCol1=1, expCol1=c(12), groupLabel1="ip",
    geneExpMatrix2=countMatrix, geneCol2=1, expCol2=c(13), groupLabel2="input",
    method="MARS", outputDir=paste(directory, "shM3", sep="/"))

DEGexp(geneExpMatrix1=countMatrix, geneCol1=1, expCol1=c(14), groupLabel1="ip",
    geneExpMatrix2=countMatrix, geneCol2=1, expCol2=c(15), groupLabel2="input",
    method="MARS", outputDir=paste(directory, "shWTAP", sep="/"))
    
```