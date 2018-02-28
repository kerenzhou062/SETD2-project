#!/usr/bin/bash

echo mES_ChIP-seq_Ctrl_D0_ChIP-rep1.fastq
bowtie -t -p 16 -v 2 -m 1 --best --strata --sam /public/genomes/mouse/mm10/bowtieIndex_UCSC/mm10 mES_ChIP-seq_Ctrl_D0_ChIP-rep1.fastq mES_ChIP-seq_Ctrl_D0_ChIP-rep1.sam
samtools view -h -bS -F 4 -@ 32 mES_ChIP-seq_Ctrl_D0_ChIP-rep1.sam -o mES_ChIP-seq_Ctrl_D0_ChIP-rep1.bam
rm mES_ChIP-seq_Ctrl_D0_ChIP-rep1.sam
samtools sort --threads 16 -m 2G -O bam -o mES_ChIP-seq_Ctrl_D0_ChIP-rep1.sorted.bam mES_ChIP-seq_Ctrl_D0_ChIP-rep1.bam
rm mES_ChIP-seq_Ctrl_D0_ChIP-rep1.bam
samtools index -b mES_ChIP-seq_Ctrl_D0_ChIP-rep1.sorted.bam

echo mES_ChIP-seq_Ctrl_D0_ChIP-rep2.fastq
bowtie -t -p 16 -v 2 -m 1 --best --strata --sam /public/genomes/mouse/mm10/bowtieIndex_UCSC/mm10 mES_ChIP-seq_Ctrl_D0_ChIP-rep2.fastq mES_ChIP-seq_Ctrl_D0_ChIP-rep2.sam
samtools view -h -bS -F 4 -@ 32 mES_ChIP-seq_Ctrl_D0_ChIP-rep2.sam -o mES_ChIP-seq_Ctrl_D0_ChIP-rep2.bam
rm mES_ChIP-seq_Ctrl_D0_ChIP-rep2.sam
samtools sort --threads 16 -m 2G -O bam -o mES_ChIP-seq_Ctrl_D0_ChIP-rep2.sorted.bam mES_ChIP-seq_Ctrl_D0_ChIP-rep2.bam
rm mES_ChIP-seq_Ctrl_D0_ChIP-rep2.bam
samtools index -b mES_ChIP-seq_Ctrl_D0_ChIP-rep2.sorted.bam

echo mES_ChIP-seq_Ctrl_D0_input-rep1.fastq
bowtie -t -p 16 -v 2 -m 1 --best --strata --sam /public/genomes/mouse/mm10/bowtieIndex_UCSC/mm10 mES_ChIP-seq_Ctrl_D0_input-rep1.fastq mES_ChIP-seq_Ctrl_D0_input-rep1.sam
samtools view -h -bS -F 4 -@ 32 mES_ChIP-seq_Ctrl_D0_input-rep1.sam -o mES_ChIP-seq_Ctrl_D0_input-rep1.bam
rm mES_ChIP-seq_Ctrl_D0_input-rep1.sam
samtools sort --threads 16 -m 2G -O bam -o mES_ChIP-seq_Ctrl_D0_input-rep1.sorted.bam mES_ChIP-seq_Ctrl_D0_input-rep1.bam
rm mES_ChIP-seq_Ctrl_D0_input-rep1.bam
samtools index -b mES_ChIP-seq_Ctrl_D0_input-rep1.sorted.bam

echo mES_ChIP-seq_Ctrl_D0_input-rep2.fastq
bowtie -t -p 16 -v 2 -m 1 --best --strata --sam /public/genomes/mouse/mm10/bowtieIndex_UCSC/mm10 mES_ChIP-seq_Ctrl_D0_input-rep2.fastq mES_ChIP-seq_Ctrl_D0_input-rep2.sam
samtools view -h -bS -F 4 -@ 32 mES_ChIP-seq_Ctrl_D0_input-rep2.sam -o mES_ChIP-seq_Ctrl_D0_input-rep2.bam
rm mES_ChIP-seq_Ctrl_D0_input-rep2.sam
samtools sort --threads 16 -m 2G -O bam -o mES_ChIP-seq_Ctrl_D0_input-rep2.sorted.bam mES_ChIP-seq_Ctrl_D0_input-rep2.bam
rm mES_ChIP-seq_Ctrl_D0_input-rep2.bam
samtools index -b mES_ChIP-seq_Ctrl_D0_input-rep2.sorted.bam

echo mES_ChIP-seq_SetD2-KD_D0_ChIP-rep1.fastq
bowtie -t -p 16 -v 2 -m 1 --best --strata --sam /public/genomes/mouse/mm10/bowtieIndex_UCSC/mm10 mES_ChIP-seq_SetD2-KD_D0_ChIP-rep1.fastq mES_ChIP-seq_SetD2-KD_D0_ChIP-rep1.sam
samtools view -h -bS -F 4 -@ 32 mES_ChIP-seq_SetD2-KD_D0_ChIP-rep1.sam -o mES_ChIP-seq_SetD2-KD_D0_ChIP-rep1.bam
rm mES_ChIP-seq_SetD2-KD_D0_ChIP-rep1.sam
samtools sort --threads 16 -m 2G -O bam -o mES_ChIP-seq_SetD2-KD_D0_ChIP-rep1.sorted.bam mES_ChIP-seq_SetD2-KD_D0_ChIP-rep1.bam
rm mES_ChIP-seq_SetD2-KD_D0_ChIP-rep1.bam
samtools index -b mES_ChIP-seq_SetD2-KD_D0_ChIP-rep1.sorted.bam

echo mES_ChIP-seq_SetD2-KD_D0_ChIP-rep2.fastq
bowtie -t -p 16 -v 2 -m 1 --best --strata --sam /public/genomes/mouse/mm10/bowtieIndex_UCSC/mm10 mES_ChIP-seq_SetD2-KD_D0_ChIP-rep2.fastq mES_ChIP-seq_SetD2-KD_D0_ChIP-rep2.sam
samtools view -h -bS -F 4 -@ 32 mES_ChIP-seq_SetD2-KD_D0_ChIP-rep2.sam -o mES_ChIP-seq_SetD2-KD_D0_ChIP-rep2.bam
rm mES_ChIP-seq_SetD2-KD_D0_ChIP-rep2.sam
samtools sort --threads 16 -m 2G -O bam -o mES_ChIP-seq_SetD2-KD_D0_ChIP-rep2.sorted.bam mES_ChIP-seq_SetD2-KD_D0_ChIP-rep2.bam
rm mES_ChIP-seq_SetD2-KD_D0_ChIP-rep2.bam
samtools index -b mES_ChIP-seq_SetD2-KD_D0_ChIP-rep2.sorted.bam

echo mES_ChIP-seq_SetD2-KD_D0_input-rep1.fastq
bowtie -t -p 16 -v 2 -m 1 --best --strata --sam /public/genomes/mouse/mm10/bowtieIndex_UCSC/mm10 mES_ChIP-seq_SetD2-KD_D0_input-rep1.fastq mES_ChIP-seq_SetD2-KD_D0_input-rep1.sam
samtools view -h -bS -F 4 -@ 32 mES_ChIP-seq_SetD2-KD_D0_input-rep1.sam -o mES_ChIP-seq_SetD2-KD_D0_input-rep1.bam
rm mES_ChIP-seq_SetD2-KD_D0_input-rep1.sam
samtools sort --threads 16 -m 2G -O bam -o mES_ChIP-seq_SetD2-KD_D0_input-rep1.sorted.bam mES_ChIP-seq_SetD2-KD_D0_input-rep1.bam
rm mES_ChIP-seq_SetD2-KD_D0_input-rep1.bam
samtools index -b mES_ChIP-seq_SetD2-KD_D0_input-rep1.sorted.bam

echo mES_ChIP-seq_SetD2-KD_D0_input-rep2.fastq
bowtie -t -p 16 -v 2 -m 1 --best --strata --sam /public/genomes/mouse/mm10/bowtieIndex_UCSC/mm10 mES_ChIP-seq_SetD2-KD_D0_input-rep2.fastq mES_ChIP-seq_SetD2-KD_D0_input-rep2.sam
samtools view -h -bS -F 4 -@ 32 mES_ChIP-seq_SetD2-KD_D0_input-rep2.sam -o mES_ChIP-seq_SetD2-KD_D0_input-rep2.bam
rm mES_ChIP-seq_SetD2-KD_D0_input-rep2.sam
samtools sort --threads 16 -m 2G -O bam -o mES_ChIP-seq_SetD2-KD_D0_input-rep2.sorted.bam mES_ChIP-seq_SetD2-KD_D0_input-rep2.bam
rm mES_ChIP-seq_SetD2-KD_D0_input-rep2.bam
samtools index -b mES_ChIP-seq_SetD2-KD_D0_input-rep2.sorted.bam

exit
