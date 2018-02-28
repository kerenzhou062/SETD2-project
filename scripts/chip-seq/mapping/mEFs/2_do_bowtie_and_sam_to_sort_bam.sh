echo mEF_ChIP-seq_SetD2-KO_NA_ChIP-rep1.fastq
bowtie -t -p 16 -v 2 -m 1 --best --strata --sam /public/genomes/mouse/mm10/bowtieIndex_UCSC/mm10 mEF_ChIP-seq_SetD2-KO_NA_ChIP-rep1.fastq mEF_ChIP-seq_SetD2-KO_NA_ChIP-rep1.sam
samtools view -h -bS -F 4 --threads 32 mEF_ChIP-seq_SetD2-KO_NA_ChIP-rep1.sam -o mEF_ChIP-seq_SetD2-KO_NA_ChIP-rep1.bam
rm mEF_ChIP-seq_SetD2-KO_NA_ChIP-rep1.sam
samtools sort --threads 16 -m 2G -O bam -o mEF_ChIP-seq_SetD2-KO_NA_ChIP-rep1.sorted.bam mEF_ChIP-seq_SetD2-KO_NA_ChIP-rep1.bam
rm mEF_ChIP-seq_SetD2-KO_NA_ChIP-rep1.bam
samtools index -b mEF_ChIP-seq_SetD2-KO_NA_ChIP-rep1.sorted.bam

echo mEF_ChIP-seq_SetD2-KO_NA_ChIP-rep2.fastq
bowtie -t -p 16 -v 2 -m 1 --best --strata --sam /public/genomes/mouse/mm10/bowtieIndex_UCSC/mm10 mEF_ChIP-seq_SetD2-KO_NA_ChIP-rep2.fastq mEF_ChIP-seq_SetD2-KO_NA_ChIP-rep2.sam
samtools view -h -bS -F 4 --threads 32 mEF_ChIP-seq_SetD2-KO_NA_ChIP-rep2.sam -o mEF_ChIP-seq_SetD2-KO_NA_ChIP-rep2.bam
rm mEF_ChIP-seq_SetD2-KO_NA_ChIP-rep2.sam
samtools sort --threads 16 -m 2G -O bam -o mEF_ChIP-seq_SetD2-KO_NA_ChIP-rep2.sorted.bam mEF_ChIP-seq_SetD2-KO_NA_ChIP-rep2.bam
rm mEF_ChIP-seq_SetD2-KO_NA_ChIP-rep2.bam
samtools index -b mEF_ChIP-seq_SetD2-KO_NA_ChIP-rep2.sorted.bam

echo mEF_ChIP-seq_SetD2-KO_NA_input-rep1.fastq
bowtie -t -p 16 -v 2 -m 1 --best --strata --sam /public/genomes/mouse/mm10/bowtieIndex_UCSC/mm10 mEF_ChIP-seq_SetD2-KO_NA_input-rep1.fastq mEF_ChIP-seq_SetD2-KO_NA_input-rep1.sam
samtools view -h -bS -F 4 --threads 32 mEF_ChIP-seq_SetD2-KO_NA_input-rep1.sam -o mEF_ChIP-seq_SetD2-KO_NA_input-rep1.bam
rm mEF_ChIP-seq_SetD2-KO_NA_input-rep1.sam
samtools sort --threads 16 -m 2G -O bam -o mEF_ChIP-seq_SetD2-KO_NA_input-rep1.sorted.bam mEF_ChIP-seq_SetD2-KO_NA_input-rep1.bam
rm mEF_ChIP-seq_SetD2-KO_NA_input-rep1.bam
samtools index -b mEF_ChIP-seq_SetD2-KO_NA_input-rep1.sorted.bam

echo mEF_ChIP-seq_SetD2-KO_NA_input-rep2.fastq
bowtie -t -p 16 -v 2 -m 1 --best --strata --sam /public/genomes/mouse/mm10/bowtieIndex_UCSC/mm10 mEF_ChIP-seq_SetD2-KO_NA_input-rep2.fastq mEF_ChIP-seq_SetD2-KO_NA_input-rep2.sam
samtools view -h -bS -F 4 --threads 32 mEF_ChIP-seq_SetD2-KO_NA_input-rep2.sam -o mEF_ChIP-seq_SetD2-KO_NA_input-rep2.bam
rm mEF_ChIP-seq_SetD2-KO_NA_input-rep2.sam
samtools sort --threads 16 -m 2G -O bam -o mEF_ChIP-seq_SetD2-KO_NA_input-rep2.sorted.bam mEF_ChIP-seq_SetD2-KO_NA_input-rep2.bam
rm mEF_ChIP-seq_SetD2-KO_NA_input-rep2.bam
samtools index -b mEF_ChIP-seq_SetD2-KO_NA_input-rep2.sorted.bam

echo mEF_ChIP-seq_SetD2-WT_NA_ChIP-rep1.fastq
bowtie -t -p 16 -v 2 -m 1 --best --strata --sam /public/genomes/mouse/mm10/bowtieIndex_UCSC/mm10 mEF_ChIP-seq_SetD2-WT_NA_ChIP-rep1.fastq mEF_ChIP-seq_SetD2-WT_NA_ChIP-rep1.sam
samtools view -h -bS -F 4 --threads 32 mEF_ChIP-seq_SetD2-WT_NA_ChIP-rep1.sam -o mEF_ChIP-seq_SetD2-WT_NA_ChIP-rep1.bam
rm mEF_ChIP-seq_SetD2-WT_NA_ChIP-rep1.sam
samtools sort --threads 16 -m 2G -O bam -o mEF_ChIP-seq_SetD2-WT_NA_ChIP-rep1.sorted.bam mEF_ChIP-seq_SetD2-WT_NA_ChIP-rep1.bam
rm mEF_ChIP-seq_SetD2-WT_NA_ChIP-rep1.bam
samtools index -b mEF_ChIP-seq_SetD2-WT_NA_ChIP-rep1.sorted.bam

echo mEF_ChIP-seq_SetD2-WT_NA_ChIP-rep2.fastq
bowtie -t -p 16 -v 2 -m 1 --best --strata --sam /public/genomes/mouse/mm10/bowtieIndex_UCSC/mm10 mEF_ChIP-seq_SetD2-WT_NA_ChIP-rep2.fastq mEF_ChIP-seq_SetD2-WT_NA_ChIP-rep2.sam
samtools view -h -bS -F 4 --threads 32 mEF_ChIP-seq_SetD2-WT_NA_ChIP-rep2.sam -o mEF_ChIP-seq_SetD2-WT_NA_ChIP-rep2.bam
rm mEF_ChIP-seq_SetD2-WT_NA_ChIP-rep2.sam
samtools sort --threads 16 -m 2G -O bam -o mEF_ChIP-seq_SetD2-WT_NA_ChIP-rep2.sorted.bam mEF_ChIP-seq_SetD2-WT_NA_ChIP-rep2.bam
rm mEF_ChIP-seq_SetD2-WT_NA_ChIP-rep2.bam
samtools index -b mEF_ChIP-seq_SetD2-WT_NA_ChIP-rep2.sorted.bam

echo mEF_ChIP-seq_SetD2-WT_NA_input-rep1.fastq
bowtie -t -p 16 -v 2 -m 1 --best --strata --sam /public/genomes/mouse/mm10/bowtieIndex_UCSC/mm10 mEF_ChIP-seq_SetD2-WT_NA_input-rep1.fastq mEF_ChIP-seq_SetD2-WT_NA_input-rep1.sam
samtools view -h -bS -F 4 --threads 32 mEF_ChIP-seq_SetD2-WT_NA_input-rep1.sam -o mEF_ChIP-seq_SetD2-WT_NA_input-rep1.bam
rm mEF_ChIP-seq_SetD2-WT_NA_input-rep1.sam
samtools sort --threads 16 -m 2G -O bam -o mEF_ChIP-seq_SetD2-WT_NA_input-rep1.sorted.bam mEF_ChIP-seq_SetD2-WT_NA_input-rep1.bam
rm mEF_ChIP-seq_SetD2-WT_NA_input-rep1.bam
samtools index -b mEF_ChIP-seq_SetD2-WT_NA_input-rep1.sorted.bam

echo mEF_ChIP-seq_SetD2-WT_NA_input-rep2.fastq
bowtie -t -p 16 -v 2 -m 1 --best --strata --sam /public/genomes/mouse/mm10/bowtieIndex_UCSC/mm10 mEF_ChIP-seq_SetD2-WT_NA_input-rep2.fastq mEF_ChIP-seq_SetD2-WT_NA_input-rep2.sam
samtools view -h -bS -F 4 --threads 32 mEF_ChIP-seq_SetD2-WT_NA_input-rep2.sam -o mEF_ChIP-seq_SetD2-WT_NA_input-rep2.bam
rm mEF_ChIP-seq_SetD2-WT_NA_input-rep2.sam
samtools sort --threads 16 -m 2G -O bam -o mEF_ChIP-seq_SetD2-WT_NA_input-rep2.sorted.bam mEF_ChIP-seq_SetD2-WT_NA_input-rep2.bam
rm mEF_ChIP-seq_SetD2-WT_NA_input-rep2.bam
samtools index -b mEF_ChIP-seq_SetD2-WT_NA_input-rep2.sorted.bam

