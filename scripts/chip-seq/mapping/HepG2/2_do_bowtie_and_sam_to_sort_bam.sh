echo HepG2_ChIP-seq_shCont_input.fastq
bowtie -t -p 16 -v 2 -m 1 --best --strata --sam /data/sunwj/genomeInfo/human/hg19/bowtieIndex_UCSC/hg19 HepG2_ChIP-seq_shCont_input.fastq HepG2_ChIP-seq_shCont_input.sam
samtools view -h -bS -F 4 --threads 32 HepG2_ChIP-seq_shCont_input.sam -o HepG2_ChIP-seq_shCont_input.bam
rm HepG2_ChIP-seq_shCont_input.sam
samtools sort --threads 16 -m 2G -O bam -o HepG2_ChIP-seq_shCont_input.sorted.bam HepG2_ChIP-seq_shCont_input.bam
rm HepG2_ChIP-seq_shCont_input.bam
samtools index -b HepG2_ChIP-seq_shCont_input.sorted.bam

echo HepG2_ChIP-seq_shCont_IP.fastq
bowtie -t -p 16 -v 2 -m 1 --best --strata --sam /data/sunwj/genomeInfo/human/hg19/bowtieIndex_UCSC/hg19 HepG2_ChIP-seq_shCont_IP.fastq HepG2_ChIP-seq_shCont_IP.sam
samtools view -h -bS -F 4 --threads 32 HepG2_ChIP-seq_shCont_IP.sam -o HepG2_ChIP-seq_shCont_IP.bam
rm HepG2_ChIP-seq_shCont_IP.sam
samtools sort --threads 16 -m 2G -O bam -o HepG2_ChIP-seq_shCont_IP.sorted.bam HepG2_ChIP-seq_shCont_IP.bam
rm HepG2_ChIP-seq_shCont_IP.bam
samtools index -b HepG2_ChIP-seq_shCont_IP.sorted.bam

echo HepG2_ChIP-seq_shM14_input.fastq
bowtie -t -p 16 -v 2 -m 1 --best --strata --sam /data/sunwj/genomeInfo/human/hg19/bowtieIndex_UCSC/hg19 HepG2_ChIP-seq_shM14_input.fastq HepG2_ChIP-seq_shM14_input.sam
samtools view -h -bS -F 4 --threads 32 HepG2_ChIP-seq_shM14_input.sam -o HepG2_ChIP-seq_shM14_input.bam
rm HepG2_ChIP-seq_shM14_input.sam
samtools sort --threads 16 -m 2G -O bam -o HepG2_ChIP-seq_shM14_input.sorted.bam HepG2_ChIP-seq_shM14_input.bam
rm HepG2_ChIP-seq_shM14_input.bam
samtools index -b HepG2_ChIP-seq_shM14_input.sorted.bam

echo HepG2_ChIP-seq_shM14_IP.fastq
bowtie -t -p 16 -v 2 -m 1 --best --strata --sam /data/sunwj/genomeInfo/human/hg19/bowtieIndex_UCSC/hg19 HepG2_ChIP-seq_shM14_IP.fastq HepG2_ChIP-seq_shM14_IP.sam
samtools view -h -bS -F 4 --threads 32 HepG2_ChIP-seq_shM14_IP.sam -o HepG2_ChIP-seq_shM14_IP.bam
rm HepG2_ChIP-seq_shM14_IP.sam
samtools sort --threads 16 -m 2G -O bam -o HepG2_ChIP-seq_shM14_IP.sorted.bam HepG2_ChIP-seq_shM14_IP.bam
rm HepG2_ChIP-seq_shM14_IP.bam
samtools index -b HepG2_ChIP-seq_shM14_IP.sorted.bam

echo HepG2_ChIP-seq_shM3_input.fastq
bowtie -t -p 16 -v 2 -m 1 --best --strata --sam /data/sunwj/genomeInfo/human/hg19/bowtieIndex_UCSC/hg19 HepG2_ChIP-seq_shM3_input.fastq HepG2_ChIP-seq_shM3_input.sam
samtools view -h -bS -F 4 --threads 32 HepG2_ChIP-seq_shM3_input.sam -o HepG2_ChIP-seq_shM3_input.bam
rm HepG2_ChIP-seq_shM3_input.sam
samtools sort --threads 16 -m 2G -O bam -o HepG2_ChIP-seq_shM3_input.sorted.bam HepG2_ChIP-seq_shM3_input.bam
rm HepG2_ChIP-seq_shM3_input.bam
samtools index -b HepG2_ChIP-seq_shM3_input.sorted.bam

echo HepG2_ChIP-seq_shM3_IP.fastq
bowtie -t -p 16 -v 2 -m 1 --best --strata --sam /data/sunwj/genomeInfo/human/hg19/bowtieIndex_UCSC/hg19 HepG2_ChIP-seq_shM3_IP.fastq HepG2_ChIP-seq_shM3_IP.sam
samtools view -h -bS -F 4 --threads 32 HepG2_ChIP-seq_shM3_IP.sam -o HepG2_ChIP-seq_shM3_IP.bam
rm HepG2_ChIP-seq_shM3_IP.sam
samtools sort --threads 16 -m 2G -O bam -o HepG2_ChIP-seq_shM3_IP.sorted.bam HepG2_ChIP-seq_shM3_IP.bam
rm HepG2_ChIP-seq_shM3_IP.bam
samtools index -b HepG2_ChIP-seq_shM3_IP.sorted.bam

echo HepG2_ChIP-seq_shSetD2_input.fastq
bowtie -t -p 16 -v 2 -m 1 --best --strata --sam /data/sunwj/genomeInfo/human/hg19/bowtieIndex_UCSC/hg19 HepG2_ChIP-seq_shSetD2_input.fastq HepG2_ChIP-seq_shSetD2_input.sam
samtools view -h -bS -F 4 --threads 32 HepG2_ChIP-seq_shSetD2_input.sam -o HepG2_ChIP-seq_shSetD2_input.bam
rm HepG2_ChIP-seq_shSetD2_input.sam
samtools sort --threads 16 -m 2G -O bam -o HepG2_ChIP-seq_shSetD2_input.sorted.bam HepG2_ChIP-seq_shSetD2_input.bam
rm HepG2_ChIP-seq_shSetD2_input.bam
samtools index -b HepG2_ChIP-seq_shSetD2_input.sorted.bam

echo HepG2_ChIP-seq_shSetD2_IP.fastq
bowtie -t -p 16 -v 2 -m 1 --best --strata --sam /data/sunwj/genomeInfo/human/hg19/bowtieIndex_UCSC/hg19 HepG2_ChIP-seq_shSetD2_IP.fastq HepG2_ChIP-seq_shSetD2_IP.sam
samtools view -h -bS -F 4 --threads 32 HepG2_ChIP-seq_shSetD2_IP.sam -o HepG2_ChIP-seq_shSetD2_IP.bam
rm HepG2_ChIP-seq_shSetD2_IP.sam
samtools sort --threads 16 -m 2G -O bam -o HepG2_ChIP-seq_shSetD2_IP.sorted.bam HepG2_ChIP-seq_shSetD2_IP.bam
rm HepG2_ChIP-seq_shSetD2_IP.bam
samtools index -b HepG2_ChIP-seq_shSetD2_IP.sorted.bam

echo HepG2_ChIP-seq_shWTAP_input.fastq
bowtie -t -p 16 -v 2 -m 1 --best --strata --sam /data/sunwj/genomeInfo/human/hg19/bowtieIndex_UCSC/hg19 HepG2_ChIP-seq_shWTAP_input.fastq HepG2_ChIP-seq_shWTAP_input.sam
samtools view -h -bS -F 4 --threads 32 HepG2_ChIP-seq_shWTAP_input.sam -o HepG2_ChIP-seq_shWTAP_input.bam
rm HepG2_ChIP-seq_shWTAP_input.sam
samtools sort --threads 16 -m 2G -O bam -o HepG2_ChIP-seq_shWTAP_input.sorted.bam HepG2_ChIP-seq_shWTAP_input.bam
rm HepG2_ChIP-seq_shWTAP_input.bam
samtools index -b HepG2_ChIP-seq_shWTAP_input.sorted.bam

echo HepG2_ChIP-seq_shWTAP_IP.fastq
bowtie -t -p 16 -v 2 -m 1 --best --strata --sam /data/sunwj/genomeInfo/human/hg19/bowtieIndex_UCSC/hg19 HepG2_ChIP-seq_shWTAP_IP.fastq HepG2_ChIP-seq_shWTAP_IP.sam
samtools view -h -bS -F 4 --threads 32 HepG2_ChIP-seq_shWTAP_IP.sam -o HepG2_ChIP-seq_shWTAP_IP.bam
rm HepG2_ChIP-seq_shWTAP_IP.sam
samtools sort --threads 16 -m 2G -O bam -o HepG2_ChIP-seq_shWTAP_IP.sorted.bam HepG2_ChIP-seq_shWTAP_IP.bam
rm HepG2_ChIP-seq_shWTAP_IP.bam
samtools index -b HepG2_ChIP-seq_shWTAP_IP.sorted.bam

