echo mEF_m6A-seq_SetD2-KO_NA_input.fastq
tophat -p 32 -g 1 --library-type=fr-firststrand -G /data/sunwj/genomeInfo/mouse/mm10/annotation/gencode.vM9.annotation.gtf -o ./SetD2-KO_NA/input/ /public/genomes/mouse/mm10/bowtie2Index_UCSC/mm10 mEF_m6A-seq_SetD2-KO_NA_input.fastq
rename ./SetD2-KO_NA/input/accepted_hits.bam ./SetD2-KO_NA/input/mEF_m6A-seq_SetD2-KO_NA_input.bam ./SetD2-KO_NA/input/accepted_hits.bam
mv ./SetD2-KO_NA/input/mEF_m6A-seq_SetD2-KO_NA_input.bam ./SetD2-KO_NA/
samtools sort --threads 16 -m 2G -O bam -o ./SetD2-KO_NA/mEF_m6A-seq_SetD2-KO_NA_input.sorted.bam ./SetD2-KO_NA/mEF_m6A-seq_SetD2-KO_NA_input.bam
rm ./SetD2-KO_NA/mEF_m6A-seq_SetD2-KO_NA_input.bam
samtools index -b ./SetD2-KO_NA/mEF_m6A-seq_SetD2-KO_NA_input.sorted.bam

echo mEF_m6A-seq_SetD2-KO_NA_IP.fastq
tophat -p 32 -g 1 --library-type=fr-firststrand -G /data/sunwj/genomeInfo/mouse/mm10/annotation/gencode.vM9.annotation.gtf -o ./SetD2-KO_NA/IP/ /public/genomes/mouse/mm10/bowtie2Index_UCSC/mm10 mEF_m6A-seq_SetD2-KO_NA_IP.fastq
rename ./SetD2-KO_NA/IP/accepted_hits.bam ./SetD2-KO_NA/IP/mEF_m6A-seq_SetD2-KO_NA_IP.bam ./SetD2-KO_NA/IP/accepted_hits.bam
mv ./SetD2-KO_NA/IP/mEF_m6A-seq_SetD2-KO_NA_IP.bam ./SetD2-KO_NA/
samtools sort --threads 16 -m 2G -O bam -o ./SetD2-KO_NA/mEF_m6A-seq_SetD2-KO_NA_IP.sorted.bam ./SetD2-KO_NA/mEF_m6A-seq_SetD2-KO_NA_IP.bam
rm ./SetD2-KO_NA/mEF_m6A-seq_SetD2-KO_NA_IP.bam
samtools index -b ./SetD2-KO_NA/mEF_m6A-seq_SetD2-KO_NA_IP.sorted.bam

echo mEF_m6A-seq_SetD2-WT_NA_input.fastq
tophat -p 32 -g 1 --library-type=fr-firststrand -G /data/sunwj/genomeInfo/mouse/mm10/annotation/gencode.vM9.annotation.gtf -o ./SetD2-WT_NA/input/ /public/genomes/mouse/mm10/bowtie2Index_UCSC/mm10 mEF_m6A-seq_SetD2-WT_NA_input.fastq
rename ./SetD2-WT_NA/input/accepted_hits.bam ./SetD2-WT_NA/input/mEF_m6A-seq_SetD2-WT_NA_input.bam ./SetD2-WT_NA/input/accepted_hits.bam
mv ./SetD2-WT_NA/input/mEF_m6A-seq_SetD2-WT_NA_input.bam ./SetD2-WT_NA/
samtools sort --threads 16 -m 2G -O bam -o ./SetD2-WT_NA/mEF_m6A-seq_SetD2-WT_NA_input.sorted.bam ./SetD2-WT_NA/mEF_m6A-seq_SetD2-WT_NA_input.bam
rm ./SetD2-WT_NA/mEF_m6A-seq_SetD2-WT_NA_input.bam
samtools index -b ./SetD2-WT_NA/mEF_m6A-seq_SetD2-WT_NA_input.sorted.bam

echo mEF_m6A-seq_SetD2-WT_NA_IP.fastq
tophat -p 32 -g 1 --library-type=fr-firststrand -G /data/sunwj/genomeInfo/mouse/mm10/annotation/gencode.vM9.annotation.gtf -o ./SetD2-WT_NA/IP/ /public/genomes/mouse/mm10/bowtie2Index_UCSC/mm10 mEF_m6A-seq_SetD2-WT_NA_IP.fastq
rename ./SetD2-WT_NA/IP/accepted_hits.bam ./SetD2-WT_NA/IP/mEF_m6A-seq_SetD2-WT_NA_IP.bam ./SetD2-WT_NA/IP/accepted_hits.bam
mv ./SetD2-WT_NA/IP/mEF_m6A-seq_SetD2-WT_NA_IP.bam ./SetD2-WT_NA/
samtools sort --threads 16 -m 2G -O bam -o ./SetD2-WT_NA/mEF_m6A-seq_SetD2-WT_NA_IP.sorted.bam ./SetD2-WT_NA/mEF_m6A-seq_SetD2-WT_NA_IP.bam
rm ./SetD2-WT_NA/mEF_m6A-seq_SetD2-WT_NA_IP.bam
samtools index -b ./SetD2-WT_NA/mEF_m6A-seq_SetD2-WT_NA_IP.sorted.bam

