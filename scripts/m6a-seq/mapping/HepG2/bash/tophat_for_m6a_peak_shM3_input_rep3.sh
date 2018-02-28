echo "\n\n\n"
echo HepG2_m6A-seq_shM3_input_rep3.fastq
tophat -p 2 -g 1 --library-type=fr-firststrand -G /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.gtf -o /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/shM3/input_rep3/ /public/genomes/human/hg19/bowtie2Index_UCSC/hg19 /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/HepG2_m6A-seq_shM3_input_rep3.fastq
rename /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/shM3/input_rep3/accepted_hits.bam /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/shM3/input_rep3/HepG2_m6A-seq_shM3_input_rep3.fastq.bam /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/shM3/input_rep3/accepted_hits.bam
mv /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/shM3/input_rep3/HepG2_m6A-seq_shM3_input_rep3.fastq.bam /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/shM3/
samtools sort --threads 16 -m 2G -O bam -o /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/shM3/HepG2_m6A-seq_shM3_input_rep3.fastq.sorted.bam /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/shM3/HepG2_m6A-seq_shM3_input_rep3.fastq.bam
rm /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/shM3/HepG2_m6A-seq_shM3_input_rep3.fastq.bam
samtools index -b /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/shM3/HepG2_m6A-seq_shM3_input_rep3.fastq.sorted.bam

