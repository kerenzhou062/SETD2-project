echo "\n\n\n"
echo Hela_m6A-seq_shM14_input_rep1.fastq
tophat -p 2 -g 1 --library-type=fr-firststrand -G /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.gtf -o /data/zhoukr/hhl_setd2_m6a/Hela_m6A-seq/shM14/input_rep1/ /public/genomes/human/hg19/bowtie2Index_UCSC/hg19 /data/zhoukr/hhl_setd2_m6a/Hela_m6A-seq/Hela_m6A-seq_shM14_input_rep1.fastq
rename /data/zhoukr/hhl_setd2_m6a/Hela_m6A-seq/shM14/input_rep1/accepted_hits.bam /data/zhoukr/hhl_setd2_m6a/Hela_m6A-seq/shM14/input_rep1/Hela_m6A-seq_shM14_input_rep1.fastq.bam /data/zhoukr/hhl_setd2_m6a/Hela_m6A-seq/shM14/input_rep1/accepted_hits.bam
mv /data/zhoukr/hhl_setd2_m6a/Hela_m6A-seq/shM14/input_rep1/Hela_m6A-seq_shM14_input_rep1.fastq.bam /data/zhoukr/hhl_setd2_m6a/Hela_m6A-seq/shM14/
samtools sort --threads 16 -m 2G -O bam -o /data/zhoukr/hhl_setd2_m6a/Hela_m6A-seq/shM14/Hela_m6A-seq_shM14_input_rep1.fastq.sorted.bam /data/zhoukr/hhl_setd2_m6a/Hela_m6A-seq/shM14/Hela_m6A-seq_shM14_input_rep1.fastq.bam
rm /data/zhoukr/hhl_setd2_m6a/Hela_m6A-seq/shM14/Hela_m6A-seq_shM14_input_rep1.fastq.bam
samtools index -b /data/zhoukr/hhl_setd2_m6a/Hela_m6A-seq/shM14/Hela_m6A-seq_shM14_input_rep1.fastq.sorted.bam

