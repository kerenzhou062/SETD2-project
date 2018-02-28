echo "\n\n\n"
echo Hela_m6A-seq_shCont_input_rep2.fastq
tophat -p 2 -g 1 --library-type=fr-firststrand -G /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.gtf -o /data/zhoukr/hhl_setd2_m6a/Hela_m6A-seq/shCont/input_rep2/ /public/genomes/human/hg19/bowtie2Index_UCSC/hg19 /data/zhoukr/hhl_setd2_m6a/Hela_m6A-seq/Hela_m6A-seq_shCont_input_rep2.fastq
rename /data/zhoukr/hhl_setd2_m6a/Hela_m6A-seq/shCont/input_rep2/accepted_hits.bam /data/zhoukr/hhl_setd2_m6a/Hela_m6A-seq/shCont/input_rep2/Hela_m6A-seq_shCont_input_rep2.fastq.bam /data/zhoukr/hhl_setd2_m6a/Hela_m6A-seq/shCont/input_rep2/accepted_hits.bam
mv /data/zhoukr/hhl_setd2_m6a/Hela_m6A-seq/shCont/input_rep2/Hela_m6A-seq_shCont_input_rep2.fastq.bam /data/zhoukr/hhl_setd2_m6a/Hela_m6A-seq/shCont/
samtools sort --threads 16 -m 2G -O bam -o /data/zhoukr/hhl_setd2_m6a/Hela_m6A-seq/shCont/Hela_m6A-seq_shCont_input_rep2.fastq.sorted.bam /data/zhoukr/hhl_setd2_m6a/Hela_m6A-seq/shCont/Hela_m6A-seq_shCont_input_rep2.fastq.bam
rm /data/zhoukr/hhl_setd2_m6a/Hela_m6A-seq/shCont/Hela_m6A-seq_shCont_input_rep2.fastq.bam
samtools index -b /data/zhoukr/hhl_setd2_m6a/Hela_m6A-seq/shCont/Hela_m6A-seq_shCont_input_rep2.fastq.sorted.bam

