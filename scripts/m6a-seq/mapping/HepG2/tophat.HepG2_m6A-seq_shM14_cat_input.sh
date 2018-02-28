echo HepG2_m6A-seq_shM14_input.fastq
tophat -p 10 -g 1 --library-type=fr-firststrand -G /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.gtf -o ./shM14/input/ /public/genomes/human/hg19/bowtie2Index_UCSC/hg19 /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/HepG2_m6A-seq_shM14_input.fastq
rename ./shM14/input/accepted_hits.bam ./shM14/input/HepG2_m6A-seq_shM14_input.fastq.bam ./shM14/input/accepted_hits.bam
mv ./shM14/input/HepG2_m6A-seq_shM14_input.fastq.bam ./shM14/
samtools sort --threads 16 -m 2G -O bam -o ./shM14/HepG2_m6A-seq_shM14_input.fastq.sorted.bam ./shM14/HepG2_m6A-seq_shM14_input.fastq.bam
rm ./shM14/HepG2_m6A-seq_shM14_input.fastq.bam
samtools index -b ./shM14/HepG2_m6A-seq_shM14_input.fastq.sorted.bam
