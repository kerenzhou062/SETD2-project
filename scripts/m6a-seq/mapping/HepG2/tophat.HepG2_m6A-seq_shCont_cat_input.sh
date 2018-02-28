echo HepG2_m6A-seq_shCont_input.fastq
tophat -p 10 -g 1 --library-type=fr-firststrand -G /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.gtf -o ./shCont/input/ /public/genomes/human/hg19/bowtie2Index_UCSC/hg19 /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/HepG2_m6A-seq_shCont_input.fastq
rename ./shCont/input/accepted_hits.bam ./shCont/input/HepG2_m6A-seq_shCont_input.fastq.bam ./shCont/input/accepted_hits.bam
mv ./shCont/input/HepG2_m6A-seq_shCont_input.fastq.bam ./shCont/
samtools sort --threads 16 -m 2G -O bam -o ./shCont/HepG2_m6A-seq_shCont_input.fastq.sorted.bam ./shCont/HepG2_m6A-seq_shCont_input.fastq.bam
rm ./shCont/HepG2_m6A-seq_shCont_input.fastq.bam
samtools index -b ./shCont/HepG2_m6A-seq_shCont_input.fastq.sorted.bam
