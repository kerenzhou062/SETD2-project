echo HepG2_m6A-seq_shWTAP_input.fastq
tophat -p 10 -g 1 --library-type=fr-firststrand -G /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.gtf -o ./shWTAP/input/ /public/genomes/human/hg19/bowtie2Index_UCSC/hg19 /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/HepG2_m6A-seq_shWTAP_input.fastq
rename ./shWTAP/input/accepted_hits.bam ./shWTAP/input/HepG2_m6A-seq_shWTAP_input.fastq.bam ./shWTAP/input/accepted_hits.bam
mv ./shWTAP/input/HepG2_m6A-seq_shWTAP_input.fastq.bam ./shWTAP/
samtools sort --threads 16 -m 2G -O bam -o ./shWTAP/HepG2_m6A-seq_shWTAP_input.fastq.sorted.bam ./shWTAP/HepG2_m6A-seq_shWTAP_input.fastq.bam
rm ./shWTAP/HepG2_m6A-seq_shWTAP_input.fastq.bam
samtools index -b ./shWTAP/HepG2_m6A-seq_shWTAP_input.fastq.sorted.bam
