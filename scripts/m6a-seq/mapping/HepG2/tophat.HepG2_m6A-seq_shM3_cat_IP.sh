echo HepG2_m6A-seq_shM3_IP.fastq
tophat -p 10 -g 1 --library-type=fr-firststrand -G /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.gtf -o ./shM3/IP/ /public/genomes/human/hg19/bowtie2Index_UCSC/hg19 /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/HepG2_m6A-seq_shM3_IP.fastq
rename ./shM3/IP/accepted_hits.bam ./shM3/IP/HepG2_m6A-seq_shM3_IP.fastq.bam ./shM3/IP/accepted_hits.bam
mv ./shM3/IP/HepG2_m6A-seq_shM3_IP.fastq.bam ./shM3/
samtools sort --threads 16 -m 2G -O bam -o ./shM3/HepG2_m6A-seq_shM3_IP.fastq.sorted.bam ./shM3/HepG2_m6A-seq_shM3_IP.fastq.bam
rm ./shM3/HepG2_m6A-seq_shM3_IP.fastq.bam
samtools index -b ./shM3/HepG2_m6A-seq_shM3_IP.fastq.sorted.bam
