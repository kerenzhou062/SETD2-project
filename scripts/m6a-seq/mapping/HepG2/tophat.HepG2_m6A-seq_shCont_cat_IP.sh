echo HepG2_m6A-seq_shCont_IP.fastq
tophat -p 10 -g 1 --library-type=fr-firststrand -G /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.gtf -o ./shCont/IP/ /public/genomes/human/hg19/bowtie2Index_UCSC/hg19 /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/HepG2_m6A-seq_shCont_IP.fastq
rename ./shCont/IP/accepted_hits.bam ./shCont/IP/HepG2_m6A-seq_shCont_IP.fastq.bam ./shCont/IP/accepted_hits.bam
mv ./shCont/IP/HepG2_m6A-seq_shCont_IP.fastq.bam ./shCont/
samtools sort --threads 16 -m 2G -O bam -o ./shCont/HepG2_m6A-seq_shCont_IP.fastq.sorted.bam ./shCont/HepG2_m6A-seq_shCont_IP.fastq.bam
rm ./shCont/HepG2_m6A-seq_shCont_IP.fastq.bam
samtools index -b ./shCont/HepG2_m6A-seq_shCont_IP.fastq.sorted.bam
