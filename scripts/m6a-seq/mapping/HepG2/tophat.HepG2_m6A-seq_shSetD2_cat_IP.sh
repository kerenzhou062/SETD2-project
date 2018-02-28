echo HepG2_m6A-seq_shSetD2_IP.fastq
tophat -p 10 -g 1 --library-type=fr-firststrand -G /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.gtf -o ./shSetD2/IP/ /public/genomes/human/hg19/bowtie2Index_UCSC/hg19 /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/HepG2_m6A-seq_shSetD2_IP.fastq
rename ./shSetD2/IP/accepted_hits.bam ./shSetD2/IP/HepG2_m6A-seq_shSetD2_IP.fastq.bam ./shSetD2/IP/accepted_hits.bam
mv ./shSetD2/IP/HepG2_m6A-seq_shSetD2_IP.fastq.bam ./shSetD2/
samtools sort --threads 16 -m 2G -O bam -o ./shSetD2/HepG2_m6A-seq_shSetD2_IP.fastq.sorted.bam ./shSetD2/HepG2_m6A-seq_shSetD2_IP.fastq.bam
rm ./shSetD2/HepG2_m6A-seq_shSetD2_IP.fastq.bam
samtools index -b ./shSetD2/HepG2_m6A-seq_shSetD2_IP.fastq.sorted.bam
