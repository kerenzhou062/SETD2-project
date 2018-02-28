
mv CHe-BZ-JJ-3_S3_L007_R1_001.fastq.gz HepG2_ChIP-seq_shCont_input_rep2.fastq.gz
mv CHe-BZ-JJ-4_S4_L007_R1_001.fastq.gz HepG2_ChIP-seq_shSetD2_input_rep2.fastq.gz
mv CHe-BZ-JJ-7_S7_L007_R1_001.fastq.gz HepG2_ChIP-seq_shCont_IP_rep2.fastq.gz
mv CHe-BZ-JJ-8_S8_L007_R1_001.fastq.gz HepG2_ChIP-seq_shSetD2_IP_rep2.fastq.gz

gzip -dq HepG2_ChIP-seq_shCont_input_rep2.fastq.gz &
gzip -dq HepG2_ChIP-seq_shSetD2_input_rep2.fastq.gz &
gzip -dq HepG2_ChIP-seq_shCont_IP_rep2.fastq.gz &
gzip -dq HepG2_ChIP-seq_shSetD2_IP_rep2.fastq.gz &
