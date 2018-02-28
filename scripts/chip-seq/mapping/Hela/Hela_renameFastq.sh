
mv CHe-BZ-JJ-1_S1_L007_R1_001.fastq.gz Hela_ChIP-seq_shCont_input_rep2.fastq.gz
mv CHe-BZ-JJ-2_S2_L007_R1_001.fastq.gz Hela_ChIP-seq_shSetD2_input_rep2.fastq.gz
mv CHe-BZ-JJ-5_S5_L007_R1_001.fastq.gz Hela_ChIP-seq_shCont_IP_rep2.fastq.gz
mv CHe-BZ-JJ-6_S6_L007_R1_001.fastq.gz Hela_ChIP-seq_shSetD2_IP_rep2.fastq.gz

gzip -dq Hela_ChIP-seq_shCont_input_rep2.fastq.gz &
gzip -dq Hela_ChIP-seq_shSetD2_input_rep2.fastq.gz &
gzip -dq Hela_ChIP-seq_shCont_IP_rep2.fastq.gz &
gzip -dq Hela_ChIP-seq_shSetD2_IP_rep2.fastq.gz &

