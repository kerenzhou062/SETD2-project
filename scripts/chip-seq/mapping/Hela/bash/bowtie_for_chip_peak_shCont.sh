echo -e "Hela_ChIP-seq_shCont_input.fastq"
bowtie -t -p 16 -v 2 -m 1 --best --strata --sam /public/genomes/human/hg19/bowtieIndex_UCSC/hg19 /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/Hela_ChIP-seq_shCont_input.fastq /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shCont/input/Hela_ChIP-seq_shCont_input.fastq.sam
samtools view -h -bS -F 4 --threads 32 /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shCont/input/Hela_ChIP-seq_shCont_input.fastq.sam -o /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shCont/input/Hela_ChIP-seq_shCont_input.fastq.bam
rm /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shCont/input/Hela_ChIP-seq_shCont_input.fastq.sam
samtools sort --threads 16 -m 2G -O bam -o /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shCont/input/Hela_ChIP-seq_shCont_input.fastq.sorted.bam /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shCont/input/Hela_ChIP-seq_shCont_input.fastq.bam
rm /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shCont/input/Hela_ChIP-seq_shCont_input.fastq.bam
samtools index -b /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shCont/input/Hela_ChIP-seq_shCont_input.fastq.sorted.bam

echo -e "Hela_ChIP-seq_shCont_IP.fastq"
bowtie -t -p 16 -v 2 -m 1 --best --strata --sam /public/genomes/human/hg19/bowtieIndex_UCSC/hg19 /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/Hela_ChIP-seq_shCont_IP.fastq /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shCont/IP/Hela_ChIP-seq_shCont_IP.fastq.sam
samtools view -h -bS -F 4 --threads 32 /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shCont/IP/Hela_ChIP-seq_shCont_IP.fastq.sam -o /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shCont/IP/Hela_ChIP-seq_shCont_IP.fastq.bam
rm /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shCont/IP/Hela_ChIP-seq_shCont_IP.fastq.sam
samtools sort --threads 16 -m 2G -O bam -o /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shCont/IP/Hela_ChIP-seq_shCont_IP.fastq.sorted.bam /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shCont/IP/Hela_ChIP-seq_shCont_IP.fastq.bam
rm /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shCont/IP/Hela_ChIP-seq_shCont_IP.fastq.bam
samtools index -b /data/zhoukr/hhl_setd2_m6a/Hela_ChIP-seq/shCont/IP/Hela_ChIP-seq_shCont_IP.fastq.sorted.bam


