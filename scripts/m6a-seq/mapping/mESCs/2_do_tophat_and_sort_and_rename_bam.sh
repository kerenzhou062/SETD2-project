echo mES_m6A-seq_Ctrl_D0_input.fastq
tophat -p 32 -g 1 --library-type=fr-firststrand -G /data/sunwj/genomeInfo/mouse/mm10/annotation/gencode.vM9.annotation.gtf -o ./Ctrl_D0/input/ /public/genomes/mouse/mm10/bowtie2Index_UCSC/mm10 mES_m6A-seq_Ctrl_D0_input.fastq
rename ./Ctrl_D0/input/accepted_hits.bam ./Ctrl_D0/input/mES_m6A-seq_Ctrl_D0_input.bam ./Ctrl_D0/input/accepted_hits.bam
mv ./Ctrl_D0/input/mES_m6A-seq_Ctrl_D0_input.bam ./Ctrl_D0/
samtools sort --threads 16 -m 2G -O bam -o ./Ctrl_D0/mES_m6A-seq_Ctrl_D0_input.sorted.bam ./Ctrl_D0/mES_m6A-seq_Ctrl_D0_input.bam
rm ./Ctrl_D0/mES_m6A-seq_Ctrl_D0_input.bam
samtools index -b ./Ctrl_D0/mES_m6A-seq_Ctrl_D0_input.sorted.bam

echo mES_m6A-seq_Ctrl_D0_IP.fastq
tophat -p 32 -g 1 --library-type=fr-firststrand -G /data/sunwj/genomeInfo/mouse/mm10/annotation/gencode.vM9.annotation.gtf -o ./Ctrl_D0/IP/ /public/genomes/mouse/mm10/bowtie2Index_UCSC/mm10 mES_m6A-seq_Ctrl_D0_IP.fastq
rename ./Ctrl_D0/IP/accepted_hits.bam ./Ctrl_D0/IP/mES_m6A-seq_Ctrl_D0_IP.bam ./Ctrl_D0/IP/accepted_hits.bam
mv ./Ctrl_D0/IP/mES_m6A-seq_Ctrl_D0_IP.bam ./Ctrl_D0/
samtools sort --threads 16 -m 2G -O bam -o ./Ctrl_D0/mES_m6A-seq_Ctrl_D0_IP.sorted.bam ./Ctrl_D0/mES_m6A-seq_Ctrl_D0_IP.bam
rm ./Ctrl_D0/mES_m6A-seq_Ctrl_D0_IP.bam
samtools index -b ./Ctrl_D0/mES_m6A-seq_Ctrl_D0_IP.sorted.bam

echo mES_m6A-seq_Ctrl_D6_input.fastq
tophat -p 32 -g 1 --library-type=fr-firststrand -G /data/sunwj/genomeInfo/mouse/mm10/annotation/gencode.vM9.annotation.gtf -o ./Ctrl_D6/input/ /public/genomes/mouse/mm10/bowtie2Index_UCSC/mm10 mES_m6A-seq_Ctrl_D6_input.fastq
rename ./Ctrl_D6/input/accepted_hits.bam ./Ctrl_D6/input/mES_m6A-seq_Ctrl_D6_input.bam ./Ctrl_D6/input/accepted_hits.bam
mv ./Ctrl_D6/input/mES_m6A-seq_Ctrl_D6_input.bam ./Ctrl_D6/
samtools sort --threads 16 -m 2G -O bam -o ./Ctrl_D6/mES_m6A-seq_Ctrl_D6_input.sorted.bam ./Ctrl_D6/mES_m6A-seq_Ctrl_D6_input.bam
rm ./Ctrl_D6/mES_m6A-seq_Ctrl_D6_input.bam
samtools index -b ./Ctrl_D6/mES_m6A-seq_Ctrl_D6_input.sorted.bam

echo mES_m6A-seq_Ctrl_D6_IP.fastq
tophat -p 32 -g 1 --library-type=fr-firststrand -G /data/sunwj/genomeInfo/mouse/mm10/annotation/gencode.vM9.annotation.gtf -o ./Ctrl_D6/IP/ /public/genomes/mouse/mm10/bowtie2Index_UCSC/mm10 mES_m6A-seq_Ctrl_D6_IP.fastq
rename ./Ctrl_D6/IP/accepted_hits.bam ./Ctrl_D6/IP/mES_m6A-seq_Ctrl_D6_IP.bam ./Ctrl_D6/IP/accepted_hits.bam
mv ./Ctrl_D6/IP/mES_m6A-seq_Ctrl_D6_IP.bam ./Ctrl_D6/
samtools sort --threads 16 -m 2G -O bam -o ./Ctrl_D6/mES_m6A-seq_Ctrl_D6_IP.sorted.bam ./Ctrl_D6/mES_m6A-seq_Ctrl_D6_IP.bam
rm ./Ctrl_D6/mES_m6A-seq_Ctrl_D6_IP.bam
samtools index -b ./Ctrl_D6/mES_m6A-seq_Ctrl_D6_IP.sorted.bam

echo mES_m6A-seq_SetD2-KD_D0_input.fastq
tophat -p 32 -g 1 --library-type=fr-firststrand -G /data/sunwj/genomeInfo/mouse/mm10/annotation/gencode.vM9.annotation.gtf -o ./SetD2-KD_D0/input/ /public/genomes/mouse/mm10/bowtie2Index_UCSC/mm10 mES_m6A-seq_SetD2-KD_D0_input.fastq
rename ./SetD2-KD_D0/input/accepted_hits.bam ./SetD2-KD_D0/input/mES_m6A-seq_SetD2-KD_D0_input.bam ./SetD2-KD_D0/input/accepted_hits.bam
mv ./SetD2-KD_D0/input/mES_m6A-seq_SetD2-KD_D0_input.bam ./SetD2-KD_D0/
samtools sort --threads 16 -m 2G -O bam -o ./SetD2-KD_D0/mES_m6A-seq_SetD2-KD_D0_input.sorted.bam ./SetD2-KD_D0/mES_m6A-seq_SetD2-KD_D0_input.bam
rm ./SetD2-KD_D0/mES_m6A-seq_SetD2-KD_D0_input.bam
samtools index -b ./SetD2-KD_D0/mES_m6A-seq_SetD2-KD_D0_input.sorted.bam

echo mES_m6A-seq_SetD2-KD_D0_IP.fastq
tophat -p 32 -g 1 --library-type=fr-firststrand -G /data/sunwj/genomeInfo/mouse/mm10/annotation/gencode.vM9.annotation.gtf -o ./SetD2-KD_D0/IP/ /public/genomes/mouse/mm10/bowtie2Index_UCSC/mm10 mES_m6A-seq_SetD2-KD_D0_IP.fastq
rename ./SetD2-KD_D0/IP/accepted_hits.bam ./SetD2-KD_D0/IP/mES_m6A-seq_SetD2-KD_D0_IP.bam ./SetD2-KD_D0/IP/accepted_hits.bam
mv ./SetD2-KD_D0/IP/mES_m6A-seq_SetD2-KD_D0_IP.bam ./SetD2-KD_D0/
samtools sort --threads 16 -m 2G -O bam -o ./SetD2-KD_D0/mES_m6A-seq_SetD2-KD_D0_IP.sorted.bam ./SetD2-KD_D0/mES_m6A-seq_SetD2-KD_D0_IP.bam
rm ./SetD2-KD_D0/mES_m6A-seq_SetD2-KD_D0_IP.bam
samtools index -b ./SetD2-KD_D0/mES_m6A-seq_SetD2-KD_D0_IP.sorted.bam

echo mES_m6A-seq_SetD2-KD_D6_input.fastq
tophat -p 32 -g 1 --library-type=fr-firststrand -G /data/sunwj/genomeInfo/mouse/mm10/annotation/gencode.vM9.annotation.gtf -o ./SetD2-KD_D6/input/ /public/genomes/mouse/mm10/bowtie2Index_UCSC/mm10 mES_m6A-seq_SetD2-KD_D6_input.fastq
rename ./SetD2-KD_D6/input/accepted_hits.bam ./SetD2-KD_D6/input/mES_m6A-seq_SetD2-KD_D6_input.bam ./SetD2-KD_D6/input/accepted_hits.bam
mv ./SetD2-KD_D6/input/mES_m6A-seq_SetD2-KD_D6_input.bam ./SetD2-KD_D6/
samtools sort --threads 16 -m 2G -O bam -o ./SetD2-KD_D6/mES_m6A-seq_SetD2-KD_D6_input.sorted.bam ./SetD2-KD_D6/mES_m6A-seq_SetD2-KD_D6_input.bam
rm ./SetD2-KD_D6/mES_m6A-seq_SetD2-KD_D6_input.bam
samtools index -b ./SetD2-KD_D6/mES_m6A-seq_SetD2-KD_D6_input.sorted.bam

echo mES_m6A-seq_SetD2-KD_D6_IP.fastq
tophat -p 32 -g 1 --library-type=fr-firststrand -G /data/sunwj/genomeInfo/mouse/mm10/annotation/gencode.vM9.annotation.gtf -o ./SetD2-KD_D6/IP/ /public/genomes/mouse/mm10/bowtie2Index_UCSC/mm10 mES_m6A-seq_SetD2-KD_D6_IP.fastq
rename ./SetD2-KD_D6/IP/accepted_hits.bam ./SetD2-KD_D6/IP/mES_m6A-seq_SetD2-KD_D6_IP.bam ./SetD2-KD_D6/IP/accepted_hits.bam
mv ./SetD2-KD_D6/IP/mES_m6A-seq_SetD2-KD_D6_IP.bam ./SetD2-KD_D6/
samtools sort --threads 16 -m 2G -O bam -o ./SetD2-KD_D6/mES_m6A-seq_SetD2-KD_D6_IP.sorted.bam ./SetD2-KD_D6/mES_m6A-seq_SetD2-KD_D6_IP.bam
rm ./SetD2-KD_D6/mES_m6A-seq_SetD2-KD_D6_IP.bam
samtools index -b ./SetD2-KD_D6/mES_m6A-seq_SetD2-KD_D6_IP.sorted.bam

