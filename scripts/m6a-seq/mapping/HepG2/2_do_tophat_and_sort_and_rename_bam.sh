echo "\n\n\n"
echo HepG2_m6A-seq_shCont_IP_rep1.fastq
tophat -p 48 -g 1 --library-type=fr-firststrand -G /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.gtf -o ./shCont/IP_rep1/ /public/genomes/human/hg19/bowtie2Index_UCSC/hg19 /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/HepG2_m6A-seq_shCont_IP_rep1.fastq
rename ./shCont/IP_rep1/accepted_hits.bam ./shCont/IP_rep1/HepG2_m6A-seq_shCont_IP_rep1.fastq.bam ./shCont/IP_rep1/accepted_hits.bam
mv ./shCont/IP_rep1/HepG2_m6A-seq_shCont_IP_rep1.fastq.bam ./shCont/
samtools sort --threads 16 -m 2G -O bam -o ./shCont/HepG2_m6A-seq_shCont_IP_rep1.fastq.sorted.bam ./shCont/HepG2_m6A-seq_shCont_IP_rep1.fastq.bam
rm ./shCont/HepG2_m6A-seq_shCont_IP_rep1.fastq.bam
samtools index -b ./shCont/HepG2_m6A-seq_shCont_IP_rep1.fastq.sorted.bam

echo "\n\n\n"
echo HepG2_m6A-seq_shCont_IP_rep2.fastq
tophat -p 48 -g 1 --library-type=fr-firststrand -G /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.gtf -o ./shCont/IP_rep2/ /public/genomes/human/hg19/bowtie2Index_UCSC/hg19 /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/HepG2_m6A-seq_shCont_IP_rep2.fastq
rename ./shCont/IP_rep2/accepted_hits.bam ./shCont/IP_rep2/HepG2_m6A-seq_shCont_IP_rep2.fastq.bam ./shCont/IP_rep2/accepted_hits.bam
mv ./shCont/IP_rep2/HepG2_m6A-seq_shCont_IP_rep2.fastq.bam ./shCont/
samtools sort --threads 16 -m 2G -O bam -o ./shCont/HepG2_m6A-seq_shCont_IP_rep2.fastq.sorted.bam ./shCont/HepG2_m6A-seq_shCont_IP_rep2.fastq.bam
rm ./shCont/HepG2_m6A-seq_shCont_IP_rep2.fastq.bam
samtools index -b ./shCont/HepG2_m6A-seq_shCont_IP_rep2.fastq.sorted.bam

echo "\n\n\n"
echo HepG2_m6A-seq_shCont_IP_rep3.fastq
tophat -p 48 -g 1 --library-type=fr-firststrand -G /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.gtf -o ./shCont/IP_rep3/ /public/genomes/human/hg19/bowtie2Index_UCSC/hg19 /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/HepG2_m6A-seq_shCont_IP_rep3.fastq
rename ./shCont/IP_rep3/accepted_hits.bam ./shCont/IP_rep3/HepG2_m6A-seq_shCont_IP_rep3.fastq.bam ./shCont/IP_rep3/accepted_hits.bam
mv ./shCont/IP_rep3/HepG2_m6A-seq_shCont_IP_rep3.fastq.bam ./shCont/
samtools sort --threads 16 -m 2G -O bam -o ./shCont/HepG2_m6A-seq_shCont_IP_rep3.fastq.sorted.bam ./shCont/HepG2_m6A-seq_shCont_IP_rep3.fastq.bam
rm ./shCont/HepG2_m6A-seq_shCont_IP_rep3.fastq.bam
samtools index -b ./shCont/HepG2_m6A-seq_shCont_IP_rep3.fastq.sorted.bam

echo "\n\n\n"
echo HepG2_m6A-seq_shCont_input_rep1.fastq
tophat -p 48 -g 1 --library-type=fr-firststrand -G /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.gtf -o ./shCont/input_rep1/ /public/genomes/human/hg19/bowtie2Index_UCSC/hg19 /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/HepG2_m6A-seq_shCont_input_rep1.fastq
rename ./shCont/input_rep1/accepted_hits.bam ./shCont/input_rep1/HepG2_m6A-seq_shCont_input_rep1.fastq.bam ./shCont/input_rep1/accepted_hits.bam
mv ./shCont/input_rep1/HepG2_m6A-seq_shCont_input_rep1.fastq.bam ./shCont/
samtools sort --threads 16 -m 2G -O bam -o ./shCont/HepG2_m6A-seq_shCont_input_rep1.fastq.sorted.bam ./shCont/HepG2_m6A-seq_shCont_input_rep1.fastq.bam
rm ./shCont/HepG2_m6A-seq_shCont_input_rep1.fastq.bam
samtools index -b ./shCont/HepG2_m6A-seq_shCont_input_rep1.fastq.sorted.bam

echo "\n\n\n"
echo HepG2_m6A-seq_shCont_input_rep2.fastq
tophat -p 48 -g 1 --library-type=fr-firststrand -G /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.gtf -o ./shCont/input_rep2/ /public/genomes/human/hg19/bowtie2Index_UCSC/hg19 /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/HepG2_m6A-seq_shCont_input_rep2.fastq
rename ./shCont/input_rep2/accepted_hits.bam ./shCont/input_rep2/HepG2_m6A-seq_shCont_input_rep2.fastq.bam ./shCont/input_rep2/accepted_hits.bam
mv ./shCont/input_rep2/HepG2_m6A-seq_shCont_input_rep2.fastq.bam ./shCont/
samtools sort --threads 16 -m 2G -O bam -o ./shCont/HepG2_m6A-seq_shCont_input_rep2.fastq.sorted.bam ./shCont/HepG2_m6A-seq_shCont_input_rep2.fastq.bam
rm ./shCont/HepG2_m6A-seq_shCont_input_rep2.fastq.bam
samtools index -b ./shCont/HepG2_m6A-seq_shCont_input_rep2.fastq.sorted.bam

echo "\n\n\n"
echo HepG2_m6A-seq_shCont_input_rep3.fastq
tophat -p 48 -g 1 --library-type=fr-firststrand -G /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.gtf -o ./shCont/input_rep3/ /public/genomes/human/hg19/bowtie2Index_UCSC/hg19 /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/HepG2_m6A-seq_shCont_input_rep3.fastq
rename ./shCont/input_rep3/accepted_hits.bam ./shCont/input_rep3/HepG2_m6A-seq_shCont_input_rep3.fastq.bam ./shCont/input_rep3/accepted_hits.bam
mv ./shCont/input_rep3/HepG2_m6A-seq_shCont_input_rep3.fastq.bam ./shCont/
samtools sort --threads 16 -m 2G -O bam -o ./shCont/HepG2_m6A-seq_shCont_input_rep3.fastq.sorted.bam ./shCont/HepG2_m6A-seq_shCont_input_rep3.fastq.bam
rm ./shCont/HepG2_m6A-seq_shCont_input_rep3.fastq.bam
samtools index -b ./shCont/HepG2_m6A-seq_shCont_input_rep3.fastq.sorted.bam

echo "\n\n\n"
echo HepG2_m6A-seq_shM14_IP_rep1.fastq
tophat -p 48 -g 1 --library-type=fr-firststrand -G /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.gtf -o ./shM14/IP_rep1/ /public/genomes/human/hg19/bowtie2Index_UCSC/hg19 /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/HepG2_m6A-seq_shM14_IP_rep1.fastq
rename ./shM14/IP_rep1/accepted_hits.bam ./shM14/IP_rep1/HepG2_m6A-seq_shM14_IP_rep1.fastq.bam ./shM14/IP_rep1/accepted_hits.bam
mv ./shM14/IP_rep1/HepG2_m6A-seq_shM14_IP_rep1.fastq.bam ./shM14/
samtools sort --threads 16 -m 2G -O bam -o ./shM14/HepG2_m6A-seq_shM14_IP_rep1.fastq.sorted.bam ./shM14/HepG2_m6A-seq_shM14_IP_rep1.fastq.bam
rm ./shM14/HepG2_m6A-seq_shM14_IP_rep1.fastq.bam
samtools index -b ./shM14/HepG2_m6A-seq_shM14_IP_rep1.fastq.sorted.bam

echo "\n\n\n"
echo HepG2_m6A-seq_shM14_IP_rep2.fastq
tophat -p 48 -g 1 --library-type=fr-firststrand -G /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.gtf -o ./shM14/IP_rep2/ /public/genomes/human/hg19/bowtie2Index_UCSC/hg19 /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/HepG2_m6A-seq_shM14_IP_rep2.fastq
rename ./shM14/IP_rep2/accepted_hits.bam ./shM14/IP_rep2/HepG2_m6A-seq_shM14_IP_rep2.fastq.bam ./shM14/IP_rep2/accepted_hits.bam
mv ./shM14/IP_rep2/HepG2_m6A-seq_shM14_IP_rep2.fastq.bam ./shM14/
samtools sort --threads 16 -m 2G -O bam -o ./shM14/HepG2_m6A-seq_shM14_IP_rep2.fastq.sorted.bam ./shM14/HepG2_m6A-seq_shM14_IP_rep2.fastq.bam
rm ./shM14/HepG2_m6A-seq_shM14_IP_rep2.fastq.bam
samtools index -b ./shM14/HepG2_m6A-seq_shM14_IP_rep2.fastq.sorted.bam

echo "\n\n\n"
echo HepG2_m6A-seq_shM14_IP_rep3.fastq
tophat -p 48 -g 1 --library-type=fr-firststrand -G /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.gtf -o ./shM14/IP_rep3/ /public/genomes/human/hg19/bowtie2Index_UCSC/hg19 /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/HepG2_m6A-seq_shM14_IP_rep3.fastq
rename ./shM14/IP_rep3/accepted_hits.bam ./shM14/IP_rep3/HepG2_m6A-seq_shM14_IP_rep3.fastq.bam ./shM14/IP_rep3/accepted_hits.bam
mv ./shM14/IP_rep3/HepG2_m6A-seq_shM14_IP_rep3.fastq.bam ./shM14/
samtools sort --threads 16 -m 2G -O bam -o ./shM14/HepG2_m6A-seq_shM14_IP_rep3.fastq.sorted.bam ./shM14/HepG2_m6A-seq_shM14_IP_rep3.fastq.bam
rm ./shM14/HepG2_m6A-seq_shM14_IP_rep3.fastq.bam
samtools index -b ./shM14/HepG2_m6A-seq_shM14_IP_rep3.fastq.sorted.bam

echo "\n\n\n"
echo HepG2_m6A-seq_shM14_input_rep1.fastq
tophat -p 48 -g 1 --library-type=fr-firststrand -G /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.gtf -o ./shM14/input_rep1/ /public/genomes/human/hg19/bowtie2Index_UCSC/hg19 /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/HepG2_m6A-seq_shM14_input_rep1.fastq
rename ./shM14/input_rep1/accepted_hits.bam ./shM14/input_rep1/HepG2_m6A-seq_shM14_input_rep1.fastq.bam ./shM14/input_rep1/accepted_hits.bam
mv ./shM14/input_rep1/HepG2_m6A-seq_shM14_input_rep1.fastq.bam ./shM14/
samtools sort --threads 16 -m 2G -O bam -o ./shM14/HepG2_m6A-seq_shM14_input_rep1.fastq.sorted.bam ./shM14/HepG2_m6A-seq_shM14_input_rep1.fastq.bam
rm ./shM14/HepG2_m6A-seq_shM14_input_rep1.fastq.bam
samtools index -b ./shM14/HepG2_m6A-seq_shM14_input_rep1.fastq.sorted.bam

echo "\n\n\n"
echo HepG2_m6A-seq_shM14_input_rep2.fastq
tophat -p 48 -g 1 --library-type=fr-firststrand -G /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.gtf -o ./shM14/input_rep2/ /public/genomes/human/hg19/bowtie2Index_UCSC/hg19 /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/HepG2_m6A-seq_shM14_input_rep2.fastq
rename ./shM14/input_rep2/accepted_hits.bam ./shM14/input_rep2/HepG2_m6A-seq_shM14_input_rep2.fastq.bam ./shM14/input_rep2/accepted_hits.bam
mv ./shM14/input_rep2/HepG2_m6A-seq_shM14_input_rep2.fastq.bam ./shM14/
samtools sort --threads 16 -m 2G -O bam -o ./shM14/HepG2_m6A-seq_shM14_input_rep2.fastq.sorted.bam ./shM14/HepG2_m6A-seq_shM14_input_rep2.fastq.bam
rm ./shM14/HepG2_m6A-seq_shM14_input_rep2.fastq.bam
samtools index -b ./shM14/HepG2_m6A-seq_shM14_input_rep2.fastq.sorted.bam

echo "\n\n\n"
echo HepG2_m6A-seq_shM14_input_rep3.fastq
tophat -p 48 -g 1 --library-type=fr-firststrand -G /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.gtf -o ./shM14/input_rep3/ /public/genomes/human/hg19/bowtie2Index_UCSC/hg19 /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/HepG2_m6A-seq_shM14_input_rep3.fastq
rename ./shM14/input_rep3/accepted_hits.bam ./shM14/input_rep3/HepG2_m6A-seq_shM14_input_rep3.fastq.bam ./shM14/input_rep3/accepted_hits.bam
mv ./shM14/input_rep3/HepG2_m6A-seq_shM14_input_rep3.fastq.bam ./shM14/
samtools sort --threads 16 -m 2G -O bam -o ./shM14/HepG2_m6A-seq_shM14_input_rep3.fastq.sorted.bam ./shM14/HepG2_m6A-seq_shM14_input_rep3.fastq.bam
rm ./shM14/HepG2_m6A-seq_shM14_input_rep3.fastq.bam
samtools index -b ./shM14/HepG2_m6A-seq_shM14_input_rep3.fastq.sorted.bam

echo "\n\n\n"
echo HepG2_m6A-seq_shM3_IP_rep1.fastq
tophat -p 48 -g 1 --library-type=fr-firststrand -G /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.gtf -o ./shM3/IP_rep1/ /public/genomes/human/hg19/bowtie2Index_UCSC/hg19 /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/HepG2_m6A-seq_shM3_IP_rep1.fastq
rename ./shM3/IP_rep1/accepted_hits.bam ./shM3/IP_rep1/HepG2_m6A-seq_shM3_IP_rep1.fastq.bam ./shM3/IP_rep1/accepted_hits.bam
mv ./shM3/IP_rep1/HepG2_m6A-seq_shM3_IP_rep1.fastq.bam ./shM3/
samtools sort --threads 16 -m 2G -O bam -o ./shM3/HepG2_m6A-seq_shM3_IP_rep1.fastq.sorted.bam ./shM3/HepG2_m6A-seq_shM3_IP_rep1.fastq.bam
rm ./shM3/HepG2_m6A-seq_shM3_IP_rep1.fastq.bam
samtools index -b ./shM3/HepG2_m6A-seq_shM3_IP_rep1.fastq.sorted.bam

echo "\n\n\n"
echo HepG2_m6A-seq_shM3_IP_rep2.fastq
tophat -p 48 -g 1 --library-type=fr-firststrand -G /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.gtf -o ./shM3/IP_rep2/ /public/genomes/human/hg19/bowtie2Index_UCSC/hg19 /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/HepG2_m6A-seq_shM3_IP_rep2.fastq
rename ./shM3/IP_rep2/accepted_hits.bam ./shM3/IP_rep2/HepG2_m6A-seq_shM3_IP_rep2.fastq.bam ./shM3/IP_rep2/accepted_hits.bam
mv ./shM3/IP_rep2/HepG2_m6A-seq_shM3_IP_rep2.fastq.bam ./shM3/
samtools sort --threads 16 -m 2G -O bam -o ./shM3/HepG2_m6A-seq_shM3_IP_rep2.fastq.sorted.bam ./shM3/HepG2_m6A-seq_shM3_IP_rep2.fastq.bam
rm ./shM3/HepG2_m6A-seq_shM3_IP_rep2.fastq.bam
samtools index -b ./shM3/HepG2_m6A-seq_shM3_IP_rep2.fastq.sorted.bam

echo "\n\n\n"
echo HepG2_m6A-seq_shM3_IP_rep3.fastq
tophat -p 48 -g 1 --library-type=fr-firststrand -G /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.gtf -o ./shM3/IP_rep3/ /public/genomes/human/hg19/bowtie2Index_UCSC/hg19 /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/HepG2_m6A-seq_shM3_IP_rep3.fastq
rename ./shM3/IP_rep3/accepted_hits.bam ./shM3/IP_rep3/HepG2_m6A-seq_shM3_IP_rep3.fastq.bam ./shM3/IP_rep3/accepted_hits.bam
mv ./shM3/IP_rep3/HepG2_m6A-seq_shM3_IP_rep3.fastq.bam ./shM3/
samtools sort --threads 16 -m 2G -O bam -o ./shM3/HepG2_m6A-seq_shM3_IP_rep3.fastq.sorted.bam ./shM3/HepG2_m6A-seq_shM3_IP_rep3.fastq.bam
rm ./shM3/HepG2_m6A-seq_shM3_IP_rep3.fastq.bam
samtools index -b ./shM3/HepG2_m6A-seq_shM3_IP_rep3.fastq.sorted.bam

echo "\n\n\n"
echo HepG2_m6A-seq_shM3_input_rep1.fastq
tophat -p 48 -g 1 --library-type=fr-firststrand -G /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.gtf -o ./shM3/input_rep1/ /public/genomes/human/hg19/bowtie2Index_UCSC/hg19 /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/HepG2_m6A-seq_shM3_input_rep1.fastq
rename ./shM3/input_rep1/accepted_hits.bam ./shM3/input_rep1/HepG2_m6A-seq_shM3_input_rep1.fastq.bam ./shM3/input_rep1/accepted_hits.bam
mv ./shM3/input_rep1/HepG2_m6A-seq_shM3_input_rep1.fastq.bam ./shM3/
samtools sort --threads 16 -m 2G -O bam -o ./shM3/HepG2_m6A-seq_shM3_input_rep1.fastq.sorted.bam ./shM3/HepG2_m6A-seq_shM3_input_rep1.fastq.bam
rm ./shM3/HepG2_m6A-seq_shM3_input_rep1.fastq.bam
samtools index -b ./shM3/HepG2_m6A-seq_shM3_input_rep1.fastq.sorted.bam

echo "\n\n\n"
echo HepG2_m6A-seq_shM3_input_rep2.fastq
tophat -p 48 -g 1 --library-type=fr-firststrand -G /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.gtf -o ./shM3/input_rep2/ /public/genomes/human/hg19/bowtie2Index_UCSC/hg19 /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/HepG2_m6A-seq_shM3_input_rep2.fastq
rename ./shM3/input_rep2/accepted_hits.bam ./shM3/input_rep2/HepG2_m6A-seq_shM3_input_rep2.fastq.bam ./shM3/input_rep2/accepted_hits.bam
mv ./shM3/input_rep2/HepG2_m6A-seq_shM3_input_rep2.fastq.bam ./shM3/
samtools sort --threads 16 -m 2G -O bam -o ./shM3/HepG2_m6A-seq_shM3_input_rep2.fastq.sorted.bam ./shM3/HepG2_m6A-seq_shM3_input_rep2.fastq.bam
rm ./shM3/HepG2_m6A-seq_shM3_input_rep2.fastq.bam
samtools index -b ./shM3/HepG2_m6A-seq_shM3_input_rep2.fastq.sorted.bam

echo "\n\n\n"
echo HepG2_m6A-seq_shM3_input_rep3.fastq
tophat -p 48 -g 1 --library-type=fr-firststrand -G /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.gtf -o ./shM3/input_rep3/ /public/genomes/human/hg19/bowtie2Index_UCSC/hg19 /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/HepG2_m6A-seq_shM3_input_rep3.fastq
rename ./shM3/input_rep3/accepted_hits.bam ./shM3/input_rep3/HepG2_m6A-seq_shM3_input_rep3.fastq.bam ./shM3/input_rep3/accepted_hits.bam
mv ./shM3/input_rep3/HepG2_m6A-seq_shM3_input_rep3.fastq.bam ./shM3/
samtools sort --threads 16 -m 2G -O bam -o ./shM3/HepG2_m6A-seq_shM3_input_rep3.fastq.sorted.bam ./shM3/HepG2_m6A-seq_shM3_input_rep3.fastq.bam
rm ./shM3/HepG2_m6A-seq_shM3_input_rep3.fastq.bam
samtools index -b ./shM3/HepG2_m6A-seq_shM3_input_rep3.fastq.sorted.bam

echo "\n\n\n"
echo HepG2_m6A-seq_shSetD2_IP_rep1.fastq
tophat -p 48 -g 1 --library-type=fr-firststrand -G /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.gtf -o ./shSetD2/IP_rep1/ /public/genomes/human/hg19/bowtie2Index_UCSC/hg19 /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/HepG2_m6A-seq_shSetD2_IP_rep1.fastq
rename ./shSetD2/IP_rep1/accepted_hits.bam ./shSetD2/IP_rep1/HepG2_m6A-seq_shSetD2_IP_rep1.fastq.bam ./shSetD2/IP_rep1/accepted_hits.bam
mv ./shSetD2/IP_rep1/HepG2_m6A-seq_shSetD2_IP_rep1.fastq.bam ./shSetD2/
samtools sort --threads 16 -m 2G -O bam -o ./shSetD2/HepG2_m6A-seq_shSetD2_IP_rep1.fastq.sorted.bam ./shSetD2/HepG2_m6A-seq_shSetD2_IP_rep1.fastq.bam
rm ./shSetD2/HepG2_m6A-seq_shSetD2_IP_rep1.fastq.bam
samtools index -b ./shSetD2/HepG2_m6A-seq_shSetD2_IP_rep1.fastq.sorted.bam

echo "\n\n\n"
echo HepG2_m6A-seq_shSetD2_IP_rep2.fastq
tophat -p 48 -g 1 --library-type=fr-firststrand -G /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.gtf -o ./shSetD2/IP_rep2/ /public/genomes/human/hg19/bowtie2Index_UCSC/hg19 /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/HepG2_m6A-seq_shSetD2_IP_rep2.fastq
rename ./shSetD2/IP_rep2/accepted_hits.bam ./shSetD2/IP_rep2/HepG2_m6A-seq_shSetD2_IP_rep2.fastq.bam ./shSetD2/IP_rep2/accepted_hits.bam
mv ./shSetD2/IP_rep2/HepG2_m6A-seq_shSetD2_IP_rep2.fastq.bam ./shSetD2/
samtools sort --threads 16 -m 2G -O bam -o ./shSetD2/HepG2_m6A-seq_shSetD2_IP_rep2.fastq.sorted.bam ./shSetD2/HepG2_m6A-seq_shSetD2_IP_rep2.fastq.bam
rm ./shSetD2/HepG2_m6A-seq_shSetD2_IP_rep2.fastq.bam
samtools index -b ./shSetD2/HepG2_m6A-seq_shSetD2_IP_rep2.fastq.sorted.bam

echo "\n\n\n"
echo HepG2_m6A-seq_shSetD2_IP_rep3.fastq
tophat -p 48 -g 1 --library-type=fr-firststrand -G /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.gtf -o ./shSetD2/IP_rep3/ /public/genomes/human/hg19/bowtie2Index_UCSC/hg19 /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/HepG2_m6A-seq_shSetD2_IP_rep3.fastq
rename ./shSetD2/IP_rep3/accepted_hits.bam ./shSetD2/IP_rep3/HepG2_m6A-seq_shSetD2_IP_rep3.fastq.bam ./shSetD2/IP_rep3/accepted_hits.bam
mv ./shSetD2/IP_rep3/HepG2_m6A-seq_shSetD2_IP_rep3.fastq.bam ./shSetD2/
samtools sort --threads 16 -m 2G -O bam -o ./shSetD2/HepG2_m6A-seq_shSetD2_IP_rep3.fastq.sorted.bam ./shSetD2/HepG2_m6A-seq_shSetD2_IP_rep3.fastq.bam
rm ./shSetD2/HepG2_m6A-seq_shSetD2_IP_rep3.fastq.bam
samtools index -b ./shSetD2/HepG2_m6A-seq_shSetD2_IP_rep3.fastq.sorted.bam

echo "\n\n\n"
echo HepG2_m6A-seq_shSetD2_input_rep1.fastq
tophat -p 48 -g 1 --library-type=fr-firststrand -G /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.gtf -o ./shSetD2/input_rep1/ /public/genomes/human/hg19/bowtie2Index_UCSC/hg19 /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/HepG2_m6A-seq_shSetD2_input_rep1.fastq
rename ./shSetD2/input_rep1/accepted_hits.bam ./shSetD2/input_rep1/HepG2_m6A-seq_shSetD2_input_rep1.fastq.bam ./shSetD2/input_rep1/accepted_hits.bam
mv ./shSetD2/input_rep1/HepG2_m6A-seq_shSetD2_input_rep1.fastq.bam ./shSetD2/
samtools sort --threads 16 -m 2G -O bam -o ./shSetD2/HepG2_m6A-seq_shSetD2_input_rep1.fastq.sorted.bam ./shSetD2/HepG2_m6A-seq_shSetD2_input_rep1.fastq.bam
rm ./shSetD2/HepG2_m6A-seq_shSetD2_input_rep1.fastq.bam
samtools index -b ./shSetD2/HepG2_m6A-seq_shSetD2_input_rep1.fastq.sorted.bam

echo "\n\n\n"
echo HepG2_m6A-seq_shSetD2_input_rep2.fastq
tophat -p 48 -g 1 --library-type=fr-firststrand -G /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.gtf -o ./shSetD2/input_rep2/ /public/genomes/human/hg19/bowtie2Index_UCSC/hg19 /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/HepG2_m6A-seq_shSetD2_input_rep2.fastq
rename ./shSetD2/input_rep2/accepted_hits.bam ./shSetD2/input_rep2/HepG2_m6A-seq_shSetD2_input_rep2.fastq.bam ./shSetD2/input_rep2/accepted_hits.bam
mv ./shSetD2/input_rep2/HepG2_m6A-seq_shSetD2_input_rep2.fastq.bam ./shSetD2/
samtools sort --threads 16 -m 2G -O bam -o ./shSetD2/HepG2_m6A-seq_shSetD2_input_rep2.fastq.sorted.bam ./shSetD2/HepG2_m6A-seq_shSetD2_input_rep2.fastq.bam
rm ./shSetD2/HepG2_m6A-seq_shSetD2_input_rep2.fastq.bam
samtools index -b ./shSetD2/HepG2_m6A-seq_shSetD2_input_rep2.fastq.sorted.bam

echo "\n\n\n"
echo HepG2_m6A-seq_shSetD2_input_rep3.fastq
tophat -p 48 -g 1 --library-type=fr-firststrand -G /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.gtf -o ./shSetD2/input_rep3/ /public/genomes/human/hg19/bowtie2Index_UCSC/hg19 /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/HepG2_m6A-seq_shSetD2_input_rep3.fastq
rename ./shSetD2/input_rep3/accepted_hits.bam ./shSetD2/input_rep3/HepG2_m6A-seq_shSetD2_input_rep3.fastq.bam ./shSetD2/input_rep3/accepted_hits.bam
mv ./shSetD2/input_rep3/HepG2_m6A-seq_shSetD2_input_rep3.fastq.bam ./shSetD2/
samtools sort --threads 16 -m 2G -O bam -o ./shSetD2/HepG2_m6A-seq_shSetD2_input_rep3.fastq.sorted.bam ./shSetD2/HepG2_m6A-seq_shSetD2_input_rep3.fastq.bam
rm ./shSetD2/HepG2_m6A-seq_shSetD2_input_rep3.fastq.bam
samtools index -b ./shSetD2/HepG2_m6A-seq_shSetD2_input_rep3.fastq.sorted.bam

echo "\n\n\n"
echo HepG2_m6A-seq_shWTAP_IP_rep1.fastq
tophat -p 48 -g 1 --library-type=fr-firststrand -G /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.gtf -o ./shWTAP/IP_rep1/ /public/genomes/human/hg19/bowtie2Index_UCSC/hg19 /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/HepG2_m6A-seq_shWTAP_IP_rep1.fastq
rename ./shWTAP/IP_rep1/accepted_hits.bam ./shWTAP/IP_rep1/HepG2_m6A-seq_shWTAP_IP_rep1.fastq.bam ./shWTAP/IP_rep1/accepted_hits.bam
mv ./shWTAP/IP_rep1/HepG2_m6A-seq_shWTAP_IP_rep1.fastq.bam ./shWTAP/
samtools sort --threads 16 -m 2G -O bam -o ./shWTAP/HepG2_m6A-seq_shWTAP_IP_rep1.fastq.sorted.bam ./shWTAP/HepG2_m6A-seq_shWTAP_IP_rep1.fastq.bam
rm ./shWTAP/HepG2_m6A-seq_shWTAP_IP_rep1.fastq.bam
samtools index -b ./shWTAP/HepG2_m6A-seq_shWTAP_IP_rep1.fastq.sorted.bam

echo "\n\n\n"
echo HepG2_m6A-seq_shWTAP_IP_rep2.fastq
tophat -p 48 -g 1 --library-type=fr-firststrand -G /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.gtf -o ./shWTAP/IP_rep2/ /public/genomes/human/hg19/bowtie2Index_UCSC/hg19 /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/HepG2_m6A-seq_shWTAP_IP_rep2.fastq
rename ./shWTAP/IP_rep2/accepted_hits.bam ./shWTAP/IP_rep2/HepG2_m6A-seq_shWTAP_IP_rep2.fastq.bam ./shWTAP/IP_rep2/accepted_hits.bam
mv ./shWTAP/IP_rep2/HepG2_m6A-seq_shWTAP_IP_rep2.fastq.bam ./shWTAP/
samtools sort --threads 16 -m 2G -O bam -o ./shWTAP/HepG2_m6A-seq_shWTAP_IP_rep2.fastq.sorted.bam ./shWTAP/HepG2_m6A-seq_shWTAP_IP_rep2.fastq.bam
rm ./shWTAP/HepG2_m6A-seq_shWTAP_IP_rep2.fastq.bam
samtools index -b ./shWTAP/HepG2_m6A-seq_shWTAP_IP_rep2.fastq.sorted.bam

echo "\n\n\n"
echo HepG2_m6A-seq_shWTAP_IP_rep3.fastq
tophat -p 48 -g 1 --library-type=fr-firststrand -G /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.gtf -o ./shWTAP/IP_rep3/ /public/genomes/human/hg19/bowtie2Index_UCSC/hg19 /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/HepG2_m6A-seq_shWTAP_IP_rep3.fastq
rename ./shWTAP/IP_rep3/accepted_hits.bam ./shWTAP/IP_rep3/HepG2_m6A-seq_shWTAP_IP_rep3.fastq.bam ./shWTAP/IP_rep3/accepted_hits.bam
mv ./shWTAP/IP_rep3/HepG2_m6A-seq_shWTAP_IP_rep3.fastq.bam ./shWTAP/
samtools sort --threads 16 -m 2G -O bam -o ./shWTAP/HepG2_m6A-seq_shWTAP_IP_rep3.fastq.sorted.bam ./shWTAP/HepG2_m6A-seq_shWTAP_IP_rep3.fastq.bam
rm ./shWTAP/HepG2_m6A-seq_shWTAP_IP_rep3.fastq.bam
samtools index -b ./shWTAP/HepG2_m6A-seq_shWTAP_IP_rep3.fastq.sorted.bam

echo "\n\n\n"
echo HepG2_m6A-seq_shWTAP_input_rep1.fastq
tophat -p 48 -g 1 --library-type=fr-firststrand -G /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.gtf -o ./shWTAP/input_rep1/ /public/genomes/human/hg19/bowtie2Index_UCSC/hg19 /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/HepG2_m6A-seq_shWTAP_input_rep1.fastq
rename ./shWTAP/input_rep1/accepted_hits.bam ./shWTAP/input_rep1/HepG2_m6A-seq_shWTAP_input_rep1.fastq.bam ./shWTAP/input_rep1/accepted_hits.bam
mv ./shWTAP/input_rep1/HepG2_m6A-seq_shWTAP_input_rep1.fastq.bam ./shWTAP/
samtools sort --threads 16 -m 2G -O bam -o ./shWTAP/HepG2_m6A-seq_shWTAP_input_rep1.fastq.sorted.bam ./shWTAP/HepG2_m6A-seq_shWTAP_input_rep1.fastq.bam
rm ./shWTAP/HepG2_m6A-seq_shWTAP_input_rep1.fastq.bam
samtools index -b ./shWTAP/HepG2_m6A-seq_shWTAP_input_rep1.fastq.sorted.bam

echo "\n\n\n"
echo HepG2_m6A-seq_shWTAP_input_rep2.fastq
tophat -p 48 -g 1 --library-type=fr-firststrand -G /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.gtf -o ./shWTAP/input_rep2/ /public/genomes/human/hg19/bowtie2Index_UCSC/hg19 /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/HepG2_m6A-seq_shWTAP_input_rep2.fastq
rename ./shWTAP/input_rep2/accepted_hits.bam ./shWTAP/input_rep2/HepG2_m6A-seq_shWTAP_input_rep2.fastq.bam ./shWTAP/input_rep2/accepted_hits.bam
mv ./shWTAP/input_rep2/HepG2_m6A-seq_shWTAP_input_rep2.fastq.bam ./shWTAP/
samtools sort --threads 16 -m 2G -O bam -o ./shWTAP/HepG2_m6A-seq_shWTAP_input_rep2.fastq.sorted.bam ./shWTAP/HepG2_m6A-seq_shWTAP_input_rep2.fastq.bam
rm ./shWTAP/HepG2_m6A-seq_shWTAP_input_rep2.fastq.bam
samtools index -b ./shWTAP/HepG2_m6A-seq_shWTAP_input_rep2.fastq.sorted.bam

echo "\n\n\n"
echo HepG2_m6A-seq_shWTAP_input_rep3.fastq
tophat -p 48 -g 1 --library-type=fr-firststrand -G /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.gtf -o ./shWTAP/input_rep3/ /public/genomes/human/hg19/bowtie2Index_UCSC/hg19 /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/HepG2_m6A-seq_shWTAP_input_rep3.fastq
rename ./shWTAP/input_rep3/accepted_hits.bam ./shWTAP/input_rep3/HepG2_m6A-seq_shWTAP_input_rep3.fastq.bam ./shWTAP/input_rep3/accepted_hits.bam
mv ./shWTAP/input_rep3/HepG2_m6A-seq_shWTAP_input_rep3.fastq.bam ./shWTAP/
samtools sort --threads 16 -m 2G -O bam -o ./shWTAP/HepG2_m6A-seq_shWTAP_input_rep3.fastq.sorted.bam ./shWTAP/HepG2_m6A-seq_shWTAP_input_rep3.fastq.bam
rm ./shWTAP/HepG2_m6A-seq_shWTAP_input_rep3.fastq.bam
samtools index -b ./shWTAP/HepG2_m6A-seq_shWTAP_input_rep3.fastq.sorted.bam

