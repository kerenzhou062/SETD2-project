echo "HepG2_m6A-seq_shCont_input.fastq"
hisat2 -p 3 --rna-strandness R --dta -x /data/zhoukr/reference/genome/Index/Hisat2_index/hg19 -q /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/HepG2_m6A-seq_shCont_input.fastq -S /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shCont/HepG2_m6A-seq_shCont_input.sam
echo -e "transcripts alignment done. sam to index bam by using samtools...\n"
samtools sort --threads 8 -O bam -o /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shCont/HepG2_m6A-seq_shCont_input.sorted.bam /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shCont/HepG2_m6A-seq_shCont_input.sam
samtools index -b /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shCont/HepG2_m6A-seq_shCont_input.sorted.bam
echo -e "stringtie: start assbemling transcriptome...\n"
stringtie -p 8 --rf -G /data/zhoukr/reference/genome/gtf/gencode.v24lift37.annotation.gtf -o /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shCont/HepG2_m6A-seq_shCont_input.gtf /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shCont/HepG2_m6A-seq_shCont_input.sorted.bam
sed -i 's/STRG/HepG2_m6A-seq_shCont_input/g' /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shCont/HepG2_m6A-seq_shCont_input.gtf
rm /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shCont/HepG2_m6A-seq_shCont_input.sam

gffcompare -r /data/zhoukr/reference/genome/gtf/gencode.v24lift37.annotation.gtf -o /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shCont/shCont_input_compare /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shCont/HepG2_m6A-seq_shCont_input.gtf

echo "HepG2_m6A-seq_shSetD2_input.fastq"
hisat2 -p 3 --rna-strandness R --dta -x /data/zhoukr/reference/genome/Index/Hisat2_index/hg19 -q /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/HepG2_m6A-seq_shSetD2_input.fastq -S /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shSetD2/HepG2_m6A-seq_shSetD2_input.sam
echo -e "transcripts alignment done. sam to index bam by using samtools...\n"
samtools sort --threads 8 -O bam -o /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shSetD2/HepG2_m6A-seq_shSetD2_input.sorted.bam /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shSetD2/HepG2_m6A-seq_shSetD2_input.sam
samtools index -b /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shSetD2/HepG2_m6A-seq_shSetD2_input.sorted.bam
echo -e "stringtie: start assbemling transcriptome...\n"
stringtie -p 8 --rf -G /data/zhoukr/reference/genome/gtf/gencode.v24lift37.annotation.gtf -o /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shSetD2/HepG2_m6A-seq_shSetD2_input.gtf /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shSetD2/HepG2_m6A-seq_shSetD2_input.sorted.bam
sed -i 's/STRG/HepG2_m6A-seq_shSetD2_input/g' /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shSetD2/HepG2_m6A-seq_shSetD2_input.gtf
rm /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shSetD2/HepG2_m6A-seq_shSetD2_input.sam

gffcompare -r /data/zhoukr/reference/genome/gtf/gencode.v24lift37.annotation.gtf -o /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shSetD2/shSetD2_input_compare /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shSetD2/HepG2_m6A-seq_shSetD2_input.gtf

echo "HepG2_m6A-seq_shM3_input.fastq"
hisat2 -p 3 --rna-strandness R --dta -x /data/zhoukr/reference/genome/Index/Hisat2_index/hg19 -q /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/HepG2_m6A-seq_shM3_input.fastq -S /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shM3/HepG2_m6A-seq_shM3_input.sam
echo -e "transcripts alignment done. sam to index bam by using samtools...\n"
samtools sort --threads 8 -O bam -o /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shM3/HepG2_m6A-seq_shM3_input.sorted.bam /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shM3/HepG2_m6A-seq_shM3_input.sam
samtools index -b /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shM3/HepG2_m6A-seq_shM3_input.sorted.bam
echo -e "stringtie: start assbemling transcriptome...\n"
stringtie -p 8 --rf -G /data/zhoukr/reference/genome/gtf/gencode.v24lift37.annotation.gtf -o /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shM3/HepG2_m6A-seq_shM3_input.gtf /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shM3/HepG2_m6A-seq_shM3_input.sorted.bam
sed -i 's/STRG/HepG2_m6A-seq_shM3_input/g' /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shM3/HepG2_m6A-seq_shM3_input.gtf
rm /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shM3/HepG2_m6A-seq_shM3_input.sam

gffcompare -r /data/zhoukr/reference/genome/gtf/gencode.v24lift37.annotation.gtf -o /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shM3/shM3_input_compare /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shM3/HepG2_m6A-seq_shM3_input.gtf

echo "HepG2_m6A-seq_shM14_input.fastq"
hisat2 -p 3 --rna-strandness R --dta -x /data/zhoukr/reference/genome/Index/Hisat2_index/hg19 -q /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/HepG2_m6A-seq_shM14_input.fastq -S /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shM14/HepG2_m6A-seq_shM14_input.sam
echo -e "transcripts alignment done. sam to index bam by using samtools...\n"
samtools sort --threads 8 -O bam -o /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shM14/HepG2_m6A-seq_shM14_input.sorted.bam /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shM14/HepG2_m6A-seq_shM14_input.sam
samtools index -b /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shM14/HepG2_m6A-seq_shM14_input.sorted.bam
echo -e "stringtie: start assbemling transcriptome...\n"
stringtie -p 8 --rf -G /data/zhoukr/reference/genome/gtf/gencode.v24lift37.annotation.gtf -o /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shM14/HepG2_m6A-seq_shM14_input.gtf /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shM14/HepG2_m6A-seq_shM14_input.sorted.bam
sed -i 's/STRG/HepG2_m6A-seq_shM14_input/g' /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shM14/HepG2_m6A-seq_shM14_input.gtf
rm /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shM14/HepG2_m6A-seq_shM14_input.sam

gffcompare -r /data/zhoukr/reference/genome/gtf/gencode.v24lift37.annotation.gtf -o /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shM14/shM14_input_compare /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shM14/HepG2_m6A-seq_shM14_input.gtf

echo "HepG2_m6A-seq_shWTAP_input.fastq"
hisat2 -p 3 --rna-strandness R --dta -x /data/zhoukr/reference/genome/Index/Hisat2_index/hg19 -q /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/HepG2_m6A-seq_shWTAP_input.fastq -S /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shWTAP/HepG2_m6A-seq_shWTAP_input.sam
echo -e "transcripts alignment done. sam to index bam by using samtools...\n"
samtools sort --threads 8 -O bam -o /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shWTAP/HepG2_m6A-seq_shWTAP_input.sorted.bam /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shWTAP/HepG2_m6A-seq_shWTAP_input.sam
samtools index -b /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shWTAP/HepG2_m6A-seq_shWTAP_input.sorted.bam
echo -e "stringtie: start assbemling transcriptome...\n"
stringtie -p 8 --rf -G /data/zhoukr/reference/genome/gtf/gencode.v24lift37.annotation.gtf -o /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shWTAP/HepG2_m6A-seq_shWTAP_input.gtf /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shWTAP/HepG2_m6A-seq_shWTAP_input.sorted.bam
sed -i 's/STRG/HepG2_m6A-seq_shWTAP_input/g' /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shWTAP/HepG2_m6A-seq_shWTAP_input.gtf
rm /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shWTAP/HepG2_m6A-seq_shWTAP_input.sam

gffcompare -r /data/zhoukr/reference/genome/gtf/gencode.v24lift37.annotation.gtf -o /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shWTAP/shWTAP_input_compare /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shWTAP/HepG2_m6A-seq_shWTAP_input.gtf

echo -e "jobs have been done.\n"
