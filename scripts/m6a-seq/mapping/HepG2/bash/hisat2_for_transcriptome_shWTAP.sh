echo "HepG2_m6A-seq_shWTAP_input_rep1.fastq"
hisat2 -p 3 --rna-strandness R --dta -x /data/zhoukr/reference/genome/Index/Hisat2_index/hg19 -q /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/HepG2_m6A-seq_shWTAP_input_rep1.fastq -S /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shWTAP/HepG2_m6A-seq_shWTAP_input_rep1.sam
echo -e "transcripts alignment done. sam to index bam by using samtools...\n"
samtools sort --threads 8 -O bam -o /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shWTAP/HepG2_m6A-seq_shWTAP_input_rep1.sorted.bam /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shWTAP/HepG2_m6A-seq_shWTAP_input_rep1.sam
samtools index -b /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shWTAP/HepG2_m6A-seq_shWTAP_input_rep1.sorted.bam
echo -e "stringtie: start assbemling transcriptome...\n"
stringtie -p 8 --rf -G /data/zhoukr/reference/genome/gtf/gencode.v24lift37.annotation.gtf -o /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shWTAP/HepG2_m6A-seq_shWTAP_input_rep1.gtf /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shWTAP/HepG2_m6A-seq_shWTAP_input_rep1.sorted.bam
sed -i 's/STRG/HepG2_m6A-seq_shWTAP_input_rep1/g' /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shWTAP/HepG2_m6A-seq_shWTAP_input_rep1.gtf
rm /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shWTAP/HepG2_m6A-seq_shWTAP_input_rep1.sam
echo "HepG2_m6A-seq_shWTAP_input_rep2.fastq"
hisat2 -p 3 --rna-strandness R --dta -x /data/zhoukr/reference/genome/Index/Hisat2_index/hg19 -q /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/HepG2_m6A-seq_shWTAP_input_rep2.fastq -S /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shWTAP/HepG2_m6A-seq_shWTAP_input_rep2.sam
echo -e "transcripts alignment done. sam to index bam by using samtools...\n"
samtools sort --threads 8 -O bam -o /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shWTAP/HepG2_m6A-seq_shWTAP_input_rep2.sorted.bam /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shWTAP/HepG2_m6A-seq_shWTAP_input_rep2.sam
samtools index -b /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shWTAP/HepG2_m6A-seq_shWTAP_input_rep2.sorted.bam
echo -e "stringtie: start assbemling transcriptome...\n"
stringtie -p 8 --rf -G /data/zhoukr/reference/genome/gtf/gencode.v24lift37.annotation.gtf -o /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shWTAP/HepG2_m6A-seq_shWTAP_input_rep2.gtf /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shWTAP/HepG2_m6A-seq_shWTAP_input_rep2.sorted.bam
sed -i 's/STRG/HepG2_m6A-seq_shWTAP_input_rep2/g' /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shWTAP/HepG2_m6A-seq_shWTAP_input_rep2.gtf
rm /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shWTAP/HepG2_m6A-seq_shWTAP_input_rep2.sam
echo "HepG2_m6A-seq_shWTAP_input_rep3.fastq"
hisat2 -p 3 --rna-strandness R --dta -x /data/zhoukr/reference/genome/Index/Hisat2_index/hg19 -q /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/HepG2_m6A-seq_shWTAP_input_rep3.fastq -S /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shWTAP/HepG2_m6A-seq_shWTAP_input_rep3.sam
echo -e "transcripts alignment done. sam to index bam by using samtools...\n"
samtools sort --threads 8 -O bam -o /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shWTAP/HepG2_m6A-seq_shWTAP_input_rep3.sorted.bam /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shWTAP/HepG2_m6A-seq_shWTAP_input_rep3.sam
samtools index -b /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shWTAP/HepG2_m6A-seq_shWTAP_input_rep3.sorted.bam
echo -e "stringtie: start assbemling transcriptome...\n"
stringtie -p 8 --rf -G /data/zhoukr/reference/genome/gtf/gencode.v24lift37.annotation.gtf -o /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shWTAP/HepG2_m6A-seq_shWTAP_input_rep3.gtf /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shWTAP/HepG2_m6A-seq_shWTAP_input_rep3.sorted.bam
sed -i 's/STRG/HepG2_m6A-seq_shWTAP_input_rep3/g' /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shWTAP/HepG2_m6A-seq_shWTAP_input_rep3.gtf
rm /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shWTAP/HepG2_m6A-seq_shWTAP_input_rep3.sam
echo "/data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shWTAP/HepG2_m6A-seq_shWTAP_input_rep1.gtf
/data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shWTAP/HepG2_m6A-seq_shWTAP_input_rep2.gtf
/data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shWTAP/HepG2_m6A-seq_shWTAP_input_rep3.gtf
" > /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shWTAP/mergeGtfList.txt
echo -e "stringtie: start merging the assembled transcriptome...\n"
stringtie --merge -p 8 -G /data/zhoukr/reference/genome/gtf/gencode.v24lift37.annotation.gtf -o /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shWTAP_transcriptome_merge.gtf /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shWTAP/mergeGtfList.txt
echo -e "gffcompare: start comparing to the reference gtf annotation...\n"
gffcompare -r /data/zhoukr/reference/genome/gtf/gencode.v24lift37.annotation.gtf -o /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shWTAP/shWTAP_merged_compare /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/transcriptome_assembly/shWTAP_transcriptome_merge.gtf
echo -e "jobs have been done.\n"
