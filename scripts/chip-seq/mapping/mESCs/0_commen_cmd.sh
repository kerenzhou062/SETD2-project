
# 测序质量检测
nohup fastqc *.fastq >1_fastqc_runlog.txt 2>&1 &

mkdir 0_fastqc_output
mv *.html ./0_fastqc_output/
mv *.zip ./0_fastqc_output/

# 运行下面代码块，用以生成 2_do_bowtie_and_sam_to_sort_bam.sh 脚本
for i in `ls -1 *.fastq`
do 
	echo 'echo '$i;
	echo 'bowtie -t -p 16 -v 2 -m 1 --best --strata --sam /public/genomes/mouse/mm10/bowtieIndex_UCSC/mm10 '$i' '${i%.fastq}.sam;
	echo 'samtools view -h -bS -F 4 --threads 32 '${i%.fastq}.sam' -o '${i%.fastq}.bam;
	echo 'rm '${i%.fastq}.sam;
	echo 'samtools sort --threads 16 -m 2G -O bam -o '${i%.fastq}.sorted.bam' '${i%.fastq}.bam
	echo 'rm '${i%.fastq}.bam;
	echo 'samtools index -b '${i%.fastq}.sorted.bam;
	echo ''
done >2_do_bowtie_and_sam_to_sort_bam.sh

# 比对并输出排序的 bam及其索引文件
nohup bash 2_do_bowtie_and_sam_to_sort_bam.sh >2_do_bowtie_and_sam_to_sort_bam_runlog.txt 2>&1 &


mkdir Ctrl_D0
mkdir SetD2-KD_D0

mv *Ctrl_D0*.sorted.bam ./Ctrl_D0/
mv *SetD2-KD_D0*.sorted.bam ./SetD2-KD_D0/

mv *Ctrl_D0*.sorted.bam.bai ./Ctrl_D0/
mv *SetD2-KD_D0*.sorted.bam.bai ./SetD2-KD_D0/


# call peak
cd ./Ctrl_D0/
nohup macs14 -t mES_ChIP-seq_Ctrl_D0_ChIP-rep1.sorted.bam mES_ChIP-seq_Ctrl_D0_ChIP-rep2.sorted.bam -c mES_ChIP-seq_Ctrl_D0_input-rep1.sorted.bam mES_ChIP-seq_Ctrl_D0_input-rep2.sorted.bam --name=Ctrl_D0_macs --format="BAM" --g mm --tsize=50 --nomodel --shiftsize=125 -w -S --call-subpeaks >macs_runlog_for_Ctrl_D0.txt 2>&1 &
cd ../

cd ./SetD2-KD_D0/
nohup macs14 -t mES_ChIP-seq_SetD2-KD_D0_input-rep1.sorted.bam mES_ChIP-seq_SetD2-KD_D0_input-rep2.sorted.bam -c mES_ChIP-seq_SetD2-KD_D0_ChIP-rep2.sorted.bam mES_ChIP-seq_SetD2-KD_D0_ChIP-rep1.sorted.bam --name=SetD2-KD_D0_macs --format="BAM" --g mm --tsize=50 --nomodel --shiftsize=125 -w -S --call-subpeaks >macs_runlog_for_SetD2-KD_D0.txt 2>&1 &
cd ../






