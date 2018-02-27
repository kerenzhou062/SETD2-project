#!/public/zhoukr/softwares/R/bin/Rscript
library("exomePeak")

setwd("/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/R")
dir.create("./shSetD2")
setwd("./shSetD2")

gtf = "/data/zhoukr/reference/genome/gtf/gencode.v24lift37.annotation.gtf"
shContIP_rep1 = "/data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/shCont/HepG2_m6A-seq_shCont_IP_rep1.fastq.sorted.bam"
shContIP_rep2 = "/data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/shCont/HepG2_m6A-seq_shCont_IP_rep2.fastq.sorted.bam"
shContIP_rep3 = "/data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/shCont/HepG2_m6A-seq_shCont_IP_rep3.fastq.sorted.bam"

shContInput_rep1 = "/data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/shCont/HepG2_m6A-seq_shCont_input_rep1.fastq.sorted.bam"
shContInput_rep2 = "/data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/shCont/HepG2_m6A-seq_shCont_input_rep2.fastq.sorted.bam"
shContInput_rep3 = "/data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/shCont/HepG2_m6A-seq_shCont_input_rep3.fastq.sorted.bam"

shSetD2IP_rep1 = "/data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/shSetD2/HepG2_m6A-seq_shSetD2_IP_rep1.fastq.sorted.bam"
shSetD2IP_rep2 = "/data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/shSetD2/HepG2_m6A-seq_shSetD2_IP_rep2.fastq.sorted.bam"
shSetD2IP_rep3 = "/data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/shSetD2/HepG2_m6A-seq_shSetD2_IP_rep3.fastq.sorted.bam"

shSetD2Input_rep1 = "/data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/shSetD2/HepG2_m6A-seq_shSetD2_input_rep1.fastq.sorted.bam"
shSetD2Input_rep2 = "/data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/shSetD2/HepG2_m6A-seq_shSetD2_input_rep2.fastq.sorted.bam"
shSetD2Input_rep3 = "/data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/shSetD2/HepG2_m6A-seq_shSetD2_input_rep3.fastq.sorted.bam"

result = exomepeak(GENE_ANNO_GTF=gtf, IP_BAM=c(shContIP_rep1, shContIP_rep2, shContIP_rep3), INPUT_BAM=c(shContInput_rep1, shContInput_rep2, shContInput_rep3), TREATED_IP_BAM=c(shSetD2IP_rep1, shSetD2IP_rep2, shSetD2IP_rep3), TREATED_INPUT_BAM=c(shSetD2Input_rep1, shSetD2Input_rep2, shSetD2Input_rep3))

##rename and remove files
file.rename("./exomePeak_output/con_sig_diff_peak.bed", "./HepG2_shSetD2_con_sig_diff_peak.bed")
file.rename("./exomePeak_output/con_sig_diff_peak.xls", "./HepG2_shSetD2_con_sig_diff_peak.xls")
file.rename("./exomePeak_output/diff_peak.bed", "./HepG2_shSetD2_diff_peak.bed")
file.rename("./exomePeak_output/diff_peak.xls", "./HepG2_shSetD2_diff_peak.xls")
file.rename("./exomePeak_output/sig_diff_peak.bed", "./HepG2_shSetD2_sig_diff_peak.bed")
file.rename("./exomePeak_output/sig_diff_peak.xls", "./HepG2_shSetD2_sig_diff_peak.xls")
file.rename("./exomePeak_output/exomePeak.Rdata", "./HepG2_shSetD2_exomePeak.Rdata")

file.remove("./exomePeak_output")
