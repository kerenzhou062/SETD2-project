#!/public/zhoukr/softwares/R/bin/Rscript
library("exomePeak")

setwd("/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/Hela/R")
dir.create("./shM3")
setwd("./shM3")

gtf = "/data/zhoukr/reference/genome/gtf/gencode.v24lift37.annotation.gtf"
shContIP_rep1 = "/data/zhoukr/hhl_setd2_m6a/Hela_m6A-seq/shCont/Hela_m6A-seq_shCont_IP_rep1.fastq.sorted.bam"
shContIP_rep2 = "/data/zhoukr/hhl_setd2_m6a/Hela_m6A-seq/shCont/Hela_m6A-seq_shCont_IP_rep2.fastq.sorted.bam"
shContIP_rep3 = "/data/zhoukr/hhl_setd2_m6a/Hela_m6A-seq/shCont/Hela_m6A-seq_shCont_IP_rep3.fastq.sorted.bam"

shContInput_rep1 = "/data/zhoukr/hhl_setd2_m6a/Hela_m6A-seq/shCont/Hela_m6A-seq_shCont_input_rep1.fastq.sorted.bam"
shContInput_rep2 = "/data/zhoukr/hhl_setd2_m6a/Hela_m6A-seq/shCont/Hela_m6A-seq_shCont_input_rep2.fastq.sorted.bam"
shContInput_rep3 = "/data/zhoukr/hhl_setd2_m6a/Hela_m6A-seq/shCont/Hela_m6A-seq_shCont_input_rep3.fastq.sorted.bam"

shM3IP_rep1 = "/data/zhoukr/hhl_setd2_m6a/Hela_m6A-seq/shM3/Hela_m6A-seq_shM3_IP_rep1.fastq.sorted.bam"
shM3IP_rep2 = "/data/zhoukr/hhl_setd2_m6a/Hela_m6A-seq/shM3/Hela_m6A-seq_shM3_IP_rep2.fastq.sorted.bam"
shM3IP_rep3 = "/data/zhoukr/hhl_setd2_m6a/Hela_m6A-seq/shM3/Hela_m6A-seq_shM3_IP_rep3.fastq.sorted.bam"

shM3Input_rep1 = "/data/zhoukr/hhl_setd2_m6a/Hela_m6A-seq/shM3/Hela_m6A-seq_shM3_input_rep1.fastq.sorted.bam"
shM3Input_rep2 = "/data/zhoukr/hhl_setd2_m6a/Hela_m6A-seq/shM3/Hela_m6A-seq_shM3_input_rep2.fastq.sorted.bam"
shM3Input_rep3 = "/data/zhoukr/hhl_setd2_m6a/Hela_m6A-seq/shM3/Hela_m6A-seq_shM3_input_rep3.fastq.sorted.bam"

result = exomepeak(GENE_ANNO_GTF=gtf, IP_BAM=c(shContIP_rep1, shContIP_rep2, shContIP_rep3), INPUT_BAM=c(shContInput_rep1, shContInput_rep2, shContInput_rep3), TREATED_IP_BAM=c(shM3IP_rep1, shM3IP_rep2, shM3IP_rep3), TREATED_INPUT_BAM=c(shM3Input_rep1, shM3Input_rep2, shM3Input_rep3))

##rename and remove files
file.rename("./exomePeak_output/con_sig_diff_peak.bed", "./Hela_shM3_con_sig_diff_peak.bed")
file.rename("./exomePeak_output/con_sig_diff_peak.xls", "./Hela_shM3_con_sig_diff_peak.xls")
file.rename("./exomePeak_output/diff_peak.bed", "./Hela_shM3_diff_peak.bed")
file.rename("./exomePeak_output/diff_peak.xls", "./Hela_shM3_diff_peak.xls")
file.rename("./exomePeak_output/sig_diff_peak.bed", "./Hela_shM3_sig_diff_peak.bed")
file.rename("./exomePeak_output/sig_diff_peak.xls", "./Hela_shM3_sig_diff_peak.xls")
file.rename("./exomePeak_output/exomePeak.Rdata", "./Hela_shM3_exomePeak.Rdata")

file.remove("./exomePeak_output")
