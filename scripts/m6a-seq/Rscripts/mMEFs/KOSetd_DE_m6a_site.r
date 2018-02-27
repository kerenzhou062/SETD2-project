#!/public/zhoukr/softwares/R/bin/Rscript
library("exomePeak")

setwd("/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/mMEFs")
dir.create("./KO")
setwd("./KO")

gtf = "/data/zhoukr/reference/genome/gtf/gencode.v24lift37.annotation.gtf"
f1 = "/data/zhoukr/hhl_setd2_m6a/mouse_MEFs_m6A-seq/SetD2-WT_NA/mMEFs_m6A-seq_SetD2-WT_NA_IP.sorted.bam"
f2 = "/data/zhoukr/hhl_setd2_m6a/mouse_MEFs_m6A-seq/SetD2-WT_NA/mMEFs_m6A-seq_SetD2-WT_NA_input.sorted.bam"
f3 = "/data/zhoukr/hhl_setd2_m6a/mouse_MEFs_m6A-seq/SetD2-KO_NA/mMEFs_m6A-seq_SetD2-KO_NA_IP.sorted.bam"
f4 = "/data/zhoukr/hhl_setd2_m6a/mouse_MEFs_m6A-seq/SetD2-KO_NA/mMEFs_m6A-seq_SetD2-KO_NA_input.sorted.bam"

result = exomepeak(GENE_ANNO_GTF=gtf, IP_BAM=f1, INPUT_BAM=f2, TREATED_IP_BAM=f3, TREATED_INPUT_BAM=f4)

##rename and remove files
file.rename("./exomePeak_output/con_sig_diff_peak.bed", "./../mMEFs_KOSetD2_diff_peak.bed")
file.rename("./exomePeak_output/con_sig_diff_peak.xls", "./../mMEFs_KOSetD2_diff_peak.xls")

file.remove("./exomePeak_output/exomePeak.Rdata", "./exomePeak_output/diff_peak.bed", "./exomePeak_output/sig_diff_peak.bed",
	"./exomePeak_output/diff_peak.xls", "./exomePeak_output/sig_diff_peak.xls", "./exomePeak_output")
setwd("./../")
file.remove("./KO")
