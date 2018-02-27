#!/public/zhoukr/softwares/R/bin/Rscript
library(Guitar)
setwd("/data/zhoukr/hhl_setd2_m6a/analysis/")
txdb_file = "/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.sqlite"
hg19_toy_txdb = loadDb(txdb_file)
gc_hg19_txdb = makeGuitarCoordsFromTxDb(hg19_toy_txdb, noBins=100)

shCont_M3_m6a_bed12 = "./m6a-seq/HepG2/HepG2_m6a_shCont_M3_dependent.bed12"
shSetD2_M3_m6a_bed12 = "./m6a-seq/HepG2/HepG2_m6a_shSetD2_M3_dependent.bed12"
shCont_M3_m6a = BED12toGRangesList(shCont_M3_m6a_bed12)
shSetD2_M3_m6a = BED12toGRangesList(shSetD2_M3_m6a_bed12)
feature_M3 = list(shCont_M3_m6a, shSetD2_M3_m6a)
names(feature_M3) = c("shCont_M3", "shSetD2_M3")

GuitarPlot(gfeatures = feature_M3, GuitarCoordsFromTxDb = gc_hg19_txdb, rescaleComponent=FALSE, saveToPDFprefix = "HepG2_m6a_M3")

shCont_m6a_bed12 = "./m6a-seq/HepG2/HepG2_shCont_m6a.bed12"
shSetD2_m6a_bed12 = "./m6a-seq/HepG2/HepG2_shSetD2_m6a.bed12"
shCont_m6a = BED12toGRangesList(shCont_m6a_bed12)
shSetD2_m6a = BED12toGRangesList(shSetD2_m6a_bed12)
feature_HepG2 = list(shCont_m6a, shSetD2_m6a)
names(feature_HepG2) = c("shCont", "shSetD2")
GuitarPlot(gfeatures = feature_HepG2,GuitarCoordsFromTxDb = gc_hg19_txdb, rescaleComponent=FALSE, saveToPDFprefix = "HepG2_m6a")

shCont_M14_m6a_bed12 = "./m6a-seq/HepG2/HepG2_m6a_shCont_M14_dependent.bed12"
shSetD2_M14_m6a_bed12 = "./m6a-seq/HepG2/HepG2_m6a_shSetD2_M14_dependent.bed12"
shCont_M14_m6a = BED12toGRangesList(shCont_M14_m6a_bed12)
shSetD2_M14_m6a = BED12toGRangesList(shSetD2_M14_m6a_bed12)
feature_M14 = list(shCont_M14_m6a, shSetD2_M14_m6a)
names(feature_M14) = c("shCont_M14", "shSetD2_M14")

feature_Combine = list(shCont_M3_m6a, shSetD2_M3_m6a, shCont_M14_m6a, shSetD2_M14_m6a)
names(feature_Combine) = c("shCont_M3", "shSetD2_M3", "shCont_M14", "shSetD2_M14")

GuitarPlot(gfeatures = feature_Combine, GuitarCoordsFromTxDb = gc_hg19_txdb, rescaleComponent=FALSE, saveToPDFprefix = "HepG2_m6a_M3_M4")
