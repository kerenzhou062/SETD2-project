#!/usr/bin/env Rscript
# exomepeak Script 2
# R script
# Define parameters and load library
library("ggplot2")
library("reshape2")
library("grid")

setwd("/data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/bed6")

myData <- read.table("HepG2_macs_shCont_SetD2_dependent.region",sep="\t",header = TRUE, stringsAsFactors=FALSE)
rowNum <- nrow(myData)
for (i in 1:rowNum) {
  if (myData[i,1] == '5utr' ) {
    myData[i,1] = "5' UTR"
  }else if (myData[i,1] == 'cds' ) {
    myData[i,1] = 'CDS'
  }else if (myData[i,1] == 'stopCodon' ) {
    myData[i,1] = 'stop codon region'
  }else if (myData[i,1] == '3utr' ) {
    myData[i,1] = "3' UTR"
  }
}
myData$region <- factor(myData$region, levels=c("promoter", "5' UTR", "CDS", "stop codon region", "3' UTR"))

HepG2_macs_shCont_SetD2_dependent_enrichment_barPlot <- ggplot(data=myData, aes(x=region, y=enrichment, fill=region)) + geom_bar(stat="identity")+theme_minimal() +
  labs(y = "Enrichment of H3K36me3 peaks")+
  theme(panel.background = element_rect(fill = "white", colour = "black"), axis.title = element_text(size = 13, face = "bold"),
    axis.title.x=element_blank(), axis.text.x=element_blank())+
  theme(axis.text = element_text(size = 11), legend.text = element_text(size = 8), legend.title=element_blank())+
  scale_fill_discrete(name="")+
  scale_fill_brewer(palette="Dark2")

pdf("HepG2_macs_shCont_SetD2_dependent_enrichment_barPlot.pdf")
plot(HepG2_macs_shCont_SetD2_dependent_enrichment_barPlot)
dev.off()



myLabel = as.vector(myData$region)
myLabel = paste(myLabel, "(", round(myData$peakNumber / sum(myData$peakNumber) * 100, 2), "%)", sep = "")

HepG2_macs_shCont_SetD2_dependent_region_pieChart <- ggplot(data=myData, aes(x="", y=peakNumber, fill=region)) + geom_bar(stat="identity", width = 1)+
  coord_polar(theta = "y") + labs(x = "", y = "", title = "")+
  theme(panel.background = element_rect(fill = "white", colour = "black"))+
  theme(axis.ticks = element_blank()) +
  theme(legend.title = element_blank())+
  scale_fill_discrete(breaks = myData$region, labels = myLabel) +
  theme(axis.text.x = element_blank()) +
  theme(panel.grid=element_blank()) +
  theme(panel.border=element_blank())


pdf("HepG2_macs_shCont_SetD2_dependent_region_pieChart.pdf")
plot(HepG2_macs_shCont_SetD2_dependent_region_pieChart)
dev.off()
