#!/usr/bin/env Rscript
# exomepeak Script 2
# R script
# Define parameters and load library
library("ggplot2")
library("reshape2")
library("grid")

setwd("/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/Hela/bed12")
myData <- read.table("Hela_m6a_shCont_SetD2_dependent.gene",sep="\t",header = TRUE, stringsAsFactors=FALSE)
rowNum <- nrow(myData)
for (i in 1:rowNum) {
  if (myData[i,1] == 'lncRNA' ) {
    myData[i,1] = 'lncRNA'
  }else if (myData[i,1] == 'protein_coding' ) {
    myData[i,1] = 'protein coding'
  }else if (myData[i,1] == 'snRNA' ) {
    myData[i,1] = 'snRNA'
  }else if (myData[i,1] == 'snoRNA' ) {
    myData[i,1] = 'snoRNA'
  }else if (myData[i,1] == 'miRNA' ) {
    myData[i,1] = 'micro RNA'
  }else if (myData[i,1] == 'small ncRNA' ) {
    myData[i,1] = 'small ncRNA'
  }else if (myData[i,1] == 'pseudogene' ) {
    myData[i,1] = 'pseudogene'
  }else{
    myData[i,1] = 'others'
  }
}

Hela_m6a_shCont_SetD2_dependent_geneType_barPlot <- ggplot(data=myData, aes(x=geneType, y=peakNumber, fill=geneType)) + geom_bar(stat="identity")+
  theme_minimal() +
  labs(y = "Number of m6A peaks")+
  theme(panel.background = element_rect(fill = "white", colour = "black"), axis.title = element_text(size = 13, face = "bold"),
    axis.title.x=element_blank(), axis.text.x=element_blank())+
  theme(axis.text = element_text(size = 11), legend.text = element_text(size = 8), legend.title=element_blank())+
  scale_fill_brewer(palette="Dark2")

pdf("Hela_m6a_shCont_SetD2_dependent_geneType_barPlot.pdf")
plot(Hela_m6a_shCont_SetD2_dependent_geneType_barPlot)
dev.off()

myLabel = as.vector(myData$geneType)
myLabel = paste(myLabel, "(", round(myData$peakNumber / sum(myData$peakNumber) * 100, 2), "%)", sep = "")

Hela_m6a_shCont_SetD2_dependent_geneType_pieChart <- ggplot(data=myData, aes(x="", y=peakNumber, fill=geneType)) + geom_bar(stat="identity", width = 1)+
  coord_polar(theta = "y") + labs(x = "", y = "", title = "")+
  theme(panel.background = element_rect(fill = "white", colour = "black"))+
  theme(axis.ticks = element_blank()) +
  theme(legend.title = element_blank())+
  scale_fill_discrete(breaks = myData$geneType, labels = myLabel) +
  theme(axis.text.x = element_blank()) +
  theme(panel.grid=element_blank()) +
  theme(panel.border=element_blank())

pdf("Hela_m6a_shCont_SetD2_dependent_geneType_pieChart.pdf")
plot(Hela_m6a_shCont_SetD2_dependent_geneType_pieChart)
dev.off()


myData <- read.table("Hela_m6a_shCont_SetD2_dependent.region",sep="\t",header = TRUE, stringsAsFactors=FALSE)
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
myData$region <- factor(myData$region, levels=c("5' UTR", "CDS", "stop codon region", "3' UTR"))

Hela_m6a_shCont_SetD2_dependent_enrichment_barPlot <- ggplot(data=myData, aes(x=region, y=enrichment, fill=region)) + geom_bar(stat="identity")+theme_minimal() +
  labs(y = "Enrichment of m6A peaks")+
  theme(panel.background = element_rect(fill = "white", colour = "black"), axis.title = element_text(size = 13, face = "bold"),
    axis.title.x=element_blank(), axis.text.x=element_blank())+
  theme(axis.text = element_text(size = 11), legend.text = element_text(size = 8), legend.title=element_blank())+
  scale_fill_discrete(name="")+
  scale_fill_brewer(palette="Dark2")

pdf("Hela_m6a_shCont_SetD2_dependent_enrichment_barPlot.pdf")
plot(Hela_m6a_shCont_SetD2_dependent_enrichment_barPlot)
dev.off()



myLabel = as.vector(myData$region)
myLabel = paste(myLabel, "(", round(myData$peakNumber / sum(myData$peakNumber) * 100, 2), "%)", sep = "")

Hela_m6a_shCont_SetD2_dependent_region_pieChart <- ggplot(data=myData, aes(x="", y=peakNumber, fill=region)) + geom_bar(stat="identity", width = 1)+
  coord_polar(theta = "y") + labs(x = "", y = "", title = "")+
  theme(panel.background = element_rect(fill = "white", colour = "black"))+
  theme(axis.ticks = element_blank()) +
  theme(legend.title = element_blank())+
  scale_fill_discrete(breaks = myData$region, labels = myLabel) +
  theme(axis.text.x = element_blank()) +
  theme(panel.grid=element_blank()) +
  theme(panel.border=element_blank())


pdf("Hela_m6a_shCont_SetD2_dependent_region_pieChart.pdf")
plot(Hela_m6a_shCont_SetD2_dependent_region_pieChart)
dev.off()

setwd("/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/Hela/bed6")

myData <- read.table("Hela_m6a_shCont_SetD2_dependent.gene",sep="\t",header = TRUE, stringsAsFactors=FALSE)
rowNum <- nrow(myData)
for (i in 1:rowNum) {
  if (myData[i,1] == 'lncRNA' ) {
    myData[i,1] = 'lncRNA'
  }else if (myData[i,1] == 'protein_coding' ) {
    myData[i,1] = 'protein coding'
  }else if (myData[i,1] == 'snRNA' ) {
    myData[i,1] = 'snRNA'
  }else if (myData[i,1] == 'snoRNA' ) {
    myData[i,1] = 'snoRNA'
  }else if (myData[i,1] == 'miRNA' ) {
    myData[i,1] = 'micro RNA'
  }else if (myData[i,1] == 'small ncRNA' ) {
    myData[i,1] = 'small ncRNA'
  }else if (myData[i,1] == 'pseudogene' ) {
    myData[i,1] = 'pseudogene'
  }else{
    myData[i,1] = 'others'
  }
}

Hela_m6a_shCont_SetD2_dependent_geneType_barPlot <- ggplot(data=myData, aes(x=geneType, y=peakNumber, fill=geneType)) + geom_bar(stat="identity")+
  theme_minimal() +
  labs(y = "Number of m6A peaks")+
  theme(panel.background = element_rect(fill = "white", colour = "black"), axis.title = element_text(size = 13, face = "bold"),
    axis.title.x=element_blank(), axis.text.x=element_blank())+
  theme(axis.text = element_text(size = 11), legend.text = element_text(size = 8), legend.title=element_blank())+
  scale_fill_brewer(palette="Dark2")

pdf("Hela_m6a_shCont_SetD2_dependent_geneType_barPlot.pdf")
plot(Hela_m6a_shCont_SetD2_dependent_geneType_barPlot)
dev.off()

myLabel = as.vector(myData$geneType)
myLabel = paste(myLabel, "(", round(myData$peakNumber / sum(myData$peakNumber) * 100, 2), "%)", sep = "")

Hela_m6a_shCont_SetD2_dependent_geneType_pieChart <- ggplot(data=myData, aes(x="", y=peakNumber, fill=geneType)) + geom_bar(stat="identity", width = 1)+
  coord_polar(theta = "y") + labs(x = "", y = "", title = "")+
  theme(panel.background = element_rect(fill = "white", colour = "black"))+
  theme(axis.ticks = element_blank()) +
  theme(legend.title = element_blank())+
  scale_fill_discrete(breaks = myData$geneType, labels = myLabel) +
  theme(axis.text.x = element_blank()) +
  theme(panel.grid=element_blank()) +
  theme(panel.border=element_blank())

pdf("Hela_m6a_shCont_SetD2_dependent_geneType_pieChart.pdf")
plot(Hela_m6a_shCont_SetD2_dependent_geneType_pieChart)
dev.off()


myData <- read.table("Hela_m6a_shCont_SetD2_dependent.region",sep="\t",header = TRUE, stringsAsFactors=FALSE)
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
myData$region <- factor(myData$region, levels=c("5' UTR", "CDS", "stop codon region", "3' UTR"))

Hela_m6a_shCont_SetD2_dependent_enrichment_barPlot <- ggplot(data=myData, aes(x=region, y=enrichment, fill=region)) + geom_bar(stat="identity")+theme_minimal() +
  labs(y = "Enrichment of m6A peaks")+
  theme(panel.background = element_rect(fill = "white", colour = "black"), axis.title = element_text(size = 13, face = "bold"),
    axis.title.x=element_blank(), axis.text.x=element_blank())+
  theme(axis.text = element_text(size = 11), legend.text = element_text(size = 8), legend.title=element_blank())+
  scale_fill_discrete(name="")+
  scale_fill_brewer(palette="Dark2")

pdf("Hela_m6a_shCont_SetD2_dependent_enrichment_barPlot.pdf")
plot(Hela_m6a_shCont_SetD2_dependent_enrichment_barPlot)
dev.off()



myLabel = as.vector(myData$region)
myLabel = paste(myLabel, "(", round(myData$peakNumber / sum(myData$peakNumber) * 100, 2), "%)", sep = "")

Hela_m6a_shCont_SetD2_dependent_region_pieChart <- ggplot(data=myData, aes(x="", y=peakNumber, fill=region)) + geom_bar(stat="identity", width = 1)+
  coord_polar(theta = "y") + labs(x = "", y = "", title = "")+
  theme(panel.background = element_rect(fill = "white", colour = "black"))+
  theme(axis.ticks = element_blank()) +
  theme(legend.title = element_blank())+
  scale_fill_discrete(breaks = myData$region, labels = myLabel) +
  theme(axis.text.x = element_blank()) +
  theme(panel.grid=element_blank()) +
  theme(panel.border=element_blank())


pdf("Hela_m6a_shCont_SetD2_dependent_region_pieChart.pdf")
plot(Hela_m6a_shCont_SetD2_dependent_region_pieChart)
dev.off()
