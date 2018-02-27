#!/usr/bin/env Rscript
# exomepeak Script 2
# R script
# Define parameters and load library
library("ggplot2")
library("reshape2")
library("grid")

setwd("/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/bed6")

utr5Text = textGrob(paste("5' UTR", sep = " "))
cdsText = textGrob(paste("CDS", sep = " "))
utr3Text = textGrob(paste("3' UTR", sep = " "))

myPeak <- read.table("HepG2_m6a.bed6bin",sep="\t",header = TRUE)
reshapeData <- melt(myPeak,id=c("region", "bin"))

levels(reshapeData$variable)[levels(reshapeData$variable)=="shCont"] <- "shNS"
levels(reshapeData$variable)[levels(reshapeData$variable)=="shSetD2"] <- "shSETD2"
levels(reshapeData$variable)[levels(reshapeData$variable)=="shM14"] <- "shMETTL14"
levels(reshapeData$variable)[levels(reshapeData$variable)=="shM3"] <- "shMETTL3"
levels(reshapeData$variable)[levels(reshapeData$variable)=="shWTAP"] <- "shWTAP"

annotation_position = -17
HepG2_m6a_all_bin <- ggplot(data=reshapeData, aes(x=bin, y=value)) +
  geom_line(aes(colour=variable), size=1.5) + theme_bw() +
  geom_vline(aes(xintercept=100), colour="#BB0000", linetype="dashed") +
  geom_vline(aes(xintercept=200), colour="#BB0000", linetype="dashed") + labs(y = "Number of m6A peaks") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.title = element_text(size = 13, face = "bold"), axis.title.x=element_blank(),axis.text.x=element_blank())+
  theme(axis.text = element_text(size = 11, face = "bold"), legend.text = element_text(size = 9, face = "bold"), legend.title=element_blank(),)+
  theme(legend.key.width=unit(1,"line"),legend.key.height=unit(1,"line"), legend.justification=c(1,0), legend.position=c(0.2,0.7))+
  theme(plot.margin = unit(c(1,1,3,1), "cm"))+
  scale_y_continuous(limits=c(0,250)) +
  annotation_custom(grob = utr5Text, xmin = 50, xmax = 50, ymin = annotation_position, ymax = annotation_position)+
  annotation_custom(grob = cdsText, xmin = 150, xmax = 150, ymin = annotation_position, ymax = annotation_position)+
  annotation_custom(grob = utr3Text, xmin = 250, xmax = 250, ymin = annotation_position, ymax = annotation_position)

HepG2_m6a_all_bin <- ggplotGrob(HepG2_m6a_all_bin)
HepG2_m6a_all_bin$layout$clip[HepG2_m6a_all_bin$layout$name=="panel"] <- "off"
grid.draw(HepG2_m6a_all_bin)

pdf("HepG2_m6a_all_bin.pdf", width=7, height=7)
plot(HepG2_m6a_all_bin)
dev.off()

myPeak <- read.table("HepG2_m6a_M14_dependent.bed6bin",sep="\t",header = TRUE)
reshapeData <- melt(myPeak,id=c("region", "bin"))

levels(reshapeData$variable)[levels(reshapeData$variable)=="shCont"] <- "shNS"
levels(reshapeData$variable)[levels(reshapeData$variable)=="shSetD2"] <- "shSETD2"

annotation_position = -8
HepG2_m6a_M14_dependent_bin <- ggplot(data=reshapeData, aes(x=bin, y=value)) +
  geom_line(aes(colour=variable), size=1.5) +
  scale_color_manual(values=c("#1b6ab1", "#ca5713")) + theme_bw() +
  geom_vline(aes(xintercept=100), colour="#BB0000", linetype="dashed") +
  geom_vline(aes(xintercept=200), colour="#BB0000", linetype="dashed") + labs(y = "Number of m6A peaks") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.title = element_text(size = 13, face = "bold"), axis.title.x=element_blank(),axis.text.x=element_blank())+
  theme(axis.text = element_text(size = 11, face = "bold"), legend.text = element_text(size = 9, face = "bold"), legend.title=element_blank(),)+
  theme(legend.key.width=unit(1.5,"line"),legend.key.height=unit(1.5,"line"), legend.justification=c(1,0), legend.position=c(0.2,0.7))+
  theme(plot.margin = unit(c(1,1,3,1), "cm"))+
  scale_y_continuous(limits=c(0,110)) +
  annotation_custom(grob = utr5Text, xmin = 50, xmax = 50, ymin = annotation_position, ymax = annotation_position)+
  annotation_custom(grob = cdsText, xmin = 150, xmax = 150, ymin = annotation_position, ymax = annotation_position)+
  annotation_custom(grob = utr3Text, xmin = 250, xmax = 250, ymin = annotation_position, ymax = annotation_position)

HepG2_m6a_M14_dependent_bin <- ggplotGrob(HepG2_m6a_M14_dependent_bin)
HepG2_m6a_M14_dependent_bin$layout$clip[HepG2_m6a_M14_dependent_bin$layout$name=="panel"] <- "off"
grid.draw(HepG2_m6a_M14_dependent_bin)

pdf("HepG2_m6a_M14_dependent_bin.pdf", width=7, height=7)
plot(HepG2_m6a_M14_dependent_bin)
dev.off()

myPeak <- read.table("HepG2_m6a_M3_dependent.bed6bin",sep="\t",header = TRUE)
reshapeData <- melt(myPeak,id=c("region", "bin"))

levels(reshapeData$variable)[levels(reshapeData$variable)=="shCont"] <- "shNS"
levels(reshapeData$variable)[levels(reshapeData$variable)=="shSetD2"] <- "shSETD2"

HepG2_m6a_M3_dependent_bin <- ggplot(data=reshapeData, aes(x=bin, y=value)) +
  geom_line(aes(colour=variable), size=1.5) +
  scale_color_manual(values=c("#1b6ab1", "#ca5713")) + theme_bw() +
  geom_vline(aes(xintercept=100), colour="#BB0000", linetype="dashed") +
  geom_vline(aes(xintercept=200), colour="#BB0000", linetype="dashed") + labs(y = "Number of m6A peaks") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.title = element_text(size = 13, face = "bold"), axis.title.x=element_blank(),axis.text.x=element_blank())+
  theme(axis.text = element_text(size = 11, face = "bold"), legend.text = element_text(size = 9, face = "bold"), legend.title=element_blank(),)+
  theme(legend.key.width=unit(1.5,"line"),legend.key.height=unit(1.5,"line"), legend.justification=c(1,0), legend.position=c(0.2,0.7))+
  theme(plot.margin = unit(c(1,1,3,1), "cm"))+
  scale_y_continuous(limits=c(0,110)) +
  annotation_custom(grob = utr5Text, xmin = 50, xmax = 50, ymin = annotation_position, ymax = annotation_position)+
  annotation_custom(grob = cdsText, xmin = 150, xmax = 150, ymin = annotation_position, ymax = annotation_position)+
  annotation_custom(grob = utr3Text, xmin = 250, xmax = 250, ymin = annotation_position, ymax = annotation_position)

HepG2_m6a_M3_dependent_bin <- ggplotGrob(HepG2_m6a_M3_dependent_bin)
HepG2_m6a_M3_dependent_bin$layout$clip[HepG2_m6a_M3_dependent_bin$layout$name=="panel"] <- "off"
grid.draw(HepG2_m6a_M3_dependent_bin)

ggarrange(HepG2_m6a_M3_dependent_bin, HepG2_m6a_M3_dependent_bin)
pdf("HepG2_m6a_M3_dependent_bin.pdf", width=7, height=7)
plot(HepG2_m6a_M3_dependent_bin)
dev.off()

myPeak <- read.table("HepG2_m6a_WTAP_dependent.bed6bin",sep="\t",header = TRUE)
reshapeData <- melt(myPeak,id=c("region", "bin"))

levels(reshapeData$variable)[levels(reshapeData$variable)=="shCont"] <- "shNS"
levels(reshapeData$variable)[levels(reshapeData$variable)=="shSetD2"] <- "shSETD2"

annotation_position = -11
HepG2_m6a_WTAP_dependent_bin <- ggplot(data=reshapeData, aes(x=bin, y=value)) +
  geom_line(aes(colour=variable), size=1.5) +
  scale_color_manual(values=c("#1b6ab1", "#ca5713")) + theme_bw() +
  geom_vline(aes(xintercept=100), colour="#BB0000", linetype="dashed") +
  geom_vline(aes(xintercept=200), colour="#BB0000", linetype="dashed") + labs(y = "Number of m6A peaks") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.title = element_text(size = 13, face = "bold"), axis.title.x=element_blank(),axis.text.x=element_blank())+
  theme(axis.text = element_text(size = 11, face = "bold"), legend.text = element_text(size = 9, face = "bold"), legend.title=element_blank(),)+
  theme(legend.key.width=unit(1.5,"line"),legend.key.height=unit(1.5,"line"), legend.justification=c(1,0), legend.position=c(0.2,0.7))+
  theme(plot.margin = unit(c(1,1,3,1), "cm"))+
  scale_y_continuous(limits=c(0,150)) +
  annotation_custom(grob = utr5Text, xmin = 50, xmax = 50, ymin = annotation_position, ymax = annotation_position)+
  annotation_custom(grob = cdsText, xmin = 150, xmax = 150, ymin = annotation_position, ymax = annotation_position)+
  annotation_custom(grob = utr3Text, xmin = 250, xmax = 250, ymin = annotation_position, ymax = annotation_position)

HepG2_m6a_WTAP_dependent_bin <- ggplotGrob(HepG2_m6a_WTAP_dependent_bin)
HepG2_m6a_WTAP_dependent_bin$layout$clip[HepG2_m6a_WTAP_dependent_bin$layout$name=="panel"] <- "off"
grid.draw(HepG2_m6a_WTAP_dependent_bin)

pdf("HepG2_m6a_WTAP_dependent_bin.pdf", width=7, height=7)
plot(HepG2_m6a_WTAP_dependent_bin)
dev.off()

myPeak <- read.table("HepG2_m6a_nonWriter_dependent.bed6bin",sep="\t",header = TRUE)
reshapeData <- melt(myPeak,id=c("region", "bin"))

levels(reshapeData$variable)[levels(reshapeData$variable)=="shCont"] <- "shNS"
levels(reshapeData$variable)[levels(reshapeData$variable)=="shSetD2"] <- "shSETD2"

annotation_position = -8
HepG2_m6a_nonWriter_dependent_bin <- ggplot(data=reshapeData, aes(x=bin, y=value)) +
  geom_line(aes(colour=variable), size=1.5) +
  scale_color_manual(values=c("#1b6ab1", "#ca5713")) + theme_bw() +
  geom_vline(aes(xintercept=100), colour="#BB0000", linetype="dashed") +
  geom_vline(aes(xintercept=200), colour="#BB0000", linetype="dashed") + labs(y = "Number of m6A peaks") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.title = element_text(size = 13, face = "bold"), axis.title.x=element_blank(),axis.text.x=element_blank())+
  theme(axis.text = element_text(size = 11, face = "bold"), legend.text = element_text(size = 9, face = "bold"), legend.title=element_blank(),)+
  theme(legend.key.width=unit(1.5,"line"),legend.key.height=unit(1.5,"line"), legend.justification=c(1,0), legend.position=c(0.2,0.7))+
  theme(plot.margin = unit(c(1,1,3,1), "cm"))+
  scale_y_continuous(limits=c(0,110)) +
  annotation_custom(grob = utr5Text, xmin = 50, xmax = 50, ymin = annotation_position, ymax = annotation_position)+
  annotation_custom(grob = cdsText, xmin = 150, xmax = 150, ymin = annotation_position, ymax = annotation_position)+
  annotation_custom(grob = utr3Text, xmin = 250, xmax = 250, ymin = annotation_position, ymax = annotation_position)

HepG2_m6a_nonWriter_dependent_bin <- ggplotGrob(HepG2_m6a_nonWriter_dependent_bin)
HepG2_m6a_nonWriter_dependent_bin$layout$clip[HepG2_m6a_nonWriter_dependent_bin$layout$name=="panel"] <- "off"
grid.draw(HepG2_m6a_nonWriter_dependent_bin)

pdf("HepG2_m6a_nonWriter_dependent_bin.pdf", width=7, height=7)
plot(HepG2_m6a_nonWriter_dependent_bin)
dev.off()
