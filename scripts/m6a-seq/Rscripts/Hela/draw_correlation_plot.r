#!/usr/bin/env Rscript
# exomepeak Script 2
# R script
# Define parameters and load library
library("ggplot2")
library("corrplot")

setwd("/data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/Hela/xls/overlappedPeak")
myData <- read.table("Hela_m6a_peak_FC_log2.txt",sep="\t",header = TRUE)
log2FC <- myData[14:17]
M <- cor(log2FC)
col1 <- colorRampPalette(c("#00007F","blue","#007FFF","cyan","white", "yellow", "#FF7F00", "red","#7F0000"))
corrplot(M, method="color", col=col1(500), cl.lim=c(0,1), order = "AOE", addCoef.col="black")

panel.smooth <- function (x, y) {
  points(x, y,col="blue")
  abline(lm(y~x), col="red",lwd=2)
}
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
pairs(log2FC,lower.panel=panel.smooth, upper.panel=panel.cor)

pValFunction <- function(pval) {
  if (pval == 0) {
    return ("< 2.22e-308")
  }else if (pval == 1) {
    return ("= 1")
  }else{
    exactVal = sprintf("%.2e", pval)
    exactVal = paste("= ", exactVal, sep = "")
    return (exactVal)
  }
}

shM3 <- log2FC$shM3
shM14 <- log2FC$shM14
shWTAP <- log2FC$shWTAP
shSetD2 <- log2FC$shSetD2
cor_test_shM3 <- cor.test(shM3,shSetD2,method="pearson",exact=TRUE)
cor_test_shM14 <- cor.test(shM14,shSetD2,method="pearson",exact=TRUE)
cor_test_shWTAP <- cor.test(shWTAP,shSetD2,method="pearson",exact=TRUE)
cor_test_shM3_cor = sprintf("%.2f",cor_test_shM3$estimate)
cor_test_shM14_cor = sprintf("%.2f",cor_test_shM14$estimate)
cor_test_shWTAP_cor = sprintf("%.2f",cor_test_shWTAP$estimate)
cor_test_shM3_pval = pValFunction(cor_test_shM3$p.value)
cor_test_shM14_pval = pValFunction(cor_test_shM14$p.value)
cor_test_shWTAP_pval = pValFunction(cor_test_shWTAP$p.value)

p1 <- ggplot(log2FC, aes(x=shSetD2, y=shM14)) +
  geom_point(colour="blue") + geom_smooth(method=lm,se=FALSE,colour="red",size = 2) +
  annotate("text", x = 2.37, y = -4.7, label = paste("r = ", cor_test_shM14_cor, sep = ""),fontface="bold") +
  annotate("text", x = 2.7, y = -5.0, label = paste("p ", cor_test_shM14_pval, sep = ""),fontface="bold") +
  theme_bw() + labs(x="shSETD2 FC(log2)", y = "shMETTL14 FC(log2)") +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 13, face = "bold"))

p2 <- ggplot(log2FC, aes(x=shSetD2, y=shM3)) +
  geom_point(colour="blue") + geom_smooth(method=lm,se=FALSE,colour="red",size = 2) +
  annotate("text", x = 2.37, y = -4.7, label = paste("r = ", cor_test_shM3_cor, sep = ""),fontface="bold") +
  annotate("text", x = 2.7, y = -5.0, label = paste("p ", cor_test_shM3_pval, sep = ""),fontface="bold") +
  theme_bw() + labs(x="shSETD2 FC(log2)", y = "shMETTL3 FC(log2)") +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 13, face = "bold"))

p3 <- ggplot(log2FC, aes(x=shSetD2, y=shWTAP)) +
  geom_point(colour="blue") + geom_smooth(method=lm,se=FALSE,colour="red",size = 2) +
  annotate("text", x = 2.37, y = -4.7, label = paste("r = ", cor_test_shWTAP_cor, sep = ""),fontface="bold") +
  annotate("text", x = 2.7, y = -5.0, label = paste("p ", cor_test_shWTAP_pval, sep = ""),fontface="bold") +
  theme_bw() + labs(x="shSETD2 FC(log2)", y = "shWTAP FC(log2)") +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 13, face = "bold"))

library("gridExtra")
p4 <- grid.arrange(p1, p2, p3, ncol=3, nrow=2, widths=c(1,1,1), heights =c(1,1))

pdf("Hela_FC_correlation_shM14vs.shSetD2.pdf")
plot(p1)
dev.off()

pdf("Hela_FC_correlation_shM3vs.shSetD2.pdf")
plot(p2)
dev.off()

pdf("Hela_FC_correlation_shWTAPvs.shSetD2.pdf")
plot(p3)
dev.off()


pdf("Hela_FC_correlation_allvs.shSetD2.pdf")
plot(p4)
dev.off()


myData <- read.table("Hela_m6a_peak_FoldEnrichment_log2.txt",sep="\t",header = TRUE)
log2Erichment <- myData[14:17]
M <- cor(log2Erichment)
col1 <- colorRampPalette(c("#00007F","blue","#007FFF","cyan","white", "yellow", "#FF7F00", "red","#7F0000"))
corrplot(M, method="color", col=col1(500), cl.lim=c(0,1), order = "AOE", addCoef.col="black")

panel.smooth <- function (x, y) {
  points(x, y,col="blue")
  abline(lm(y~x), col="red",lwd=2)
}
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
pairs(log2Erichment,lower.panel=panel.smooth, upper.panel=panel.cor)



shM3 <- log2Erichment$shM3
shM14 <- log2Erichment$shM14
shWTAP <- log2Erichment$shWTAP
shSetD2 <- log2Erichment$shSetD2
cor_test_shM3 <- cor.test(shM3,shSetD2,method="pearson",exact=TRUE)
cor_test_shM14 <- cor.test(shM14,shSetD2,method="pearson",exact=TRUE)
cor_test_shWTAP <- cor.test(shWTAP,shSetD2,method="pearson",exact=TRUE)
cor_test_shM3_cor = sprintf("%.2f",cor_test_shM3$estimate)
cor_test_shM14_cor = sprintf("%.2f",cor_test_shM14$estimate)
cor_test_shWTAP_cor = sprintf("%.2f",cor_test_shWTAP$estimate)
cor_test_shM3_pval = pValFunction(cor_test_shM3$p.value)
cor_test_shM14_pval = pValFunction(cor_test_shM14$p.value)
cor_test_shWTAP_pval = pValFunction(cor_test_shWTAP$p.value)

p1 <- ggplot(log2Erichment, aes(x=shSetD2, y=shM14)) +
  geom_point(colour="blue") + geom_smooth(method=lm,se=FALSE,colour="red",size = 2) +
  annotate("text", x = 5.87, y = 0.3, label = paste("r = ", cor_test_shM14_cor, sep = ""),fontface="bold") +
  annotate("text", x = 6.1, y = 0, label = paste("p ", cor_test_shM14_pval, sep = ""),fontface="bold") +
  theme_bw() + labs(x="shSETD2 Erichment(log2)", y = "shMETTL14 Erichment(log2)") +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 13, face = "bold"))

p2 <- ggplot(log2Erichment, aes(x=shSetD2, y=shM3)) +
  geom_point(colour="blue") + geom_smooth(method=lm,se=FALSE,colour="red",size = 2) +
  annotate("text", x = 5.87, y = 0.3, label = paste("r = ", cor_test_shM3_cor, sep = ""),fontface="bold") +
  annotate("text", x = 6.1, y = 0, label = paste("p ", cor_test_shM3_pval, sep = ""),fontface="bold") +
  theme_bw() + labs(x="shSETD2 Erichment(log2)", y = "shMETTL3 Erichment(log2)") +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 13, face = "bold"))

p3 <- ggplot(log2Erichment, aes(x=shSetD2, y=shWTAP)) +
  geom_point(colour="blue") + geom_smooth(method=lm,se=FALSE,colour="red",size = 2) +
  annotate("text", x = 5.87, y = 0.3, label = paste("r = ", cor_test_shWTAP_cor, sep = ""),fontface="bold") +
  annotate("text", x = 6.1, y = 0, label = paste("p ", cor_test_shWTAP_pval, sep = ""),fontface="bold") +
  theme_bw() + labs(x="shSETD2 Erichment(log2)", y = "shWTAP Erichment(log2)") +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 13, face = "bold"))

library("gridExtra")
p4 <- grid.arrange(p1, p2, p3, ncol=3, nrow=2, widths=c(1,1,1), heights =c(1,1))

pdf("Hela_Erichment_correlation_shM14vs.shSetD2.pdf")
plot(p1)
dev.off()

pdf("Hela_Erichment_correlation_shM3vs.shSetD2.pdf")
plot(p2)
dev.off()

pdf("Hela_Erichment_correlation_shWTAPvs.shSetD2.pdf")
plot(p3)
dev.off()


pdf("Hela_Erichment_correlation_allvs.shSetD2.pdf")
plot(p4)
dev.off()
