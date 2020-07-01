#!/usr/bin/Rscript
library(ggplot2)
library(ggpubr)
library(ade4)

sampleMeta <- read.table("sample.country.metadata",head=T)

#species
dis <- as.dist(read.table("whole.spearman.dist"))
dis <- quasieuclid(dis)
pcoa <- dudi.pco(dis, scannf=F, nf=5)
PC_CONTU <- round(pcoa$eig/sum(pcoa$eig),digits=3)[1:5] * 100
PCs <- pcoa$li
#check if sample metadata and Pcs' sample name is constant
if (FALSE %in% (rownames(pcoa$li) == sampleMeta$sample)){
stop("pcoa sample name and cohort sample name not constant");
}
PCplot <- data.frame(sample=sampleMeta$sample,status=sampleMeta$status,country=sampleMeta$country,PC1=pcoa$li[,1],PC2=pcoa$li[,2])
WholePCOA = ggscatter(PCplot,size=3,"PC1","PC2",shape="country",color="status",alpha=I(0.6),xlab=paste0("PCO1(",PC_CONTU[1],"%)"),title="species",ylab=paste0("PCO2(",PC_CONTU[2],"%)")) + theme(legend.position="none")

#CAZY
dis <- as.dist(read.table("EC.spearman.dist"))
dis <- quasieuclid(dis)
pcoa <- dudi.pco(dis, scannf=F, nf=5)
PC_CONTU <- round(pcoa$eig/sum(pcoa$eig),digits=3)[1:5] * 100
PCs <- pcoa$li
#check if sample metadata and Pcs' sample name is constant
if (FALSE %in% (rownames(pcoa$li) == sampleMeta$sample)){
stop("pcoa sample name and cohort sample name not constant");
}
PCplot <- data.frame(sample=sampleMeta$sample,status=sampleMeta$status,country=sampleMeta$country,PC1=pcoa$li[,1],PC2=pcoa$li[,2])
ECPCOA = ggscatter(PCplot,size=3,"PC1","PC2",shape="country",color="status",alpha=I(0.6),xlab=paste0("PCO1(",PC_CONTU[1],"%)"),title="CAZY genes",ylab=paste0("PCO2(",PC_CONTU[2],"%)")) + theme(legend.position="none")

ggarrange(WholePCOA,ECPCOA,common.legend=TRUE,legend="right",labels=letters[1:2])
