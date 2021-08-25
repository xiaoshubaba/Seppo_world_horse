#!/usr/bin/Rscript
library(ggplot2)
library(ggpubr)
library(ade4)
args=commandArgs(T)
if (length(args) !=3){
stop("1) disance list ,head must be \"DistanceName Distance\". 2) Sample metadata. 3)Output file prefix")
}

#check if distance list exists
if (file.exists(args[1])){
disList <- read.table(args[1],head=T)
}else if (!file.exists(args[1])){
stop(paste0(args[1]," distance list not eixsts"))
}
# Check if head is right
if ( FALSE %in% (colnames(disList) == c("DistanceName","Distance"))){
stop(paste0(args[1]," head must be \"DistanceName Distance\""));
}
# CHeck if sample metadata exists
if (file.exists(args[2])){
sampleMeta <- read.table(args[2],head=T)
}else if (! file.exists(args[2])){
stop(paste0(args[2]," sample metadata not exist"))
}
#output
#PC corordinate
#	D1-P1	D1-P2	D2-P1	D2-P2
#S1
#S2
#S3
#PC x aixs and y axis
#	x	y	z	.	.
#D1
#D2
pcoaMerge <- list()
i = 1
GR <- list()
while (i <= dim(disList)[1]){
dis <- as.dist(read.table(as.character(disList$Distance[i])))
	if(!is.euclid(dis)){
	dis <- quasieuclid(dis)
	print(paste0("Warining:", disList$DistanceName[i]," distance is not euclid format"))
	}
pcoa <- dudi.pco(dis, scannf=F, nf=5)
PC_CONTU <- round(pcoa$eig/sum(pcoa$eig),digits=3)[1:5] * 100
PCs <- pcoa$li
#check if sample metadata and Pcs' sample name is constant
if (FALSE %in% (rownames(pcoa$li) == sampleMeta$sample)){
stop("pcoa sample name and cohort sample name not constant");
}
PCplot <- data.frame(sample=sampleMeta$sample,cohort=sampleMeta$cohort,PC1=pcoa$li[,1],PC2=pcoa$li[,2],status=sampleMeta$status)
GR[[i]] = ggscatter(PCplot,"PC1","PC2",color="status",shape="cohort",ggtheme=theme_minimal(),alpha=I(0.6),xlab=paste0("PC1(",PC_CONTU[1],"%)"),title=paste0(disList$DistanceName[i]),ylab=paste0("PC2(",PC_CONTU[2],"%)"),palette="jco") + theme(legend.position="none")
colnames(PCs) <- paste0(disList$DistanceName[i],".",colnames(PCs))
	if (class(pcoaMerge$PCs) == "NULL"){
	pcoaMerge$PCs = PCs
	pcoaMerge$Axis = PC_CONTU
	}else if (class(pcoaMerge$PCs) == "data.frame"){
	pcoaMerge$PCs = cbind(pcoaMerge$PCs,PCs)
	pcoaMerge$Axis = rbind(pcoaMerge$Axis,PC_CONTU)
	}
i = i + 1
}
#ggarrange(GR[[1]],GR[[2]],GR[[3]],GR[[4]],GR[[5]],GR[[6]],common.legend=TRUE,legend="bottom",labels=LETTERS[1:6])
ggarrange(GR[[1]],GR[[2]],GR[[3]],GR[[4]],GR[[5]],GR[[6]],GR[[7]],common.legend=TRUE,legend="bottom",labels=LETTERS[1:7])
rownames(pcoaMerge$Axis) = disList$DistanceName
write.table(pcoaMerge$Axis,paste0(args[3],".Axis"),quote=F,sep="\t")
write.table(pcoaMerge$PCs,paste0(args[3],".PCs"),quote=F,sep="\t")
save.image("test.R.image")
