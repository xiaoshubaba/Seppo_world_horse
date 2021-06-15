library(RColorBrewer)
library(ggsci)
library(ggpubr)
library(pheatmap)
Metadata <- read.table("Wild.specific.feature.txt",head=T,sep="\t")
Metadata$Bar = abs(log(10,abs(Metadata$coefficient)))
Metadata$Feature = Metadata$specy	
Metadata$FeatureN = NA
i = 1
while (i<=length(Metadata$Bar)){
        if (Metadata$coefficient[i] < 0){
        Metadata$Bar[i] = Metadata$Bar[i] * -1
        }
        if (Metadata$qvalue[i] < 1e-3){
        Metadata$FeatureN[i] = paste0(Metadata$Feature[i],"  (+++)")
        }
        if( (Metadata$qvalue[i] > 1e-3) & (Metadata$qvalue[i] < 1e-2) ){
        Metadata$FeatureN[i] = paste0(Metadata$Feature[i],"  (++)")
        }
        if( Metadata$qvalue[i] > 1e-2 ){
        Metadata$FeatureN[i] = paste0(Metadata$Feature[i],"  (+)")
        }
i = i + 1
}
p = ggbarplot(Metadata,x="FeatureN",y="Bar",fill="cohort",color="white",sort.val="asc",sory.by.groups=FALSE,palette=c("#B03060","#86CEEB"),x.text.angle=90,rotate=TRUE,ylab="log(Coefficient)")+ theme(axis.text.y=element_text(size=6))  +theme(legend.position="right",axis.text.y=element_text(colour=c(rep("#B03060",65),rep("#86CEEB",16))))

MA = Metadata[order(Metadata$Bar,decreasing=TRUE),]
H <- data.frame(MA$Bar,MA$qvalue)
rownames(H) = MA$FeatureN
#
row_ann = MA$kingdom
row_ann = as.character(row_ann)
row_ann = factor(row_ann)
ROW = data.frame(category=row_ann)
#pal_d3("category10")(10)
Var1 <- c("#F8766D","gray","#00BFC4","#C77CFF") ; names(Var1) = c("Eukaryota","Bacteria","Virus","Archaea")	
ann_colors <- list(category=Var1)
rownames(ROW) = MA$FeatureN
pheatmap(H,fontsize_row=5,cluster_rows=FALSE,annotation_row=ROW,annotation_colors=ann_colors)
