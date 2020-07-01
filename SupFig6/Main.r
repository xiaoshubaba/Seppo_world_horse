library(ggpubr)
library(reshape2)
ME = read.table("sample.metadata",head=T)
TW = read.table("CeSta.matrix",head=T,sep="\t")
TW$cohort = ME$cohort
TW$status = ME$status
TWN = melt(TW,id=c("sample","cohort","status"))
p = ggboxplot(TWN,palette=c("#B03060","#87CEEB"),"status","value",ylab="Relative abundance",add="jitter",color="status") + theme(legend.position="none")+stat_compare_means(comparisons=list(c("domestic","wild")),method.args=list(alternative="less"))
facet(p,facet.by="variable",scales="free")
