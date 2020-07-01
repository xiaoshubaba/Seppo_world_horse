library(ggpubr)
library(reshape2)
CA = read.table("CAZY.families.txt",head=T)
ME = read.table("sample.metadata",head=T)
CA$status = ME$status
#color = colorRampPalette(c("white","#87CEEB", "#B03060"))
#CAZY_Figure_1
CAN = CA[,-which(colnames(CA)=="Others")]
p = ggboxplot(CA_L,palette=c("#B03060","#87CEEB"),"status","value",ylab="Relative abundance",add="jitter",color="status") + theme(legend.position="none")  + stat_compare_means(comparisons=list(c("domestic","wild")),method.args=list(alternative="less"),size=2.5)
facet(p,facet.by="variable",scales="free",strip.position="bottom")
