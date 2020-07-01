#!/usr/bin/Rscript
library(ggpubr)
library(reshape2)	
#"#87CEEB", "#B03060"))
Main <- read.table("HorseMicrobiotaChara.status.txt",head=T)
Main_L <- melt(Main,id=c("sample","cohort","status"))
Main_L$variable = factor(Main_L$variable,levels=c("MappingRatio","ObservedSpecies","GeneRichness","AntibiotisResistantGeneNumber","CAZYgeneNumber"))
p = ggboxplot(Main_L,ylab="Microbial characteristics","status","value",add="jitter",color="status",palette=c("#B03060","#87CEEB")) + stat_compare_means(comparisons=list(c("domestic","wild")),method.args=list(alternative="less"),method="t.test",size=2.8)
p = p + facet_wrap(~variable,strip.position="bottom",scales="free_y")

p +  theme(legend.position="none", axis.ticks.x =element_blank(), axis.text.x=element_blank(),strip.text.x = element_text(colour = "black", size = 6))
