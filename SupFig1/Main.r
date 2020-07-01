#!/usr/bin/Rscript
library(ggpubr)
library(reshape2)	
#"#87CEEB", "#B03060"))
Main <- read.table("HorseMicrobiotaChara.status.txt",head=T)

 p = ggboxplot(Main,xlab="",ylab="Mapping Ratio","cohort","MappingRatio",add="jitter",color="cohort",ggtheme=theme_minimal())
 p = p + facet_wrap(~status,strip.position="bottom",scales="free_x")
G1 =  p +  theme(legend.position="none", axis.ticks.x =element_blank(),strip.text.x = element_text(colour = "black", size = 10))


 p = ggboxplot(Main,xlab="",ylab="Antibiotis-Resistant Gene Number","cohort","AntibiotisResistantGeneNumber",add="jitter",color="cohort",ggtheme=theme_minimal())
 p = p + facet_wrap(~status,strip.position="bottom",scales="free_x")
G2 =  p +  theme(legend.position="none", axis.ticks.x =element_blank(),strip.text.x = element_text(colour = "black", size = 10))



 p = ggboxplot(Main,xlab="",ylab="CAZY Gene Number","cohort","CAZYgeneNumber",add="jitter",color="cohort",ggtheme=theme_minimal())
 p = p + facet_wrap(~status,strip.position="bottom",scales="free_x")
G3 =  p +  theme(legend.position="none", axis.ticks.x =element_blank(),strip.text.x = element_text(colour = "black", size = 10))

 ggarrange(G1,G2,G3,nrow=3,labels=letters[1:3])
