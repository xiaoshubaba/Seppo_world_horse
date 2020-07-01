library(ggpubr)
K = read.table("kingdom.ratio",head=T)
K$sample = factor(K$sample,levels=K$sample)
cK_L = melt(K,id=c("sample","cohort","status"))
names(cK_L) = c("sample","cohort","status","kingdom","RelativeAbundance")

p = ggboxplot(cK_L,"status","RelativeAbundance") + stat_compare_means(method.args=list(alternative="less"),comparisons=list(c("domestic","wild")),size=3)
facet(p,facet.by="kingdom",scale="free")
