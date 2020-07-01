#
sm = read.table("sample.metadata",head=T);sm$cohort=factor(sm$cohort,levels=c("FinF","SpaR","SpaS","FinR","ArgCTL","SpaW","Jp","Sib","ArgW"))
aC = data.frame(cohort=sm$cohort,status=sm$status)
Var1 <- c("#86CEEB","#B03060"); names(Var1) = c("wild","domestic")
nn_colors <- list(Status=Var1)	
A = as.matrix(read.table("HorseV3.phenotype.Relative.profile"))
mat_breaks <- quantile_breaks(A, n = 11)
rownames(aC) = sm$sample
pheatmap(A, cluster_cols=FALSE, color= colorRampPalette(c("white","#87CEEB", "#B03060"))(1000),annotation_col= aC,cutree_rows=2)
