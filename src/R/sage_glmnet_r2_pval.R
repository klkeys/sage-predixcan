##------ [/media/BurchardRaid01/LabShare/Home/kkeys/gala_sage/rnaseq/glmnet] Thu Nov 30 10:45:47 2017 ------##
x = fread("sage_gtex_r2_defaultweights.txt")
head(x)
x[!is.na(x$V3),]
x[!is.na(x$V3) & !is.nan(x$V4),]
dim(x[!is.na(x$V3) & !is.nan(x$V4),])
dim(x[!is.na(x$V3),])
colnames(x) = c("gene", "nsnps", "R2", "R2.pval")
x.enet = fread("./elasticnet/sage_elasticnet_results.txt")
x.lasso = fread("./elasticnet/sage_lasso_results.txt")
x.lasso = fread("./lasso/sage_lasso_results.txt")
head(x.lasso)
colnames(x.lasso)[7:8] = c("R2.lasso", "R2.pval.lasso")
head(x.lasso)
colnames(x.enet)[7:8] = c("R2.enet", "R2.pval.enet")
head(x.enet)
cbind(x.lasso[,7:8], x.enet[,7:8]
)
cbind(x[,7:8],x.lasso[,7:8], x.enet[,7:8]
)
cbind(x[,3:4],x.lasso[,7:8], x.enet[,7:8]
)
genes.in.common = intersect(x$V1, x.enet$gene, x.lasso$gene)
genes.in.common = intersect(intersect(x$V1, x.enet$gene), x.lasso$gene)
genes.in.common
x$v1
x$V1
head(x)
genes.in.common = intersect(intersect(x$gene, x.enet$gene), x.lasso$gene)
head(x)
genes.in.common
x.common = x[x$gene %in% genes.in.common,]
x.lasso.common = x[x.lasso$gene %in% genes.in.common,]
dim(x.lasso)
length(genes.in.common)
x.lasso.common = x.lasso[x.lasso$gene %in% genes.in.common,]
x.enet.common = x.enet[x.enet$gene %in% genes.in.common,]
cbind(x.common[,3:4], x.lasso.common[,7:8], x.enet.common[,7:8])
dim(x.lasso.common)
dim(x.enet.common)
dim(x.common)
length(genes.in.common)
sum(x.gene %in% genes.in.common)
sum(x$gene %in% genes.in.common)
duplicated(x)
x[1:20,]
x[duplicated(x),]
?unique
dim(unique(x))
x = unique(x)
x.common = x[x$gene %in% genes.in.common,]
cbind(x.common[,3:4], x.lasso.common[,7:8], x.enet.common[,7:8])
x.r2.all = cbind(x.common[,3:4], x.lasso.common[,7:8], x.enet.common[,7:8])
sum(x.r2.all$R2 < x.r2.all$R2.lasso)
sum(x.r2.all$R2 < x.r2.all$R2.lasso, na.rm = T)
sum(x.r2.all$R2 < x.r2.all$R2.enet, na.rm = T)
sum(x.r2.all$R2 <= x.r2.all$R2.enet, na.rm = T)
sum(x.r2.all$R2 <= x.r2.all$R2.lasso, na.rm = T)
sum(x.r2.all$R2 >= x.r2.all$R2.lasso, na.rm = T)
sum(x.r2.all$R2 >= x.r2.all$R2.enet, na.rm = T)
sum(x.r2.all$R2.lasso >= x.r2.all$R2.enet, na.rm = T)
sum(x.r2.all$R2.lasso <= x.r2.all$R2.enet, na.rm = T)
fwrite(x=x.r2.all, file="sage_glmnet_r2_pval_results.txt", row.names=F, col.names=T, quote=F)
savehistory("sage_glmnet_r2_pval.R")
