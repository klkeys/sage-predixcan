##------ [/media/BurchardRaid01/LabShare/Home/kkeys/gala_sage/rnaseq/glmnet] Tue Nov 28 10:16:13 2017 ------##
x = fread("sage_glmnet_results.txt")
names(x)
colnames(x)
x[!is.na(x$pval),]
dim(x[!is.na(x$pval),])
dim(x[!is.na(x$pval),]) -> x.nona
head(x.nona)
x[!is.na(x$pval),] -> x.nona
head(x.nona)
x.nona[x.nona$pval < 0.05,]
dim(x.nona[x.nona$pval < 0.05,])
dim(x.nona[x.nona$R2 > 0.6,])
dim(x.nona[x.nona$R2 > 0.5,]
x.nona[x.nona$pval < 0.05,]$R2
mean(x.nona[x.nona$pval < 0.05,]$R2)
library(glmnet)
?glmnet
ls
?glmnet
?cv.glmnet
x.nona
x.nona.p05 = x.nona[x.nona$pval < 0.05,]
dim(x.nona.p05)
history()
h = fread("../gcta/sage_38_h2.txt")
head(h)
?intersect
?merge
head(h)
head(x.nona.p05)
merge(x.nona.p05, h, by.x = "gene", by.y = "Gene")
dim(h)
dim(x.nona.p05)
dim(merge(x.nona.p05, h, by.x = "gene", by.y = "Gene"))
x.h2 = merge(x.nona.p05, h, by.x = "gene", by.y = "Gene")
mean(x.h2$h2)
summary(lm(R2 ~ h2, data = x.h2))
sum(x.h2$R2 > x.h2$h2)
sum(x.h2$R2 <= x.h2$h2)
head(x.h2)
savehistory("sage_rnaseq_elasticnet_weight_analysis_28NOV2017.R")
