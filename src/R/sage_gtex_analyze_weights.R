#!/usr/bin/env Rscript --vanilla
#
### parse command line arguments
##args           = commandArgs(trailingOnly = TRUE)
##glmnet.results = args[1]
##h2.file        = args[2]
##glmnet.method  = args[3]
#
## load results 
##x = fread("./elasticnet/sage_elasticnet_results.txt")
##x = fread(glmnet.results)
#x.lasso = fread("sage_lasso_results.txt")
#x.enet  = fread("sage_elasticnet_results.txt")
#x.ridge = fread("sage_ridge_results.txt"
#
## purge missing values
##x.nona = x[!is.na(x$pval),]
#
## get genes with marginally significant p-values
##x.nona.p05 = x.nona[x.nona$pval < 0.05,]
#
## load heritability results
#h = fread("../gcta/sage_38_h2.txt")
##h = fread(h2.file)
#
## merge with parsed results
#x.h2 = merge(x.nona.p05, h, by.x = "gene", by.y = "Gene")


### subroutines
 # make function to compute cor.test over each individual column
 cor.test.twocol = function(x, i, stride = col.gap, method = "spearman"){
     return(cor.test(x[,i], x[,i + stride], method=method))
 }
 
 # function to compute coefficient of determination (R2) over columns
 #r2 = function(measured, predicted){ return(1 - (sum((measured - predicted)^2)/sum((measured - mean(measured, na.rm=T))^2)))}
 
 # compute Spearman correlation test and R2 for all columns of combined matrix 
 my.cor.test = function(x, stride = col.gap, method = "spearman"){
     p = dim(x)[2]
     q = p - stride
     results = matrix(-Inf, q, 3)
     for (i in 1:q){
         my.test = cor.test.twocol(x, i, stride=stride, method=method)
         my.r2 = summary(lm(x[,i + stride] ~ x[,i]))$r.squared
         results[i,] = c(my.r2, my.test$estimate, my.test$p.value)
     }   
     genenames = strtrim(colnames(x)[1:stride], 15) 
     df = data.frame("Gene" = genenames, "R2" = results[,1], "Spearman Rho" = results[,2], "P-value" = results[,3])
     return(df)
 } 
### end subroutines

# load expression data
exprs = fread("../data/sage_39_wgs_for_rnaseq_expression_sorted_headered.bed")

# put header on data.table exprs
my.header = "Chromosome_Name Start_Position End_Position Gene NWD159235 NWD262884 NWD424160 NWD331495 NWD511506 NWD742593 NWD154420 NWD101012 NWD843708 NWD990899 NWD373455 NWD595401 NWD299516 NWD713489 NWD468153 NWD821729 NWD167864 NWD127715 NWD952849 NWD975237 NWD212794 NWD841275 NWD359298 NWD251674 NWD620717 NWD537216 NWD863759 NWD671985 NWD923487 NWD621320 NWD345359 NWD769653 NWD152277 NWD765179 NWD546278 NWD833668 NWD437999 NWD730618 NWD554260"
strsplit(my.header, " ")[[1]]
colnames(exprs) = strsplit(my.header, " ")[[1]]

# shuffle rows to alphanumeric order of **genes**
setorder(exprs, Gene)

# transpose the matrix
# first 4 cols are (useless) BED info
# must manually put row/col names on transpose
texpr = as.matrix(t(exprs[,-c(1:4)]))
colnames(texpr) = exprs$Gene
row.names(texpr) = colnames(exprs)[-c(1:4)]

# shuffle rows of texpr to alphanumeric order of **samples**
# this should ensure ascending order L-->R for columns
# and same top --> bottom for rows 
texpr = texpr[order(row.names(texpr)),]

# now load prediction data
preds = fread("sage_elasticnet_results.txt")

# transpose and rename
tpred = as.matrix(t(preds[,-1]))
colnames(tpred) = preds$Gene
row.names(tpred) = sort(colnames(exprs)[-c(1:4)])

# there may be more measurements than predictions
# find genes not in common
# must purge them from analyses 
gene.diff = setdiff(colnames(texpr), colnames(tpred)) # "ENSG00000264063" "ENSG00000264462"
gene.diff.idx = which(colnames(texpr) %in% gene.diff) # 6848 6863
texpr.common = texpr[, -gene.diff.idx]

# combine data frames
sage.rnapred.common = cbind(texpr.common, tpred)
col.gap = dim(texpr.common)[2]
 
# compute R2
corrs = my.cor.test(sage.rnapred.common, col.gap, method = "spearman")

# get a column of Benjamini-Hochberg q-values for FDR adjustment
# do same for Bonferroni-corrected p-values
bh.alpha = 0.05
bonf.level = 0.05 / col.gap # ~ 6.78e-6
corrs$BH.qvalue = p.adjust(corrs$P.value, method = "BH")
corrs.sig = corrs[corrs$P.value < 0.05,]
corrs.bh  = corrs[corrs$BH.qvalue < bh.alpha,]
corrs.bonf = corrs[corrs$P.value < bonf.level,]

# get R2s and corresponding p-values
# now regress each predicted expression onto the measured expression
reg.values = matrix(-1, col.gap, 3)
rownames(reg.values) = colnames(texpr.common)
colnames(reg.values) = c("R2", "p.value", "matching.signs")
#row.names(reg.values) = colnames(sage.rnapred.common)[1:col.gap]
for (i in 1:col.gap){
    x.pred = sage.rnapred.common[,i + col.gap]
    x.expr = sage.rnapred.common[,i]
    my.lm = summary(lm(x.expr ~ x.pred ))
    matching.signs = sum(sign(x.pred) == sign(x.expr))
    #gene.name = colnames(sage.rnapred.common)[i]
    rsq = my.lm$r.squared
    pval = my.lm$coef[2,4] 
    #row.names(reg.values)[i] = gene.name
    new.row = c(rsq, pval, matching.signs) 
#print(new.row)
    reg.values[i,] = new.row 
}

corrs$p.value.reg = reg.values[,2]
corrs$matching.signs = reg.values[,3]

fwrite(x=corrs, file="sage_predixcan_correlations.txt", row.names = F, col.names = T, quote = F)
fwrite(x=sage.rnapred.common, file="sage_rnapred_common.txt", row.names = F, col.names = T, quote = F)
