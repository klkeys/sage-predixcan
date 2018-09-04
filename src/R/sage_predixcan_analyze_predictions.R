#!/usr/bin/env Rscript --vanilla
# =======================================================================================
# copyright Asthma Collaboratory (2018)
# coded by Kevin L. Keys
#
# This script processes and analyzes PrediXcan predictions in SAGE data 
# =======================================================================================

# =======================================================================================
# load libraries 
# =======================================================================================
library(data.table)
library(purrr)
library(broom)
library(ggplot2)

# =======================================================================================
# subroutines 
# =======================================================================================

lmtest = function(x, lm.formula) {
    x.nona = na.omit(x)
    n = dim(x.nona)[1]
    if ( n < 1 ) {
        return(data.table(Gene = x$Gene, R2 = NA, N = 0))
    } 
    my.lm = lm(formula(lm.formula), data = x.nona, na.action = na.omit)
    return(data.table(Gene = x$Gene, R2 = summary(my.lm)$r.squared, N = n))
}

# subroutine for extracting both Spearman rho, p-value from correlation test
cortest = function(x, cor.formula) {
    x.nona = na.omit(x)
    n = dim(x.nona)[1]
    if ( n < 1 ) {
        return(data.table(Gene = x$Gene, Correlation = NA, Corr.p.value = NA))
    } 
    my.cortest = cor.test(formula(cor.formula), x.nona, method = "spearman", na.action = na.omit)
    return(data.table(Gene = x.nona$Gene, Correlation = my.cortest$estimate, Corr.p.value = my.cortest$p.value))
}

# =======================================================================================
# file and directory paths 
# =======================================================================================
predixcan.dir  = file.path(Sys.getenv("HOME"), "gala_sage", "rnaseq", "predixcan")
data.dir       = file.path(Sys.getenv("HOME"), "gala_sage", "rnaseq", "data")

dgn.path       = file.path(predixcan.dir, "DGN",  "DGN_predixcan_prediction_ALLCHR_predicted_expression_melted.txt")
gtex6.path     = file.path(predixcan.dir, "GTEx", "GTEx_v6p_predixcan_prediction_ALLCHR_predicted_expression_melted.txt")
gtex7.path     = file.path(predixcan.dir, "GTEx", "GTEx_v7_predixcan_prediction_ALLCHR_predicted_expression_melted.txt")
mesa.afa.path  = file.path(predixcan.dir, "MESA", "MESA_AFA_predixcan_prediction_ALLCHR_predicted_expression_melted.txt")
mesa.afhi.path = file.path(predixcan.dir, "MESA", "MESA_AFHI_predixcan_prediction_ALLCHR_predicted_expression_melted.txt")
mesa.cau.path  = file.path(predixcan.dir, "MESA", "MESA_CAU_predixcan_prediction_ALLCHR_predicted_expression_melted.txt")
mesa.all.path  = file.path(predixcan.dir, "MESA", "MESA_ALL_predixcan_prediction_ALLCHR_predicted_expression_melted.txt")

sage.rna.path  = file.path(data.dir, "sage_39_wgs_for_rnaseq_expression_sorted_headered_nocomment.bed")

# =======================================================================================
# load data
# =======================================================================================
dgn       = fread(dgn.path)
gtex6     = fread(gtex6.path)
gtex7     = fread(gtex7.path)
mesa.afa  = fread(mesa.afa.path)
mesa.afhi = fread(mesa.afhi.path)
mesa.cau  = fread(mesa.cau.path)
mesa.all  = fread(mesa.all.path)

sage = fread(sage.rna.path, header = TRUE)

# =======================================================================================
# merge data into single data frame  
# =======================================================================================

# must rename columns of each data.table, particularly the ones with predicted expression values
# this facilitates merging them later
repos = c("DGN", "GTEx_v6p", "GTEx_v7", "MESA_AFA", "MESA_AFHI", "MESA_CAU", "MESA_ALL")
repo.results = list(dgn, gtex6, gtex7, mesa.afa, mesa.afhi, mesa.cau, mesa.all)
for (i in 1:length(repos)) {
    colnames(repo.results[[i]]) = c("SubjectID", "Gene", paste0("Predicted_Expr_", repos[i]))
}

# perform a full join of all prediction results
predixcan.all = repo.results %>% reduce(full_join, by = c("SubjectID","Gene")) %>% as.data.table

# melt measurements prior to merge
sage.melt = melt(sage[,-c(1:3)], id.vars = "Gene", variable.name = "SubjectID", value.name = "Measured_Expr")

# merge predictions with measurements
sage.predixcan.all = merge(sage.melt, predixcan.all, by = c("Gene", "SubjectID"))


# =======================================================================================
# analyze predictions 
# =======================================================================================

# seed a list to save results
### NOT DONE!!!! ####
### preallocate the list with correct names
repo.results.predvmeas = list() 


# compute genewise regressions
for (i in 1:length(repos)) {
    my.repo = repos[i]
    my.lm.formula = paste0("Predicted_Expr_", my.repo, " ~ Measured_Expr")
    my.corr.formula = paste0("~ Predicted_Expr_", my.repo, " + Measured_Expr")

    r2s = sage.predixcan.all %>%
        group_by(Gene) %>%
        #na.omit %>%
        #do(tidy(summary(lm(my.formula, data = .))$r.squared)) %>%
        do(lmtest(., my.lm.formula)) %>%
        #select(Gene, x) %>%
        as.data.table %>%
        #na.omit %>%
        unique
    colnames(r2s) = c("Gene", paste0("R2_", my.repo), paste0("Num_Pred_", my.repo))

    corrs = sage.predixcan.all %>%
        group_by(Gene) %>%
       # na.omit %>%
        do(cortest(., my.corr.formula)) %>% # use previous subroutine to perform correlation test and extract three columns (gene, correlation, p-value)
        as.data.table %>% # cast as data.table to purge duplicate rows
        #na.omit %>%
        unique # need this because inelegant subroutine prints repeated rows, 1 per sample instead of 1 per gene group
    colnames(corrs) = c("Gene", paste0("Corr_", my.repo), paste0("Corr_pval_", my.repo))
    
    my.results = merge(r2s, corrs, by = c("Gene"), all = TRUE)
    repo.results.predvmeas[[my.repo]] = my.results 
    
}

predixcan.all = sage.melt = corrs = r2s = FALSE
gc()
sage.predixcan.all.results = repo.results.predvmeas %>% reduce(full_join, by = "Gene") %>% as.data.table
repo.results.predvmeas = FALSE
gc()

# =======================================================================================
# compile testing metrics
# =======================================================================================

# compare imputation performance from all four prediction weights
# must first load GTEx v7 testing R2s from PredictDB
predixcan.gtex7.r2s.path =  file.path(Sys.getenv("HOME"), "gala_sage", "rnaseq", "predixcan", "gtex7.testR2.txt")
predixcan.gtex7.r2s = fread(predixcan.gtex7.r2s.path) 

# =======================================================================================
# plot results
# =======================================================================================

# set color palette and boxplot labels
my.colors = c("blue", "orange", "red", "gold")

# gplot2 colorblind-friendly palette with grey:
cbPalette = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
my.boxplot.labels = c("GTEx v6p", "GTEx v7", "DGN", "MESA_AFA", "MESA_AFHI", "MESA_CAU", "MESA_ALL", "GTEx7 Test")


# extract coefficients of determination (R2) from genewise summary results
# make sure to rename the columns!
all.r2.commongenes = sage.predixcan.all.results %>%
    merge(., predixcan.gtex7.r2s, by = "Gene", all = T) %>%
    select(Gene, R2_GTEx_v6p, R2_GTEx_v7, R2_DGN, R2_MESA_AFA, R2_MESA_AFHI, R2_MESA_CAU, R2_MESA_ALL, GTEx7_test_R2_avg) %>%
    as.data.table %>%
    na.omit %>%
    melt
colnames(all.r2.commongenes) = c("Gene", "Prediction.Weights", "R2")
setorderv(all.r2.commongenes, "Prediction.Weights")

# will produce three kinds of summary plots:
# (1): boxplot
# (2): violin plot
# (3): histogram + density plot 
my.boxplot.r2 = ggplot(all.r2.commongenes, aes(x = Prediction.Weights, y = R2, fill = Prediction.Weights)) +
    geom_boxplot() +
    xlab("Prediction Weight Set") +
    ylab(expression(R^{2})) +
    ggtitle(expression(paste("Distribution of ", R^{2}, " across different prediction weight sets"))) +
    scale_fill_manual(values = cbPalette) +
    scale_x_discrete(labels = my.boxplot.labels) +
    theme(legend.position = "none")

my.violinplot.r2 = ggplot(all.r2.commongenes, aes(x = Prediction.Weights, y = R2)) +
    geom_violin() +
    xlab("Prediction Weight Set") +
    ylab(expression(R^{2})) +
    ggtitle("Distribution of R2 across different prediction weight sets")

my.hist.r2 = ggplot(all.r2.commongenes, aes(x = R2)) +
    geom_histogram(aes(y = ..density..)) +
    geom_density() +
    xlab("Prediction Weight Set") +
    ylab(expression(R^{2})) +
    ggtitle("Distribution of R2 across different prediction weight sets") +
    facet_wrap(~ Prediction.Weights)
 
# save plots to file
ggsave(plot = my.boxplot.r2, filename = "all.r2.commongenes.boxplot.png", dpi = 300, type = "cairo")
ggsave(plot = my.violinplot.r2, filename = "all.r2.commongenes.violinplot.png", dpi = 300, type = "cairo")
ggsave(plot = my.hist.r2, filename = "all.r2.commongenes.histogram.png", dpi = 300, type = "cairo")

# do same, but for correlations instead of R2s
all.rho.commongenes = sage.predixcan.all.results %>%
    merge(., predixcan.gtex7.r2s, by = "Gene", all = T) %>%
    select(Gene, Corr_GTEx_v6p, Corr_GTEx_v7, Corr_DGN, Corr_MESA_AFA, Corr_MESA_AFHI, Corr_MESA_CAU, Corr_MESA_ALL, GTEx7_test_rho_avg) %>%
    as.data.table %>%
    na.omit %>%
    melt
colnames(all.rho.commongenes) = c("Gene", "Prediction.Weights", "Correlation")
setorderv(all.rho.commongenes, "Prediction.Weights")

my.hist.rho = ggplot(all.rho.commongenes, aes(x = Correlation)) +
    geom_histogram(aes(y = ..density..)) +
    geom_density() +
    xlab("Prediction Weight Set") +
    ylab("Spearman rho") +
    ggtitle("Distribution of R2 across different prediction weight sets") +
    facet_wrap(~ Prediction.Weights)
my.boxplot.rho = ggplot(all.rho.commongenes, aes(x = Prediction.Weights, y = Correlation, fill = Prediction.Weights)) +
    geom_boxplot() +
    xlab("Prediction Weight Set") +
    ylab(expression(paste("Correlations (Spearman ", rho, ")"))) +
    ggtitle("Distribution of correlations across different prediction weight sets") +
    scale_fill_manual(values = cbPalette) +
    theme(legend.position = "none") +
    scale_x_discrete(labels = my.boxplot.labels)
my.violinplot.rho = ggplot(all.rho.commongenes, aes(x = Prediction.Weights, y = Correlation)) +
    geom_violin() +
    xlab("Prediction Weight Set") +
    ylab("Correlations") +
    ggtitle("Distribution of correlations across different prediction weight sets")

ggsave(plot = my.hist.rho, filename = "all.rho.commongenes.histogram.png")
ggsave(plot = my.boxplot.rho, filename = "all.rho.commongenes.boxplot.png")
ggsave(plot = my.violinplot.rho, filename = "all.rho.commongenes.violinplot.png")
