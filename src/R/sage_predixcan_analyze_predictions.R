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
library(dplyr)

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

sage.rna.path  = file.path(data.dir, "sage_39_wgs_for_rnaseq_expression_melted.txt")

rdata.path     = file.path("sage_predixcan_allresults_allplots.Rdata")

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
#sage.melt = melt(sage[,-c(1:3)], id.vars = "Gene", variable.name = "SubjectID", value.name = "Measured_Expr")
sage.melt = sage
names(sage.melt) = c("Gene", "SubjectID", "Measured_Expr")

# merge predictions with measurements
sage.predixcan.all = merge(sage.melt, predixcan.all, by = c("Gene", "SubjectID"), all = TRUE)


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
    #colnames(r2s) = c("Gene", paste0("R2_", my.repo), paste0("Num_Pred_", my.repo))
    colnames(r2s) = c("Gene", "R2", "Num_Pred")
    r2s$Repo = my.repo

    corrs = sage.predixcan.all %>%
        group_by(Gene) %>%
       # na.omit %>%
        do(cortest(., my.corr.formula)) %>% # use previous subroutine to perform correlation test and extract three columns (gene, correlation, p-value)
        as.data.table %>% # cast as data.table to purge duplicate rows
        #na.omit %>%
        unique # need this because inelegant subroutine prints repeated rows, 1 per sample instead of 1 per gene group
    #colnames(corrs) = c("Gene", paste0("Corr_", my.repo), paste0("Corr_pval_", my.repo))
    colnames(corrs) = c("Gene", "Corr", "Corr_pval")
    corrs$Repo = my.repo
    
    my.results = merge(r2s, corrs, by = c("Gene", "Repo"), all = TRUE)
    repo.results.predvmeas[[my.repo]] = my.results 
    
}

predixcan.all = sage.melt = corrs = r2s = FALSE
gc()
#sage.predixcan.all.results = repo.results.predvmeas %>% reduce(full_join, by = "Gene") %>% as.data.table
sage.predixcan.all.results = rbindlist(repo.results.predvmeas) 
repo.results.predvmeas = FALSE
gc()


# =======================================================================================
# plot results
# =======================================================================================

# set color palette and boxplot labels
my.colors = c("blue", "orange", "red", "gold")
#my.colors = c("blue", "red", "gold")

# gplot2 colorblind-friendly palette with grey:
cbPalette = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999", "#E69F00", "#56B4E9", "#009E73") 
my.boxplot.labels = c("GTEx v6p", "GTEx v7", "DGN", "MESA_AFA", "MESA_AFHI", "MESA_CAU", "MESA_ALL", "Test: GTEx v7", "Test: MESA AFA", "Test: MESA AFHI", "Test: MESA CAU", "Test: MESA ALL")

# compare imputation performance from all four prediction weights
# must first load GTEx v7 testing R2s from PredictDB
predixcan.gtex7.metrics.path     = file.path(Sys.getenv("HOME"), "gala_sage", "rnaseq", "predixcan", "gtex7.test.metrics.txt")
predixcan.mesa.afa.metrics.path  = file.path(Sys.getenv("HOME"), "gala_sage", "rnaseq", "predixcan", "mesa.AFA.test.metrics.txt")
predixcan.mesa.afhi.metrics.path = file.path(Sys.getenv("HOME"), "gala_sage", "rnaseq", "predixcan", "mesa.AFHI.test.metrics.txt")
predixcan.mesa.cau.metrics.path  = file.path(Sys.getenv("HOME"), "gala_sage", "rnaseq", "predixcan", "mesa.CAU.test.metrics.txt")
predixcan.mesa.all.metrics.path  = file.path(Sys.getenv("HOME"), "gala_sage", "rnaseq", "predixcan", "mesa.ALL.test.metrics.txt")

predixcan.gtex7.metrics     = fread(predixcan.gtex7.metrics.path) 
predixcan.mesa.afa.metrics  = fread(predixcan.mesa.afa.metrics.path) 
predixcan.mesa.afhi.metrics = fread(predixcan.mesa.afhi.metrics.path) 
predixcan.mesa.cau.metrics  = fread(predixcan.mesa.cau.metrics.path) 
predixcan.mesa.all.metrics  = fread(predixcan.mesa.all.metrics.path) 

# must rename columns of each data.table, particularly the ones with predicted expression values
# this facilitates merging them later
repos.test = c("Test: GTEx_v7", "Test: MESA_AFA", "Test: MESA_AFHI", "Test: MESA_CAU", "Test: MESA_ALL")
repo.test = list(predixcan.gtex7.metrics, predixcan.mesa.afa.metrics, predixcan.mesa.afhi.metrics, predixcan.mesa.cau.metrics, predixcan.mesa.all.metrics)
for (i in 1:length(repos.test)) {
    #colnames(repo.test[[i]]) = c("Gene", "HUGO", paste0(repos.test[i], "_test_R2"), paste0(repos.test[i], "_test_Corr"))
    colnames(repo.test[[i]]) = c("Gene", "Repo", "R2", "Corr")
    #repo.test[[i]] = repo.test[[i]][,-2]
    repo.test[[i]]$Repo = repos.test[i]
    repo.test[[i]]$Gene = strtrim(repo.test[[i]]$Gene, 15)  ## trim .X, the transcript number to the ENSG ID
}

# perform a full join of all prediction results
#predixcan.all.test = repo.test %>% reduce(full_join, by = c("Gene")) %>% as.data.table
predixcan.all.test = rbindlist(repo.test)
repo.test = FALSE; gc()  ## recover memory

# extract coefficients of determination (R2) from genewise summary results
# make sure to rename the columns!
#all.r2.commongenes = sage.predixcan.all.results %>%
    #merge(., predixcan.all.test, by = "Gene", all = T) %>%
    #select(Gene, R2_GTEx_v6p, R2_GTEx_v7, R2_DGN, R2_MESA_AFA, R2_MESA_AFHI, R2_MESA_CAU, R2_MESA_ALL, GTEx_v7_test_R2, MESA_AFA_test_R2, MESA_AFHI_test_R2, MESA_CAU_test_R2, MESA_ALL_test_R2) %>%
all.r2 = rbindlist(list(sage.predixcan.all.results, predixcan.all.test), fill = TRUE) %>%
    select(Gene, Repo, R2) %>%
    as.data.table %>%
    na.omit #%>%
    #melt
colnames(all.r2) = c("Gene", "Prediction.Weights", "R2")
setorderv(all.r2, "Prediction.Weights")
ngenes.r2 = all.r2 %>% select(Gene) %>% unlist %>% sort %>% unique %>% length

# will produce three kinds of summary plots:
# (1): boxplot
# (2): violin plot
# (3): histogram + density plot 
my.boxplot.r2 = ggplot(all.r2, aes(x = Prediction.Weights, y = R2, fill = Prediction.Weights)) +
    geom_boxplot() +
    xlab("Prediction Weight Set") +
    ylab(expression(R^{2})) +
    ggtitle(bquote("Distribution of" ~ R^2 ~ "across different prediction weight sets over" ~ .(ngenes.r2) ~ "genes")) +
    scale_fill_manual(values = cbPalette) +
    scale_x_discrete(labels = my.boxplot.labels) +
    ylim(-0.1, 1) +
    theme(legend.position = "none")

my.violinplot.r2 = ggplot(all.r2, aes(x = Prediction.Weights, y = R2)) +
    geom_violin() +
    xlab("Prediction Weight Set") +
    ylab(expression(R^{2})) +
    ggtitle(bquote("Distribution of" ~ R^2 ~ "across different prediction weight sets over" ~ .(ngenes.r2) ~ "genes"))

my.hist.r2 = ggplot(all.r2, aes(x = R2)) +
    geom_histogram(aes(y = ..density..)) +
    geom_density() +
    xlab("Prediction Weight Set") +
    ylab(expression(R^{2})) +
    ggtitle(bquote("Distribution of" ~ R^2 ~ "across different prediction weight sets over" ~ ngenes.r2 ~ "genes")) +
    facet_wrap(~ Prediction.Weights)
 
# save plots to file
ggsave(plot = my.boxplot.r2, filename = "all.r2.boxplot.png", dpi = 300, type = "cairo", width = 15, height = 5, units = "in")
ggsave(plot = my.violinplot.r2, filename = "all.r2.violinplot.png", dpi = 300, type = "cairo")
ggsave(plot = my.hist.r2, filename = "all.r2.histogram.png", dpi = 300, type = "cairo")


# do for only common genes
r2.commongenes = all.r2 %>%
    group_by(Gene) %>%
    tally %>%
    dplyr::filter(n > 11) %>%
    select(Gene) %>%
    as.data.table
all.r2.commongenes = all.r2 %>% dplyr::filter(Gene %in% r2.commongenes$Gene)
ncommongenes.r2 = r2.commongenes %>% select(Gene) %>% unlist %>% sort %>% unique %>% length 

my.boxplot.r2.commongenes = ggplot(all.r2.commongenes, aes(x = Prediction.Weights, y = R2, fill = Prediction.Weights)) +
    geom_boxplot() +
    xlab("Prediction Weight Set") +
    ylab(expression(R^{2})) +
    ggtitle(bquote("Distribution of " ~ R^2 ~ "across different prediction weight sets over" ~ .(ncommongenes.r2) ~ "common genes")) +
    scale_fill_manual(values = cbPalette) +
    scale_x_discrete(labels = my.boxplot.labels) +
    ylim(-0.1, 1) +
    theme(legend.position = "none")

my.violinplot.r2.commongenes = ggplot(all.r2.commongenes, aes(x = Prediction.Weights, y = R2)) +
    geom_violin() +
    xlab("Prediction Weight Set") +
    ylab(expression(R^{2})) +
    ggtitle(bquote("Distribution of" ~ R^2 ~ "across different prediction weight sets over" ~ .(ncommongenes.r2) ~ "common genes"))

my.hist.r2.commongenes = ggplot(all.r2.commongenes, aes(x = R2)) +
    geom_histogram(aes(y = ..density..)) +
    geom_density() +
    xlab("Prediction Weight Set") +
    ylab(expression(R^{2})) +
    ggtitle(bquote("Distribution of" ~ R^2 ~ "across different prediction weight sets over" ~ .(ncommongenes.r2) ~ "common genes")) +
    facet_wrap(~ Prediction.Weights)

# save plots to file
ggsave(plot = my.boxplot.r2.commongenes, filename = "all.r2.commongenes.boxplot.png", dpi = 300, type = "cairo", width = 15, height = 5, units = "in")
ggsave(plot = my.violinplot.r2.commongenes, filename = "all.r2.commongenes.violinplot.png", dpi = 300, type = "cairo")
ggsave(plot = my.hist.r2.commongenes, filename = "all.r2.commongenes.histogram.png", dpi = 300, type = "cairo")


# do same, but for correlations instead of R2s
#all.rho.commongenes = sage.predixcan.all.results %>%
#    merge(., predixcan.all.test, by = "Gene", all = T) %>%
    #select(Gene, Corr_GTEx_v6p, Corr_GTEx_v7, Corr_DGN, Corr_MESA_AFA, Corr_MESA_AFHI, Corr_MESA_CAU, Corr_MESA_ALL, GTEx_v7_test_Corr, MESA_AFA_test_Corr, MESA_AFHI_test_Corr, MESA_CAU_test_Corr, MESA_ALL_test_Corr) %>%
all.rho = rbindlist(list(sage.predixcan.all.results, predixcan.all.test), fill = TRUE) %>%
    select(Gene, Repo, Corr) %>% 
    as.data.table %>%
    na.omit #%>%
    #melt
colnames(all.rho) = c("Gene", "Prediction.Weights", "Correlation")
setorderv(all.rho, "Prediction.Weights")
ngenes.rho = all.rho  %>% select(Gene) %>% unlist %>% sort %>% unique %>% length 

my.hist.rho = ggplot(all.rho, aes(x = Correlation)) +
    geom_histogram(aes(y = ..density..)) +
    geom_density() +
    xlab("Prediction Weight Set") +
    ylab("Spearman rho") +
    ggtitle(paste0("Distribution of correlations across different prediction weight sets for ", ngenes.rho, " genes")) +
    facet_wrap(~ Prediction.Weights)
my.boxplot.rho = ggplot(all.rho, aes(x = Prediction.Weights, y = Correlation, fill = Prediction.Weights)) +
    geom_boxplot() +
    xlab("Prediction Weight Set") +
    ylab(expression(paste("Correlations (Spearman ", rho, ")"))) +
    ggtitle(paste0("Distribution of correlations across different prediction weight sets for ", ngenes.rho, " genes")) +
    scale_fill_manual(values = cbPalette) +
    theme(legend.position = "none") +
    scale_x_discrete(labels = my.boxplot.labels)
my.violinplot.rho = ggplot(all.rho, aes(x = Prediction.Weights, y = Correlation)) +
    geom_violin() +
    xlab("Prediction Weight Set") +
    ylab("Correlations") +
    ggtitle(paste0("Distribution of correlations across different prediction weight sets for ", ngenes.rho, " genes"))

ggsave(plot = my.boxplot.rho, filename = "all.rho.boxplot.png", dpi = 300, type = "cairo", width = 15, height = 5, units = "in")
ggsave(plot = my.hist.rho, filename = "all.rho.histogram.png")
ggsave(plot = my.violinplot.rho, filename = "all.rho.violinplot.png")

rho.commongenes = all.rho %>%
    group_by(Gene) %>%
    tally %>%
    dplyr::filter(n > 11) %>%
    select(Gene) %>%
    as.data.table
all.rho.commongenes = all.rho %>% dplyr::filter(Gene %in% rho.commongenes$Gene)
ncommongenes.rho = all.rho.commongenes %>% select(Gene) %>% unlist %>% sort %>% unique %>% length 

my.hist.rho.commongenes = ggplot(all.rho.commongenes, aes(x = Correlation)) +
    geom_histogram(aes(y = ..density..)) +
    geom_density() +
    xlab("Prediction Weight Set") +
    ylab("Spearman rho") +
    ggtitle(paste0("Distribution of correlations across different prediction weight sets for ", ncommongenes.rho, " genes in common")) +
    facet_wrap(~ Prediction.Weights)
my.boxplot.rho.commongenes = ggplot(all.rho.commongenes, aes(x = Prediction.Weights, y = Correlation, fill = Prediction.Weights)) +
    geom_boxplot() +
    xlab("Prediction Weight Set") +
    ylab(expression(paste("Correlations (Spearman ", rho, ")"))) +
    ggtitle(paste0("Distribution of correlations across different prediction weight sets for ", ncommongenes.rho, " genes in common")) +
    scale_fill_manual(values = cbPalette) +
    theme(legend.position = "none") +
    scale_x_discrete(labels = my.boxplot.labels)
my.violinplot.rho.commongenes = ggplot(all.rho.commongenes, aes(x = Prediction.Weights, y = Correlation)) +
    geom_violin() +
    xlab("Prediction Weight Set") +
    ylab("Correlations") +
    ggtitle(paste0("Distribution of correlations across different prediction weight sets for ", ncommongenes.rho, " genes in common"))

ggsave(plot = my.boxplot.rho.commongenes, filename = "all.rho.commongenes.boxplot.png", dpi = 300, type = "cairo", width = 15, height = 5, units = "in")
ggsave(plot = my.hist.rho.commongenes, filename = "all.rho.commongenes.histogram.png", dpi = 300, type = "cairo")
ggsave(plot = my.violinplot.rho.commongenes, filename = "all.rho.commongenes.violinplot.png", dpi = 300, type = "cairo")

# compute some statistics
r2.summaries  = all.r2 %>% group_by(Prediction.Weights) %>% summarize(r2 = mean(R2, na.rm = T)) %>% as.data.table
r2.commongenes.summaries  = all.r2.commongenes %>% group_by(Prediction.Weights) %>% summarize(r2 = mean(R2, na.rm = T)) %>% as.data.table
rho.summaries = all.rho %>% group_by(Prediction.Weights) %>% summarize(rho = mean(Correlation, na.rm = T)) %>% as.data.table
rho.commongenes.summaries = all.rho.commongenes %>% group_by(Prediction.Weights) %>% summarize(rho = mean(Correlation, na.rm = T)) %>% as.data.table

# save Rdata file
save.image(rdata.path)
