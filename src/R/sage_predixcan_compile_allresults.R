#!/usr/bin/env Rscript --vanilla

library(data.table)
library(dplyr)
library(broom)
library(purrr)
library(ggplot2)
library(dunn.test)

# switches
compile.results    = FALSE
compute.statistics = TRUE

# file paths
all.results.file = "sage_predixcan_allresults_alltrainsets.txt"

# compile imputation results from system files?
# these files must come from QB3; check there first if they are missing
if (compile.results) {

    # load all results files that we want to compile into single output
    gtex6.results = fread("~/gala_sage/rnaseq/predixcan/GTEx/GTEx_v6p_predixcan_prediction_ALLCHR_predicted_expression_melted.txt")
    gtex7.results = fread("~/gala_sage/rnaseq/predixcan/GTEx/GTEx_v7_predixcan_prediction_ALLCHR_predicted_expression_melted.txt")
    dgn.results = fread("~/gala_sage/rnaseq/predixcan/DGN/DGN_predixcan_prediction_ALLCHR_predicted_expression_melted.txt")
    sage.results.df = fread("~/gala_sage/rnaseq/glmnet/elasticnet/sage39/sage_elasticnet_predictions.txt")

    # melt the sage results into a gene-major data.table
    sage.results = melt(sage.results.df, id.vars = c("Gene"))

    # rename column names as necessary
    # rename subject identifiers to SubjectID for consistency (some files have IID, others don't)
    colnames(gtex6.results)[1] = "SubjectID"
    colnames(gtex7.results)[1] = "SubjectID"
    colnames(sage.results) = c("Gene", "SubjectID", "SAGE_predicted_expr")

    # load normalized gene expression for SAGE from UCSC BED file
    # must add SAGE header back since it is commented in original file
    sage.measured.expr = fread("~/gala_sage/rnaseq/data/sage_39_wgs_for_rnaseq_expression_sorted_headered.bed")
    my.header = strsplit("Chromosome_Name Start_Position End_Position Gene NWD159235 NWD262884 NWD424160 NWD331495 NWD511506 NWD742593 NWD154420 NWD101012 NWD843708 NWD990899 NWD373455 NWD595401 NWD299516 NWD713489 NWD468153 NWD821729 NWD167864 NWD127715 NWD952849 NWD975237 NWD212794 NWD841275 NWD359298 NWD251674 NWD620717 NWD537216 NWD863759 NWD671985 NWD923487 NWD621320 NWD345359 NWD769653 NWD152277 NWD765179 NWD546278 NWD833668 NWD437999 NWD730618 NWD554260", " ")[[1]]
    colnames(sage.measured.expr) = my.header

    # ditch the chr,start,end information for each gene since we don't need it here
    sage.measured.expr = sage.measured.expr[,-c(1:3)]

    # melt the resulting data frame
    sage.measured.expr.melt = melt(sage.measured.expr, id.vars = c("Gene"))
    colnames(sage.measured.expr.melt) = c("Gene", "SubjectID", "Measured_Expr")

    # begin merging
    all.results = merge(sage.measured.expr.melt, sage.results, by = c("Gene", "SubjectID"), all = TRUE)
    all.results = merge(all.results, gtex6.results, by = c("Gene", "SubjectID"), all = TRUE)
    all.results = merge(all.results, gtex6.results, by = c("Gene", "SubjectID"), all = TRUE)
    all.results = merge(all.results, gtex7.results, by = c("Gene", "SubjectID"), all = TRUE)
    colnames(dgn.results)[1] = "SubjectID"
    all.results = merge(all.results, dgn.results, by = c("Gene", "SubjectID"), all = TRUE)

    # write merged results to file
    fwrite(x = all.results, file = all.results.file, col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
}

# compute summary statistics over imputation results?
if (compute.statistics) {

    # set color palette and boxplot labels
    my.colors = c("blue", "orange", "red", "gold")
    #my.colors = c("blue", "red", "gold")

    # gplot2 colorblind-friendly palette with grey:
    cbPalette = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    my.boxplot.labels = c("GTEx v6p", "GTEx v7", "DGN", "GTEx7 Test")

    # load output from previous step
    x = fread(all.results.file)

    # will use all imputation sets (SAGE, GTEX6, GTEX7, DGN)
    prediction.sets = c("SAGE", "GTEx_v6p", "GTEx_v7", "DGN")

    # create formulae for use when constructing linear models and correlation tests
    # online solution is to build them as lists scoped to global environment
    # for lm(), formula is Y ~ X
    # for cor.test(), formula is ~ Y + X
    # only use relevant colnames of results file as variables. first 3 are fixed as Gene, SubjectID, and Measured_Expr
    lm.formulae = lapply(colnames(x)[-c(1:3)], function(var) formula(paste0(var, " ~ Measured_Expr"), env=globalenv()))
    cor.test.formulae = lapply(colnames(x)[-c(1:3)], function(var) formula(paste0("~ ", var, " + Measured_Expr"), env=globalenv()))

    # seed a data table with the genes from the compiled results file
    summary.stats.df = data.table("Gene" = unique(sort(x$Gene)))

    # loop over imputation sets
    for (i in 1:4) {

        # extract variables for this column
        prediction.col = colnames(x)[i+3]        # name of column
        my.lm.formula  = lm.formulae[[i]]        # formula for linear model
        my.cor.formula = cor.test.formulae[[i]]  # formula for correlation test
        prediction.set = prediction.sets[i]      # name of the current prediction set

        # compute linear models and extract the coefficient of determination for each gene
        r2s = x %>%
            select(Gene, SubjectID, Measured_Expr, prediction.col) %>%
            na.omit %>%
            group_by(Gene) %>%
            do(tidy(summary(lm(my.lm.formula, data = .))$r.squared)) %>%
            select(Gene, x)
        r2s = data.table(r2s)
        colnames(r2s) = c("Gene", paste0(prediction.set, "_R2"))

        # do same as before, but here we compute Spearman correlations
        correlations = x %>%
            select(Gene, SubjectID, Measured_Expr, prediction.col) %>%
            na.omit %>%
            group_by(Gene) %>%
            do(tidy(cor.test(.[[prediction.col]], .[["Measured_Expr"]], method = "spearman")$estimate)) %>%
            select(Gene, x)
        correlations = data.table(correlations)
        colnames(correlations) = c("Gene", paste0(prediction.set, "_Spearman.rho"))

        # merge the results for this prediction set
        current.results = merge(r2s, correlations, by = "Gene", all = TRUE)

        # now merge with the overall results
        summary.stats.df = merge(summary.stats.df, current.results, by = "Gene", all = TRUE)

    }

    # save results to file
    fwrite(x = summary.stats.df, file = "sage_predixcan_summaries_genewise.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

    ### PLOT SUMMARIES OF ALL RESULTS
    # compare imputation performance from all four prediction weights
    # must first load GTEx v7 testing R2s from PredictDB
    predixcan.gtex7.r2s = fread("~/gala_sage/rnaseq/predixcan/gtex7.testR2.txt")

    # extract coefficients of determination (R2) from genewise summary results
    # make sure to rename the columns!
    all.r2.commongenes = summary.stats.df %>%
        merge(., predixcan.gtex7.r2s, by = "Gene", all = T) %>%
        select(Gene, GTEx_v6p_R2, GTEx_v7_R2, DGN_R2, GTEx7_test_R2_avg) %>%
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
    ggsave(plot = my.boxplot.r2, filename = "all.r2.commongenes.boxplot.png")
    ggsave(plot = my.violinplot.r2, filename = "all.r2.commongenes.violinplot.png")
    ggsave(plot = my.hist.r2, filename = "all.r2.commongenes.histogram.png")

    # do same, but for correlations instead of R2s
    all.rho.commongenes = summary.stats.df %>%
        merge(., predixcan.gtex7.r2s, by = "Gene", all = T) %>%
        select(Gene, GTEx_v6p_Spearman.rho, GTEx_v7_Spearman.rho, DGN_Spearman.rho, GTEx7_test_rho_avg) %>%
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

    # in addition to plots, perform statistical tests of differences in distributions
    # for pairwise tests, can use Kruskal-Wallis test
    # the Dunn test can compute all pairwise and overall KW tests
    dunn.test(all.r2.commongenes$R2, g = as.factor(all.r2.commongenes$Prediction.Weights), kw = T, method = "bh")
    #  Kruskal-Wallis rank sum test
    #
    #data: x and group
    #Kruskal-Wallis chi-squared = 781.4558, df = 3, p-value = 0
    #
    #
    #                           Comparison of x by group                            
    #                             (Benjamini-Hochberg)                              
    #Col Mean-|
    #Row Mean |     DGN_R2   GTEx_v6p   GTEx_v7_
    #---------+---------------------------------
    #GTEx_v6p |   0.431356
    #         |     0.3997
    #         |
    #GTEx_v7_ |  -0.308636  -0.739992
    #         |     0.3788     0.3445
    #         |
    #GTEx7_te |  -22.77580  -23.20715  -22.46716
    #         |    0.0000*    0.0000*    0.0000*
    #
    #alpha = 0.05
    #Reject Ho if p <= alpha/2
    dunn.test(all.rho.commongenes$Correlation, g = as.factor(all.rho.commongenes$Prediction.Weights), kw = T, method = "bh")
    #  Kruskal-Wallis rank sum test
    #
    #data: x and group
    #Kruskal-Wallis chi-squared = 1234.0534, df = 3, p-value = 0
    #
    #
    #                           Comparison of x by group                            
    #                             (Benjamini-Hochberg)                              
    #Col Mean-|
    #Row Mean |   DGN_Spea   GTEx_v6p   GTEx_v7_
    #---------+---------------------------------
    #GTEx_v6p |   0.074241
    #         |     0.5645
    #         |
    #GTEx_v7_ |  -0.016965  -0.091206
    #         |     0.4932     0.6955
    #         |
    #GTEx7_te |  -28.66358  -28.73782  -28.64662
    #         |    0.0000*    0.0000*    0.0000*
    #
    #alpha = 0.05
    #Reject Ho if p <= alpha/2

    ### FOCUS ON GTEX
    # 3 comparisons:
    # (1): extract genes with v6 predictions, and pull along any v7 predictions there. compare prediction quality
    # (2): extract genes with v7 predictions, pull along v6. compare again
    # (3): extract genes with *both* v6,v7 predictions. compare one last time
    # draw up Venn diagram with v6,v7 statistics in each case
    gtex6.genes = summary.stats.df %>% select(Gene, GTEx_v6p_R2) %>% na.omit %>% select(Gene) %>% map(sort) %>% map(unique)  ## nota bene: these objects are returned as LISTS
    gtex7.genes = summary.stats.df %>% select(Gene, GTEx_v7_R2) %>% na.omit %>% select(Gene) %>% map(sort) %>% map(unique)
    gtex.genes  = summary.stats.df %>% select(Gene, GTEx_v6p_R2, GTEx_v7_R2) %>% na.omit %>% select(Gene) %>% map(sort) %>% map(unique)

    # how many predictions in each case?
    # > sapply(c(gtex6.genes, gtex7.genes, gtex.genes), length)
    # Gene Gene Gene
    # 2382 1107  955

    # subset GTEx results
    gtex.summary = summary.stats.df %>% select(Gene, GTEx_v6p_R2, GTEx_v7_R2) %>% as.data.table

    # (1): subset to columns in gtex6.genes
    gtex6.summary = gtex.summary %>% filter(Gene %in% gtex6.genes[[1]])
    summary(gtex6.summary)
    # $GTEx_v6p_R2
    #     Min.  1st Qu.   Median  Mean  3rd Qu.     Max.
    #     0.00  0.002959 0.012330 0.028 0.036760 0.407300
    #
    # $GTEx_v7_R2
    #     Min.  1st Qu.  Median    Mean    3rd Qu.    Max.    NA's
    #     0.00  0.0031  0.0126  0.0270  0.0369  0.5497    1427

    # (2): subset to columsn in gtex7.genes
    gtex7.summary = gtex.summary %>% filter(Gene %in% gtex7.genes[[1]])
    summary(gtex7.summary)
    # $GTEx_v6p_R2
    #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
    #    0.00 0.00289 0.01220 0.02930 0.03492 0.40730     152
    #
    # $GTEx_v7_R2
    #    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
    #    0.00  0.003161 0.012870 0.027570 0.037360 0.549700

    # (3): subset to intersection of gtex6.genes, gtex.7.genes
    gtex6and7.summary = gtex.summary %>% filter(Gene %in% gtex.genes[[1]])
    summary(gtex6and7.summary)
    # $GTEx_v6p_R2
    #     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
    #     0.00  0.002892 0.012200 0.029300 0.034920 0.407300
    #
    # $GTEx_v7_R2
    #     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
    #     0.00  0.003108 0.012630 0.026980 0.036900 0.549700

    # compare this to PrediXcan R2s from GTEx v7
    gtex7.compare.r2 = gtex.summary %>%
        select(Gene, GTEx_v7_R2) %>%
        merge(., predixcan.gtex7.r2s, by = "Gene", all = T) %>%
        as.data.table %>%
        na.omit

    # save GTEx7 comparison
    fwrite(x = gtex7.compare.r2, file = "gtex7.compare.r2.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

    ### create Venn diagram to compare GTEx v6p and v7

    # plot parameters
    height = 600
    width = 600
    resolution = 1200

    # function venn.diagram() by default labels Venn circles with size of intersection
    # create new labels that instead put summary statistics
    gtex6.labels = paste(paste0("N = ", length(gtex6.genes[[1]])), paste0("v6p R2: ", round(mean(gtex6.summary$GTEx_v6p_R2, na.rm=T), digits = 3)), paste0("v7 R2: ", round(mean(gtex6.summary$GTEx_v7_R2, na.rm = T), digits = 3)), sep = "\n")
    gtex7.labels = paste(paste0("N = ", length(gtex7.genes[[1]])), paste0("v6p R2: ", round(mean(gtex7.summary$GTEx_v6p_R2, na.rm=T), digits = 3)), paste0("v7 R2: ", round(mean(gtex7.summary$GTEx_v7_R2, na.rm = T), digits = 3)), sep = "\n")
    gtex6and7.labels = paste(paste0("N = ", length(gtex.genes[[1]])), paste0("v6p R2: ", round(mean(gtex6and7.summary$GTEx_v6p_R2, na.rm=T), digits = 3)), paste0("v7 R2: ", round(mean(gtex6and7.summary$GTEx_v7_R2, na.rm = T), digits = 3)), sep="\n")

    # create Venn diagram, but do not plot it yet
    v = venn.diagram(
        list("GTEx v6p" = gtex6.genes[[1]], "GTEx v7" = gtex7.genes[[1]]),
        alpha = c(0.3, 0.3),
        fill = c("red", "blue"),
        filename = NULL,
        resolution = resolution,
        height = height,
        width = width,
        filetype = "tiff"
    )

    # change labels to reflect summary stats
    v[[5]]$label = gtex6.labels
    v[[6]]$label = gtex7.labels
    v[[7]]$label = gtex6and7.labels

    # save plot as TIFF; can load in image viewer and reexport to other formats as necessary
    grid.newpage()
    #tiff("gtex6and7.r2comparison.venn.tiff", height = height, width = width)
    png("gtex6and7.r2comparison.venn.png", height = height, width = width)
    grid.draw(v)
    dev.off()

}
