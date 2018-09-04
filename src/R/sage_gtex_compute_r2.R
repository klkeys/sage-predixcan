#!/usr/bin/env Rscript --vanilla
#
# This script computes R2 for SAGE genes using PrediXcan default weights.
# The script is modeled after sage_gtex_compute_new_weights.R but does not train 
#
# coded by Kevin L. Keys (2017)

###############################################
### Libraries
###############################################

suppressMessages(library(methods))
suppressMessages(library(data.table))

# error tracking
options(show.error.locations = TRUE)

# main script function, executed at end
# this function uses default weights from GTEx to compute R2s using the genes
# it is similar in spirit to the previous one, but it does not need to use LASSO
# instead, it loads the weights (betas) from file
compute.gtex.r2 = function(){

    # parse command line arguments
    args = commandArgs(trailingOnly=T)

    # parse command-line arguments
    x.df.path     = args[1] # points to genotype dosages
    y.df.path     = args[2] # points to BED file of expression levels
    gene.name     = args[3] # name of current gene
    ucsc.path     = args[4] # name of merged UCSC SNP / weight file 
    outpath       = args[5] # points to output file for R2 results
    genotype.data = args[6] # any nonempty value here will be rendered to TRUE

    if (!is.na(genotype.data)){
        genotype.data = TRUE
    } else{
        genotype.data = FALSE
    }

    # load data
    # data frame x.df contains the genotype dosages for the current gene
    # data frame y.df contains the expression levels
    # N.B. this script assumes that y.df is a preformatted UCSC BED with -/+ 1Mb added to gene start/end
    # we parse data X and phenotype y from these data frames
    x.df = fread(x.df.path)
    y.df = fread(y.df.path)

    # add header back to y.df
    my.header = "Chromosome_Name Start_Position End_Position Gene NWD159235 NWD262884 NWD424160 NWD331495 NWD511506 NWD742593 NWD154420 NWD101012 NWD843708 NWD990899 NWD373455 NWD595401 NWD299516 NWD713489 NWD468153 NWD821729 NWD167864 NWD127715 NWD952849 NWD975237 NWD212794 NWD841275 NWD359298 NWD251674 NWD620717 NWD537216 NWD863759 NWD671985 NWD923487 NWD621320 NWD345359 NWD769653 NWD152277 NWD765179 NWD546278 NWD833668 NWD437999 NWD730618 NWD554260"
    colnames(y.df) = strsplit(my.header, " ")[[1]] #<-- strsplit () gives 'list' output, must extract the character vector from it
    ty = t(y.df[,-c(1:4)])
    ty = ty[order(rownames(ty)),]
    colnames(ty) = y.df$Gene

    # we expect x.df to have a PLINK RAW format
    # thus, the FID/IIDs are the same,
    # and the dosages start at column 6
    # discard all other columns since we do not need them here
    tx = data.frame(x.df[,-c(1:6)])

    # must "impute" missing dosages
    # following Donglei's lead, use 2*MAF
    # use only observed dosages; this means that means are not necessarily divided by full sample size 39!
    mafs.missing = 2 * apply(tx, 2, mean, na.rm=T)
    for (i in c(1:dim(tx)[2])){
        tx[is.na(tx[,i]),i] = mafs.missing[i]
    }
    X = data.matrix(tx) # as opposed to as.matrix()...? see https://stackoverflow.com/questions/8458233/r-glmnet-as-matrix-error-message
    row.names(X) = x.df$IID

    # pull SNPs with at least 1 minor allele
    minorsnps = subset(colMeans(X), colMeans(X, na.rm=TRUE)>0)
    minorsnps = names(minorsnps)
    X = X[order(row.names(X)),minorsnps]

    # X has 38 rows, but the expression levels have 39
    # must discard the corresponding row of y.df
    ty = subset(ty, row.names(ty) %in% row.names(X))

    # which gene are we analyzing?
    cat("Analyzing gene", gene.name, "\n")

    # this will only execute if the gene has < 2 cis SNPs
    # we then proceed to extract expression levels and compute r2
    ucsc.ids     = c()
    y            = c()
    my.gene.cols = c()
    pred.mat     = matrix(0, 1, 0) 
    if(!is.null(dim(X)) & dim(X)[2] > 0){

        # extract phenotype from expression matrix
        # the grep command pulls the numerical position of the gene in the matrix
        # save that integer and use for subsetting
        y.idx = grep(gene.name, y.df$Gene) 

        if (length(y.idx) != 0){

            # this makes a numeric matrix for the phenotype
            # all missing quantities are first set to 0
            #missing.idx = as.vector(is.na(y.df[y.idx,]))
            y = as.vector(ty[,gene.name])
            y[is.na(y)] = 0

            # load SNP weight file 
            ucsc.snps = fread(ucsc.path)

            # subset the SNPs to match the current gene
#            if (genotype.data) {
#                ucsc.ids = ucsc.snps[ucsc.snps$gene == gene.name,]$id
#            } else{
                ucsc.ids = ucsc.snps[ucsc.snps$gene == gene.name,]$rsid
#            }

            if(length(ucsc.ids) > 0){
                # reformat column names of X to match UCSC SNPs
                colnames(X) = gsub("X", "", colnames(X))
                colnames(X) = gsub("_.*", "", colnames(X), perl = TRUE)
                colnames(X) = gsub(".", ":", colnames(X), fixed = TRUE)

                # get columns from UCSC snp table that appear in gene
                my.gene.cols = colnames(X) %in% ucsc.ids

                if(length(my.gene.cols) > 0){
                    # pull out genotypes corresponding to nonzero betas 
                    pred.mat = as.matrix(X[,my.gene.cols])
                } else{
                    cat(paste("no matching SNPs found for", gene.name, "\n"))
                }
            } else{
                cat(paste("no snp model found for", gene.name, "\n"))
            }
        } else{
            cat(paste("no expression data for", gene.name, "\n"))
        }
    }

    # initialize default resultsarray
    resultsarray = c(gene.name,0,NA,NA)
    # record results when glmnet finds an eQTL
    # otherwise record missing values
#    print(length(ucsc.ids))
#    print(length(ucsc.ids) > 0)
#    print(length(my.gene.cols))
#    print(length(my.gene.cols) > 0)
#    print(dim(pred.mat)[2])
#    print(dim(pred.mat)[2] > 0)
    if((length(ucsc.ids) > 0) & (length(my.gene.cols) > 0) & (dim(pred.mat)[2] > 0)){
        res  = summary(lm(y ~ pred.mat))
        rsq  = res$r.squared
        pval = res$coef[2,4]

        # print results
        resultsarray = c(gene.name, length(ucsc.ids), rsq, pval)
    }

    # write results to file
    write(resultsarray, file = outpath, ncolumns = 4, append = T, sep = "\t")

    return()
}

# shut up!
oldw <- getOption("warn")
options(warn = -1)

# run function
compute.gtex.r2()

# un-shutup
options(warn = oldw)
