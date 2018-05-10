#!/usr/bin/env Rscript --vanilla
#
# This script runs elastic net on the SAGE RNA-Seq data in order to compute new PrediXcan weights.
#
# coded by Kevin L. Keys (2017)
# heavily modified from original PrediXcan script by Heather E. Wheeler (2015-02-02)
# https://github.com/hakyim/PrediXcan/blob/master/Paper-Scripts/Heather/DGN-calc-weights/01_imputedDGN-WB_CV_elasticNet.r
##see runscripts/run_01_imputedDGN-WB_CV_elasticNet_chr*sh and qsub.txt for tarbell job submission scripts

###############################################
### Libraries
###############################################

library(methods)
library(glmnet)
library(data.table)

# error tracking
options(show.error.locations = TRUE)

# main script function
# it is executed at end
compute.new.weights = function(){

    # parse command line arguments
    args = commandArgs(trailingOnly=T)

    # script subroutines
    "%&%" = function(a,b) paste(a,b,sep="")

    ###############################################
    ### Directories & Variables
    ###############################################

    # parse command-line arguments
    x.df.path = args[1] # points to genotype dosages
    y.df.path = args[2] # points to BED file of expression levels
    gene.name = args[3] # name of current gene
    predout   = args[4] # points to output file for glmnet CV predictions 
    lambdaout = args[5] # points to output file for saving lambdas and held-out sample
    alpha     = args[6] # elastic net mixing parameter; 0 --> ridge, 1 --> lasso
    outpath   = args[7] # points to output file for saving nonzero betas from CV 
    bim.path  = args[8] # name of PLINK BIM corresponding to current gene

    # other script variables
    #alpha = 0.5  # elastic net mixing parameter; 0 --> ridge, 1 --> lasso
    #alpha = 1  # elastic net mixing parameter; 0 --> ridge, 1 --> lasso
    seed  = 2018 # fix starting position for CV

    # which gene are we analyzing?
    cat("Analyzing gene", gene.name, "\n")

    # load data
    # data frame x.df contains the genotype dosages for the current gene
    # data frame y.df contains the expression levels
    # N.B. this script assumes that y.df is a preformatted UCSC BED with -/+ 1Mb added to gene start/end
    # we parse data X and phenotype y from these data frames
    x.df = fread(x.df.path)
    y.df = fread(y.df.path)

    # forcibly reorder rows of x.df
    setorder(x.df, IID)

    # add header back to y.df
    my.header = "Chromosome_Name Start_Position End_Position Gene NWD159235 NWD262884 NWD424160 NWD331495 NWD511506 NWD742593 NWD154420 NWD101012 NWD843708 NWD990899 NWD373455 NWD595401 NWD299516 NWD713489 NWD468153 NWD821729 NWD167864 NWD127715 NWD952849 NWD975237 NWD212794 NWD841275 NWD359298 NWD251674 NWD620717 NWD537216 NWD863759 NWD671985 NWD923487 NWD621320 NWD345359 NWD769653 NWD152277 NWD765179 NWD546278 NWD833668 NWD437999 NWD730618 NWD554260"
    colnames(y.df) = strsplit(my.header, " ")[[1]] #<-- strsplit () gives 'list' output, must extract the character vector from it

    # use new column names to forcibly reorder y.df
    setorder(y.df, Gene)

    # now make a transpose
    ty = t(y.df[,-c(1:4)])
    #ty = ty[order(rownames(ty)),]
    colnames(ty) = y.df$Gene

    # add sample names to ty for reordering
    ty = ty[order(row.names(ty)),]

    # preallocate a vector to store predicted expression levels
    n = dim(ty)[1]
    y.pred = matrix(NA,n,1)

    # prepare the gene output file that will contain the prediction and lambda for each held-out sample
    # put a header on it for bookkeeping
    # the pipeline will later discard the header when forming a unified prediction file, one line per prediction
    lambdaarray = c("Gene", "Held_out_sample", "Mean_MSE", "Lambda", "Prediction")
    write(lambdaarray, file = lambdaout, ncolumns = 5, sep = "\t", append = FALSE)

    # do same for weights (betas)
    betafile = c("Gene", "Held_out_sample", "SNP", "A1", "A2", "Beta")
    write(betafile, file = outpath, ncolumns = 5, append = FALSE, sep = "\t")

    for (held.out.sample in 1:n){

        # get NWDID of held-out sample
        held.out.sample.name = x.df$IID[held.out.sample]

        # we expect x.df to have a PLINK RAW format
        # thus, the FID/IIDs are the same,
        # and the dosages start at column 6
        # discard all other columns since we do not need them here
        tx        = data.frame(x.df[-held.out.sample, -c(1:6)])
        x.heldout = as.double(x.df[held.out.sample, -c(1:6)])
        names(x.heldout) = colnames(tx)
        y.heldout = ty[held.out.sample,] 
        ty.train  = ty[-held.out.sample,]

        # must "impute" missing dosages
        # following Donglei's lead, use 2*MAF
        # use only observed dosages; this means that means are not necessarily divided by full sample size 39!
        mafs.missing = 2 * apply(tx, 2, mean, na.rm=T)
        #for (i in c(1:dim(tx)[2])){
        #    tx[is.na(tx[,i]),i] = mafs.missing[i]
        #}
        for (i in seq_along(tx)){
            set(tx, j=i, value = as.numeric(tx[[i]]))
            set(tx, i = which( is.na(tx[[i]]) ), j=i, value = mafs.missing[i])
        }
        na.heldout = is.na(x.heldout)
        x.heldout[na.heldout] = mafs.missing[na.heldout]
        X = data.matrix(tx) # as opposed to as.matrix()...? see https://stackoverflow.com/questions/8458233/r-glmnet-as-matrix-error-message
        row.names(X) = x.df$IID[-held.out.sample]

        # pull SNPs with at least 1 minor allele
        # also forcibly reorder the rows (samples) of X
        # we did this for the expression levels ty too
        minorsnps = subset(colMeans(X), colMeans(X, na.rm=TRUE)>0)
        minorsnps = names(minorsnps)
        #X = X[order(row.names(X)),minorsnps]
        X = X[,minorsnps]

        # X has 38 rows, but the expression levels have 39
        # must discard the corresponding row of y.df
        matching.samples = row.names(ty.train) %in% row.names(X)
        ty.train = subset(ty.train, matching.samples)
        X = X[matching.samples,]

        # mean-center the genotype variables
        X = scale(X, center = TRUE, scale = FALSE)

        # initialize data.frame of best betas
        bestbetas = data.frame()

        # predicted expression for held-out sample is missing until calculated otherwise
        y.heldout.expr.pred = NA

        # skip genes with < 2 cis-SNPs
        if(is.null(dim(X)) | dim(X)[2] == 0){
            bestbetas = data.frame() ###effectively skips genes with <2 cis-SNPs
        } else {
            # extract phenotype from expression matrix
            # the grep command pulls the numerical position of the gene in the matrix
            # save that integer and use for subsetting
            y.idx = grep(gene.name, y.df$Gene) 
            #y.idx = gene.name
            if (length(y.idx) != 0){

                # this makes a numeric matrix for the phenotype
                # all missing quantities are first set to 0
                #missing.idx = as.vector(is.na(y.df[y.idx,]))
                y = as.vector(ty.train[,gene.name])
                y[is.na(y)] = 0

                if (length(y) > 0){
                    # we must do leave-one-out cross-validation (LOOCV)
                    # in glmnet, we do k-fold CV, where k = number of samples
                    k = dim(X)[1]

                    # want to prevent unusual error caused by random generation of lambda values:
                    # https://github.com/lmweber/glmnet-error-example/blob/master/glmnet_error_example.R
                    # therefore, find lambda.max per glmnet default and set lambdas manually
                    #mysd = function(z) sqrt( sum( (z - mean(z))^2) / length(z) )
                    #sx   = as.matrix(scale(X, scale = apply(X, 2, mysd) ) )
                    sx   = as.matrix(scale(X, center = TRUE, scale = TRUE))
                    sy   = as.vector(scale(y, center = TRUE, scale = TRUE))
                    lam.max = max(abs(colSums(sx*sy))) / n
                    nlambda = 100
                    #lambda  = exp(seq(log(lam.max), log(0.01), length.out=nlambda))
                    lambda  = exp(seq(log(0.01), log(lam.max), length.out=nlambda))
                    penalty = rep(1, times = dim(X)[2])

                    # cross-validate!
                    cat(paste("crossvalidating sample ", held.out.sample.name, "...", sep = "")) 
                    #my.glmnet = cv.glmnet(x=X, y=y, nfolds=k, alpha=alpha, keep=TRUE, grouped=FALSE, dfmax=n, pmax=n, lambda=lambda, penalty.factor=penalty) # cannot group if <3 obs per fold, which is the case with LOOCV
                    #my.glmnet = cv.glmnet(x=X, y=y, nfolds=k, family = "gaussian", alpha=alpha, keep=TRUE, grouped=FALSE, dfmax=n, pmax=n, lambda=lambda)#, penalty.factor=penalty)
                    my.glmnet = cv.glmnet(x=X, y=y, nfolds=k, family = "gaussian", alpha=alpha, keep=TRUE, grouped = FALSE, dfmax=n, pmax=n, lambda.min.ratio = 0.01, nlambda=nlambda)#, penalty.factor=penalty)
                    cat("done\n")

                    # pull info to find best lambda
                    fit.df = data.frame("cvm" = my.glmnet$cvm, "lambda" = my.glmnet$lambda, "nrow" = 1:length(my.glmnet$cvm)) 

                    # parse results
                    best.lam    = fit.df[which.min(fit.df[,1]),]                         # use which.min or which.max depending on cv measure (MSE min, AUC max, ...)
                    cvm.best    = best.lam[,1]                                           # best crossvalidated mean squared error
                    lambda.best = best.lam[,2]                                           # corresponding λ for best MSE
                    nrow.best   = best.lam[,3]                                           # position of best λ in cv.glmnet output
                    ret         = as.data.frame(my.glmnet$glmnet.fit$beta[,nrow.best])   # get βs from best λ 
                    ret[ret == 0.0] = NA
                    bestbetas   = as.vector(ret[which(!is.na(ret)),])                    # vector of nonzero βs 
                    names(bestbetas) = rownames(ret)[which(!is.na(ret))]
                    pred.mat    = as.matrix(my.glmnet$fit.preval[,nrow.best])            # pull out predictions at best λ 

                    # we want to use the fitted glmnet model to predict the held-out sample
                    # normally we would use the lambda that provides the lowest prediction error (best.lam)
                    # but that can overfit the model if the number of nonzero coefficients exceeds dim(X)[1]
                    # to compensate somewhat, we can choose either (1) the best lambda, or (2) the smallest lambda yielding < dim(X)[1] nonzero coefficients
                    ##s = min(max(which(my.glmnet$nzero < dim(X)[1])), which.min(my.glmnet$lambda))
                    #s = min(max(which(my.glmnet$nzero < dim(X)[1])), which.min(my.glmnet$lambda))

                    # make prediction on held-out sample
                    # remember to only use SNPs that met minor allele threshold
                    # can use s = lambda.min to use best predictive λ directly
                    y.heldout.expr.pred = predict(my.glmnet, newx=t(as.matrix(x.heldout[minorsnps])), s = my.glmnet$lambda.min)

                    # save information about each internal CV fold
                    # order is gene, NWDID of held out sample, cvm of internal LOOCV, corresponding lambda 
                    lambdaarray = c(gene.name, held.out.sample.name, my.glmnet$cvm[which.min(my.glmnet$cvm)], my.glmnet$lambda.min, y.heldout.expr.pred)
                    write(lambdaarray, file = lambdaout, ncolumns = 5, sep = "\t", append = TRUE)

                    # record results when glmnet finds an eQTL and that eQTL is in our genotype set
                    # otherwise record missing values
                    betafile = c(gene.name, held.out.sample.name, NA, NA, NA, NA)
                    if(length(bestbetas) > 0 & !is.null(dim(pred.mat))){
                        # output best shrunken betas for PrediXcan
                        # output format: "gene","SNP","refAllele","effectAllele","beta"
                        # note that this makes explicit reference to rsid format from PLINK BIM from SAGE merged LATplus array
                        bestbetalist = names(bestbetas) # next lines format entries of this list to match BIM
                        bestbetalist = gsub("X", "", bestbetalist)
                        bestbetalist = gsub(".", ":", bestbetalist, fixed = TRUE)
                        bestbetalist = gsub("_.*", "", bestbetalist, perl = TRUE)
                        bimfile      = fread(bim.path)
                        bimfile$V2   = gsub("-", ":", bimfile$V2) # needed for SAGE LAT-LAT+ merged genotypes
                        bestbetainfo = bimfile[bimfile$V2 %in% bestbetalist,]
                        betatable    = as.matrix(cbind(bestbetainfo, bestbetas))
                        betafile     = cbind(gene.name, held.out.sample.name, betatable[,2], betatable[,5], betatable[,6], betatable[,7])
                    }
                    # save betas to file
                    # t() necessary for correct output from write() function
                    write(t(betafile), file = outpath, ncolumns = 6, append = T, sep = "\t")

                } else{
                    cat(paste("no expression data for", gene.name, "\nMarking missing results for ", x.df$IID[held.out.sample], "\n"))
                    lambdaarray = c(gene.name, x.df$IID[held.out.sample], NA, -1.0, NA) 
                    write(lambdaarray, file = lambdaout, ncolumns = 5, sep = "\t", append = TRUE)
                }
            } else{
                ### NOTA BENE! eventually change this to formatted output
                cat(paste("no expression data for", gene.name, "\nMarking missing results for ", x.df$IID[held.out.sample], "\n"))
                lambdaarray = c(gene.name, x.df$IID[held.out.sample], NA, -1.0, NA) 
                write(lambdaarray, file = lambdaout, ncolumns = 5, sep = "\t", append = TRUE)
            }

        # save predicted expression level
        y.pred[held.out.sample] = y.heldout.expr.pred

        }

    }
    # write y.pred to file
    predictionarray = c(gene.name, t(y.pred))
    write(predictionarray, file = predout, ncolumns = n + 1, sep = "\t", append = FALSE)

    return()
}

# shut up!
#oldw <- getOption("warn")
#options(warn = -1)

# run function
compute.new.weights()


# any warnings?
cat("any warnings?\n")
warnings()

cat("done executing sage_gtex_compute_new_weights.R\n")


# un-shutup
#options(warn = oldw)
