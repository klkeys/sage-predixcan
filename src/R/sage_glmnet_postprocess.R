#!/usr/bin/env Rscript --vanilla

suppressMessages(library(data.table))

sage.glmnet.postprocess = function(){

    # read command line arguments
    args = commandArgs(trailingOnly = TRUE)

    # parse command line arguments
    weightsfile     = args[1]
    new.weightsfile = args[2]
    discard.ratio   = as.numeric(args[3])
    num.pred.file   = args[4]
    num.samples     = as.integer(args[5])
    prediction.file = args[6]
    exprfile        = args[7]
    out.lm.file     = args[8]
    out.genelm.file = args[9]

    # open weights file
    x = fread(weightsfile)

    # ensure proper header
    my.header = c("Gene","Held_out_Sample","SNP","A1","A2","Beta")
    colnames(x) = my.header

    # subset for weights with no NA values
    x.nona = x[!is.na(x$Beta),]

    # determine how many predicted samples that each gene has
    # we will sink this into an output file and reload it as a data.table
    sink(num.pred.file)
    for (i in unique(sort(x.nona$Gene))){

        # output format is "(gene)\t(#pred samples)"
        cat(i, "\t", length(unique(sort(x.nona[x.nona$Gene == i,]$Held_out_Sample))), "\n")
    }
    sink()

    # reload the formatted file
    num.pred = fread(num.pred.file)
    colnames(num.pred) = c("Gene", "Num.Pred")

    # subset those genes with at least discard.ratio * (# samples) predicted samples
    num.pred.sub = num.pred[num.pred$Num.Pred > discard.ratio*num.samples,]

    # we can now look at the predictions for the "well-predicted" genes in num.pred.sub
    # must load prediction information first
    sage.predictions = fread(prediction.file)
    sage.predictions.sub = sage.predictions[sage.predictions$Gene %in% num.pred.sub$Gene,]

    # now load measured RNA data
    sage.rnaseq = fread(exprfile)
    my.header = c("Chromosome_Name","Start_Position","End_Position","Gene","NWD159235","NWD262884","NWD424160","NWD331495","NWD511506","NWD742593","NWD154420","NWD101012","NWD843708","NWD990899","NWD373455","NWD595401","NWD299516","NWD713489","NWD468153","NWD821729","NWD167864","NWD127715","NWD952849","NWD975237","NWD212794","NWD841275","NWD359298","NWD251674","NWD620717","NWD537216","NWD863759","NWD671985","NWD923487","NWD621320","NWD345359","NWD769653","NWD152277","NWD765179","NWD546278","NWD833668","NWD437999","NWD730618","NWD554260")
    colnames(sage.rnaseq) = my.header

    # subset the file with just the genes from sage.predictions.sub
    sage.rnaseq.sub = sage.rnaseq[sage.rnaseq$Gene %in% sage.predictions.sub$Gene,]

    # sort the genes before transposing
    setorder(sage.rnaseq.sub, Gene)
    trnaseq = data.table(t(sage.rnaseq.sub[,-c(1:4)]))

    # put column names and a column of sample names
    # make sure to sort by sample
    colnames(trnaseq) = sage.rnaseq.sub$Gene
    trnaseq = cbind(my.header[-c(1:4)], trnaseq)
    colnames(trnaseq)[1] = "SubjectID"
    setorder(trnaseq, SubjectID)

    # melt the predicted and measured expression data.tables
    pred.melted   = melt(sage.predictions.sub, id.vars = c("Gene"))
    rnaseq.melted = melt(sage.rnaseq.sub[,-c(1:3)], id.vars = "Gene")

    # standardize their column names too
    colnames(pred.melted)   = c("Gene", "SubjectID", "Predicted_Expr")
    colnames(rnaseq.melted) = c("Gene", "SubjectID", "Measured_Expr")

    # now merge the data frames
    sage.rnapred = merge(pred.melted, rnaseq.melted, by = c("Gene", "SubjectID"))
    setorder(sage.rnapred, Gene, SubjectID)

    # do two linear models
    # first regress all (Gene, SubjectID) pairs of predicted expression onto measured expression
    sink(out.lm.file)
    print(summary(lm(Predicted_Expr ~ Measured_Expr, data = sage.rnapred)))
    sink()

    # now do individual gene regressions
    sink(out.genelm.file)
    cat("Gene\tP.value\tR2\tSpearman.rho\n")
    for(gene in unique(sort(sage.rnapred$Gene))){
        sage.rnapred.sub = sage.rnapred[sage.rnapred$Gene == gene,]
        my.lm = summary(lm(Predicted_Expr ~ Measured_Expr, data = sage.rnapred.sub))
        lm.p = my.lm$coefficients[2,4]
        lm.r2 = my.lm$r.squared
        my.rho = cor.test(sage.rnapred.sub$Predicted_Expr, sage.rnapred.sub$Measured_Expr, method = "spearman")$estimate
        cat(gene, "\t", lm.p, "\t", lm.r2, "\t", my.rho, "\n")
    }
    sink()

    return()
}

sage.glmnet.postprocess()

# show any warnings
cat("any warnings?\n")
warnings()

