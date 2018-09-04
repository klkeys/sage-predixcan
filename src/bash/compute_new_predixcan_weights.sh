#!/usr/bin/env bash
#
# This script computes PrediXcan weights from SAGE transcriptome data.
#
# Call:
#
# ./compute_new_predixcan_weights.sh $ALPHA
#
# where
# -- $ALPHA = "0.0", "0.5", or "1.0", used by glmnet to determine the training algorithm
# -- $DEFAULT_R2 is any nonempty nonzero value, used to compute R2/p-values with default GTEx weights (or not)
#
# coded by Kevin L. Keys (2017)

# parse command line arguments
alpha=$1
r2_default=$2

# set method based on alpha
# only checks the three admissible values for alpha
glmmethod=""
if [[ "$alpha" == "1.0" ]]; then
    glmmethod="lasso";
fi
if [[ "$alpha" == "0.5" ]]; then
    glmmethod="elasticnet";
fi
if [[ "$alpha" == "0.0" ]]; then
    glmmethod="ridge";
fi

# make sure that glmmethod was set
# $glmmethod is empty if alpha was not given an acceptable value
if [[ "$glmmethod" == "" ]]; then
    echo -e "usage:\n\tsource compute_new_predixcan_weights.sh \$ALPHA [\$DEFAULT]\n"
    echo -e "where \$ALPHA = {0.0, 0.5, 1.0} (required)\nand DEFAULT != 0 (optional, if calulation with default weights desired)\n"
    exit 1;
fi

# script static directories
MYHOME="/netapp/home/kkeys"
rnaseqdir="${MYHOME}/gala_sage/rnaseq"
gctadir="${rnaseqdir}/gcta"
bindir="${MYHOME}/bin"
datadir="${rnaseqdir}/data"
glmnetdir="${rnaseqdir}/glmnet"
genotypedir="${MYHOME}/gala_sage/genotypes/mergedLAT-LATP/SAGE_merge"
logdir="${glmnetdir}/log"
#outdir="${glmnetdir}/${glmmethod}/genes"
outdir="/scratch/kkeys/${glmmethod}/genes"
#imputegenodir="${MYHOME}/gala_sage/genotypes/SAGE/IMPUTE_HRC_1.1_PLINK"
imputegenodir="${MYHOME}/gala_sage/genotypes/SAGE/IMPUTE2_1KGv3_incommon_with_WGS"
codedir="${MYHOME}/gala_sage/code"
resultsdir="${glmnetdir}/${glmmethod}"
resultssubdir="${resultsdir}/results"

# make output and results directories, in case they doesn't exist
mkdir -p $outdir
mkdir -p $resultssubdir

# script files
exprfile="${datadir}/sage_39_wgs_for_rnaseq_expression_sorted_headered.bed"
#predictionfile="${glmnetdir}/sage_${glmmethod}_predictions.txt"
#lambdafile="${glmnetdir}/sage_${glmmethod}_lambdas.txt"
#weightsfile="${glmnetdir}/sage_${glmmethod}_weights.txt"
ucsc_snpfile="${glmnetdir}/sage_rnaseq_UCSC_allsnps_merged_wgs_rsid.txt"
r2resultsfile="${glmnetdir}/sage_gtex_r2_defaultweights.txt"
#sagegenopfx="${genotypedir}/SAGE_mergedLAT-LATP_030816_rsID"
genelist="${gctadir}/genelist_plusmin1Mb.txt"
sage_lab2nwd="${datadir}/sage_lab2nwd.txt"
subjectids="${datadir}/sage_39_subjectids.txt"
nwdids="${datadir}/sage_39_nwdids.txt"
predictionfile="${resultsdir}/sage_${glmmethod}_predictions.txt"
lambdafile="${resultsdir}/sage_${glmmethod}_lambdas.txt"
weightsfile="${resultsdir}/sage_${glmmethod}_weights.txt"
num_pred_file="${resultsdir}/sage_${glmmethod}_numpred.txt"
newweightsfile="${resultsdir}/sage_${glmmethod}_weights_noNA.txt"
out_lm_file="${resultsdir}/sage_${glmmethod}_lm_predvmeas_results.txt"
out_genelm_file="${resultsdir}/sage_${glmmethod}_genelm_predvmeas_results.txt"

# script locations
R_compute_new_weights="${codedir}/sage_gtex_compute_new_weights.R"
R_compute_r2="${codedir}/sage_gtex_compute_r2.R"
R_glmnet_postprocess="${codedir}/sage_glmnet_postprocess.R"

# script binaries
PLINK="${bindir}/plink"
Rscript="${bindir}/Rscript"
#Rscript="/usr/bin/Rscript"

# script variables
maf="0.05"
hwe="0.03"
geno="0.03"
nthreads=1
memory_limit="2G"
memory_limit_mb="2000" # manually coordinate this with $memory_limit!!!
scratch_memory="2G"
h_rt="0:29:59"
discard_ratio="0.5" # this is the desired min % of samples with LOOCV predictions, used in postprocessing
nsamples="39"

# -------------------- #
# start script
echo "Start time: $(date)"

# -------------------- #
# make necessary output directories
if [[ ! -d "$logdir" ]]; then mkdir -p $logdir; fi
#if [[ ! -d "$outdir" ]]; then mkdir -p $outdir; fi


# -------------------- #
## do the same for default weights file
#if [[ ! -z "$r2_default" ]]; then
#    if [[ "$r2_default" != "0" ]]; then 
#        touch $r2resultsfile
#        echo -e "gene\tn.snps\tR2\tpval" > $r2resultsfile
#    fi
#fi



# -------------------- #
# compute new weights for each gene
# this entails running glmnet on a PLINK RAW genotype dosage file
# and extracting the expression level from the UCSC BED file stored at $exprfile

# how many genes are in this list?
nGenes=$(wc -l $genelist | cut -f 1 -d " ")

### for testing
#nGenes=1


# -N glmnet.${gene}    ##--- wait until array job lowComplex_byPop.${PR,AF,MX}
# -e $glmnetdir/log    ##--- where to put error files
# -o $glmnetdir/log    ##--- where to put output files
# -v ...               ##--- variables to pass to script
# -t ...               ##--- array job components, one for each gene
# -l mem_free=1G       ## -- submits on nodes with enough free memory, e.g. 1 gigabyte (required)
# -l h_rt=0:29:59      ## -- runtime limit in hours 
# 
#qsub -N glmnet.${glmmethod} \
#     -v genelist=$genelist,subjectids=$subjectids,sage_lab2nwd=$sage_lab2nwd,Rscript=$Rscript,R_compute_new_weights=$R_compute_new_weights,exprfile=$exprfile,logdir=$logdir,alpha=$alpha,gctadir=$gctadir,glmmethod=$glmmethod,outdir=$outdir,imputegenodir=$imputegenodir,maf=$maf,hwe=$hwe,nthreads=$nthreads,PLINK=$PLINK,memory_limit_mb=$memory_limit_mb,resultssubdir=$resultssubdir,resultsdir=$resultsdir,nwdids=$nwdids \
#     -t 1-$nGenes \
#     -e $logdir \
#     -o $logdir \
#     -l mem_free=$memory_limit \
#     -l scratch=$scratch_memory \
#     -l h_rt=$h_rt \
#     ${codedir}/qsub_run_sage_glmnet.sh
#     #-t 1112-1112 \
#     #-v genelist=$genelist,sagegenopfx=$sagegenopfx,subjectids=$subjectids,sage_lab2nwd=$sage_lab2nwd,Rscript=$Rscript,R_compute_new_weights=$R_compute_new_weights,exprfile=$exprfile,logdir=$logdir,alpha=$alpha,gctadir=$gctadir,glmmethod=$glmmethod,outdir=$outdir,imputegenodir=$imputegenodir,maf=$maf,hwe=$hwe,nthreads=$nthreads,PLINK=$PLINK,memory_limit_mb=$memory_limit_mb \
#
#
## collect output
#qsub -N glmnet.${glmmethod}.collect.weights \
#     -hold_jid glmnet.${glmmethod} \
#     -v weightsfile=$weightsfile,glmnetdir=$glmnetdir,glmmethod=$glmmethod,logdir=$logdir,predictionfile=$predictionfile,lambdafile=$lambdafile,outdir=$outdir,resultsdir=$resultsdir,resultssubdir=$resultssubdir \
#     -o $logdir \
#     -e $logdir \
#     -l mem_free=$memory_limit \
#     -l h_rt=$h_rt \
#     ${codedir}/sage_glmnet_collect_weights.sh

# process output file
qsub -N glmnet.${glmmethod}.postprocess \
    -hold_jid glmnet.${glmmethod},glmnet.${glmmethod}.collect.weights \
    -v weightsfile=$weightsfile,glmnetdir=$glmnetdir,newweightsfile=$newweightsfile,Rscript=$Rscript,R_glmnet_postprocess=$R_glmnet_postprocess,discard_ratio=$discard_ratio,num_pred_file=$num_pred_file,nsamples=$nsamples,predictionfile=$predictionfile,exprfile=$exprfile,out_lm_file=$out_lm_file,out_genelm_file=$out_genelm_file,logdir=$logdir \
    -o $logdir \
    -e $logdir \
     -l mem_free=$memory_limit \
     -l h_rt=$h_rt \
    ${codedir}/qsub_sage_glmnet_postprocess.sh

# end script
echo "End time: $(date)"
