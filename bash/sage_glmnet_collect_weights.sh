#!/usr/bin/env bash                # -- what is the language of this shell
#$ -S /bin/bash                    # -- the shell for the job
##$ -M kevin.keys@ucsf.edu          # -- email status of this job to this address
##$ -m bes                          # -- email on beginning, end, and suspension of job
#$ -r y                            # -- tell the system that if a job crashes, it should be restarted
#$ -j y                            # -- tell the system that the STDERR and STDOUT should be joined
#$ -l mem_free=1G                  # -- submits on nodes with enough free memory (required)
#$ -l arch=linux-x64               # -- SGE resources (CPU type)
#$ -l h_rt=0:10:00                 # -- runtime limit in hours 
#
# This script gathers results from a run of compute_new_predixcan_weights.sh.
# It produces two files:
# -- $resultsfile, which compiles prediction results for the current glmnet run;
# -- $weightsfile, which compiles nonzero prediction weights for the same glmnet run
#
# Call this script via qsub from compute_new_predixcan_weights.sh as follows:
#
# qsub -N glmnet.collect.results.${glmmethod} \
#      -o $logdir \
#      -e $logdir \
#      -v resultsfile=$resultsfile,weightfile=$weightfile,glmnetdir=$glmnetdir,glmmethod=$glmmethod,logdir=$logdir \
#      $glmnetdir/sage_glmnet_collect_weights.sh
#
# coded by Kevin L. Keys (2017)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# user limits: -c max size of core files created
date
hostname
ulimit -c 0 

# script variables from qsub
weightsfile=$weightsfile
glmnetdir=$glmnetdir
glmmethod=$glmmethod
logdir=$logdir
lambdafile=$lambdafile
predictionfile=$predictionfile
outdir=$outdir
resultsdir=$resultsdir
resultssubdir=$resultssubdir

# check variables for debugging
#echo "resultsfile = $resultsfile"
#echo "weightsfile = $weightsfile"
#echo "glmnetdir = $glmnetdir"
#echo "glmmethod = $glmmethod"
#echo "logdir = $logdir"

# put header on prediction file
# then append results to file
rm -f $predictionfile
touch $predictionfile
echo -e "Gene\tNWD101012\tNWD127715\tNWD152277\tNWD154420\tNWD159235\tNWD167864\tNWD212794\tNWD251674\tNWD262884\tNWD299516\tNWD331495\tNWD345359\tNWD359298\tNWD373455\tNWD424160\tNWD437999\tNWD468153\tNWD511506\tNWD537216\tNWD546278\tNWD554260\tNWD595401\tNWD620717\tNWD621320\tNWD671985\tNWD713489\tNWD730618\tNWD742593\tNWD765179\tNWD769653\tNWD821729\tNWD833668\tNWD841275\tNWD843708\tNWD863759\tNWD923487\tNWD952849\tNWD975237\tNWD990899" > $predictionfile
cat ${resultssubdir}/sage_predictions_${glmmethod}_ENSG00000* | sort >> $predictionfile

# check previous command
RETVAL=$?

# report on previous command
echo "exit status after making prediction file: $RETVAL"

# put header on lambda file
# then computed weights to file
# ensure when concatenating results that we discard the individual file headers
rm -f $lambdafile
touch $lambdafile
echo -e "Gene\tHeld_out_Sample\tMean_MSE\tLambda\tPredicted_Expr" > $lambdafile
cat ${resultssubdir}/sage_lambdas_${glmmethod}_ENSG00000* | grep -v "Gene" | sort >> $lambdafile

# check previous command
# compound with penultimate one
let "RETVAL+=$?"

# report on previous commands
echo "exit status after making prediction,lambda files: $RETVAL"

# finally, compile weights file
rm -f $weightsfile
touch $weightsfile
echo -e "Gene\tHeld_out_Sample\tSNP\tA1\tA2\tBeta" > $weightsfile
cat ${resultssubdir}/sage_weights_${glmmethod}_ENSG00000* | grep -v "Gene" | sort >> $weightsfile

# report on previous commands
echo "exit status after making prediction,lambda,weights files: $RETVAL"

# now it is safe to clean up scratch
rm -rf $outdir

if [ "$RETVAL" -ne "0" ];
then
  echo "ERROR" > ${logdir}/status.collectweights.${glmmethod}
else
  echo "SUCCESS" > ${logdir}/status.collectweights.${glmmethod}
fi
