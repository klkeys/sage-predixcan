#!/usr/bin/env bash                # -- what is the language of this shell?
#                                  # -- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash                    # -- the shell for the job
#$ -r y                            # -- tell the system that if a job crashes, it should be restarted
#$ -j y                            # -- tell the system that the STDERR and STDOUT should be joined
#$ -l arch=linux-x64               # -- SGE resources (CPU type)

# script variables
weightsfile=$weightsfile
glmnetdir=$glmnetdir
newweightsfile=$newweightsfile
Rscript=$Rscript
R_glmnet_postprocess=$R_glmnet_postprocess
discard_ratio=$discard_ratio
num_pred_file=$num_pred_file
nsamples=$nsamples
predictionfile=$predictionfile
exprfile=$exprfile
out_lm_file=$out_lm_file
out_genelm_file=$out_genelm_file
logdir=$logdir

# postprocess the weights file
$Rscript $R_glmnet_postprocess $weightsfile $newweightsfile $discard_ratio $num_pred_file $nsamples $predictionfile $exprfile $out_lm_file $out_genelm_file

# query return value of previous command
RETVAL=$?

# if return value is not 0, then previous command did not exit correctly
# create a status file notifying of error
# in contrary case, notify of success
if [ "$RETVAL" -ne "0" ];
then
  echo "ERROR" > ${logdir}/status.${glmmethod}.${gene}
else
  echo "SUCCESS" > ${logdir}/status.${glmmethod}.${gene}
fi
