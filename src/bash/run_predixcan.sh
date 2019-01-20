#!/usr/bin/env bash
# copyright Asthma Collaboratory (2018)
# coded by Kevin L. Keys (2018)
# with substantial help from Donglei Hu

DHU_HOME="/media/BurchardRaid01/LabShare/Home/dhu"
predixcan_dir="${HOME}/gala_sage/rnaseq/predixcan"
predixcan="${DHU_HOME}/PrediXcan-master/Software/PrediXcan.py" 
imputed_genotype_path="${DHU_HOME}/bdr_seq/1441samples_refiltered/predixcan/TW_WholeBlood" 
dosage_prefix="chr"
samples_file="${DHU_HOME}/bdr_seq/1441samples_refiltered/predixcan/TW_WholeBlood/samples_1441_refiltered.txt" 
output_prefix="${predixcan_dir}/output_predixcan_gtex_1kg"
predixcan_sfx="predixcan.genotypefile.SAGE_hg19-1kgP3v5_panel.mixedpop_phased.eagle_39.pass.x-lowComplex.dp10-gq20.SNP.id.txt" # coordinate with Cesar output

# directories
MYHOME="/netapp/home/kkeys"
galadir="${MYHOME}/gala_sage"
codedir="${galadir}/code"
bindir="${MYHOME}/bin"
outdir="/scratch/kkeys/predixcan/results"
predixcanoutdir="${galadir}/rnaseq/predixcan"
dgndir="${predixcanoutdir}/DGN"
gtexdir="${predixcanoutdir}/GTEx"
mesadir="${predixcanoutdir}/MESA"
logdir="${predixcanoutdir}/log"
predixcan_genodir="${galadir}/genotypes/SAGE/EAGLE_1kgP3v5_mixedpop_incommon_with_WGS"

# file paths
bedfile_pfx="${imputegenodir}/chr" # qsub script assumes form ${bedfile_pfx}${SGE_TASK_ID}.${bedfile_sfx} 
#bedfile_sfx="SAGE_hg19-1kgP3.phased.impute2_39.pass.x-lowComplex.dp10-gq20.SNP.id" # this is most editable of "$bedfile" variables 
bedfile_sfx="SAGE_hg19-1kgP3v5_panel.mixedpop_phased.eagle_39.pass.x-lowComplex.dp10-gq20.SNP.id"
weights_file_GTEx_v6p="/netapp/home/kkeys/git/PrediXcan/Databases/TW_Whole_Blood_0.5_1KG.db"
weights_file_GTEx_v7="/netapp/home/kkeys/git/PrediXcan/Databases/gtex_v7_Whole_Blood_imputed_europeans_tw_0.5_signif.db"
weights_file_DGN="/netapp/home/kkeys/git/PrediXcan/Databases/DGN-WB_0.5.db"
weights_file_MESA_AFA="/netapp/home/kkeys/git/PrediXcan/Databases/AFA_imputed_10_peer_3_pcs_v2.db"
weights_file_MESA_AFHI="/netapp/home/kkeys/git/PrediXcan/Databases/AFHI_imputed_10_peer_3_pcs_v2.db"
weights_file_MESA_ALL="/netapp/home/kkeys/git/PrediXcan/Databases/ALL_imputed_10_peer_3_pcs__v2.db"
weights_file_MESA_CAU="/netapp/home/kkeys/git/PrediXcan/Databases/CAU_imputed_10_peer_3_pcs_v2.db"



# binaries
PLINK="${bindir}/plink"
python2="${bindir}/python2"
Rscript="${bindir}/Rscript"

# scripts
BASH_run_predixcan="${codedir}/qsub_run_predixcan.sh"
BASH_postprocess_predixcan="${codedir}/qsub_postprocess_predixcan.sh"
R_format_predixcan="${codedir}/format_plink_traw_into_predixcan.R"
R_impute_predixcan_genos="${codedir}/impute_predixcan_genotypes.R" 
R_postprocess_predixcan="${codedir}/postprocess_predixcan_results.R" 
predixcan="/netapp/home/dhu/PrediXcan-master/Software/PrediXcan.py"

# SGE variables
memory_limit="2G"
scratch_memory="2G"
h_rt="0:29:59"

# output prefixes for PrediXcan
output_pfx_DGN="${dgndir}/DGN_predixcan_prediction"
output_pfx_GTEx_v6p="${gtexdir}/GTEx_v6p_predixcan_prediction"
output_pfx_GTEx_v7="${gtexdir}/GTEx_v7_predixcan_prediction"
output_pfx_MESA_AFA="${mesadir}/MESA_AFA_predixcan_prediction"
output_pfx_MESA_AFHI="${mesadir}/MESA_AFHI_predixcan_prediction"
output_pfx_MESA_ALL="${mesadir}/MESA_ALL_predixcan_prediction"
output_pfx_MESA_CAU="${mesadir}/MESA_CAU_predixcan_prediction"



# make directories if necessary
mkdir -p $outdir
mkdir -p $dgndir
mkdir -p $gtexdir
mkdir -p $logdir
mkdir -p $mesadir

# run predixcan pipeline
qsub -N sage.predixcan \
    -v predixcan=$predixcan,python2=$python2,Rscript=$Rscript,R_format_predixcan=$R_format_predixcan,outdir=$outdir,dgndir=$dgndir,gtexdir=$gtexdir,mesadir=$mesadir,weights_file_DGN=${weights_file_DGN},weights_file_GTEx_v6p=${weights_file_GTEx_v6p},weights_file_GTEx_v7=${weights_file_GTEx_v7},predixcan_sfx=${predixcan_sfx},predixcan_genodir=${predixcan_genodir},bedfile_pfx=${bedfile_pfx},bedfile_sfx=${bedfile_sfx},R_impute_predixcan_genos=${R_impute_predixcan_genos},weights_file_MESA_AFA=$weights_file_MESA_AFA,weights_file_MESA_AFHI=$weights_file_MESA_AFHI,weights_file_MESA_ALL=$weights_file_MESA_ALL,weights_file_MESA_CAU=$weights_file_MESA_CAU,output_pfx_DGN=$output_pfx_DGN,output_pfx_GTEx_v6p=$output_pfx_GTEx_v6p,output_pfx_GTEx_v7=$output_pfx_GTEx_v7,output_pfx_MESA_AFA=$output_pfx_MESA_AFA,output_pfx_MESA_AFHI=$output_pfx_MESA_AFHI,output_pfx_MESA_ALL=$output_pfx_MESA_ALL,output_pfx_MESA_CAU=$output_pfx_MESA_CAU \
    -o $logdir \
    -e $logdir \
    -l mem_free=$memory_limit \
    -l h_rt=$h_rt \
    -t 1-22 \
    $BASH_run_predixcan
    #${codedir}/qsub_run_predixcan.sh
    #-hold_jid "" \
    #-v "${vars}" \

# postprocess the results
qsub -N sage.predixcan.postprocess \
    -hold_jid sage.predixcan \
    -v Rscript=$Rscript,R_postprocess_predixcan=$R_postprocess_predixcan,output_pfx_DGN=$output_pfx_DGN,output_pfx_GTEx_v6p=$output_pfx_GTEx_v6p,output_pfx_GTEx_v7=$output_pfx_GTEx_v7,output_pfx_MESA_AFA=$output_pfx_MESA_AFA,output_pfx_MESA_AFHI=$output_pfx_MESA_AFHI,output_pfx_MESA_ALL=$output_pfx_MESA_ALL,output_pfx_MESA_CAU=$output_pfx_MESA_CAU \
    -o $logdir \
    -e $logdir \
    -l mem_free=$memory_limit \
    -l h_rt=$h_rt \
    $BASH_postprocess_predixcan

# don't forget to send results to Cesar! cannot easily do this in SGE setting...
