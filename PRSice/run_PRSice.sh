### This script is used to run PRSice on multiple components one after the other. To save time this script can be started multiple times with different IC/PC name files ###
### Lennart Oblong 01.03.2022 ###

#!/bin/bash

#load modules

module load 2021
module load R/4.1.0-foss-2021a
module load PRSice/2.3.3-GCCcore-10.3.0
#set directories. individual outdirs for PCs and ICs for ADHD_EVER and SUBAFF_EVER respectively

indir=/projects/0/einf2700/oblongl/polygen_comp-NeuroIMAGE/polygen_comp
target=/projects/0/einf2700/oblongl/NeuroIMAGE_genetics_imputedcopy/ricopili_imp_SELECTIVE_UNTAR_4tim/cobg_dir_genome_wide
covdir=/projects/0/einf2700/oblongl/NeuroIMAGE_genetics_imputedcopy

#outdir1=/projects/0/einf2700/oblongl/polygen_comp/PRSice-output_NeuroIMAGE-ADHD-PCs-ricopili-bgn
#outdir2=/projects/0/einf2700/oblongl/polygen_comp/PRSice-output_NeuroIMAGE-subaffADHD-PCs-ricopili-bgn
#outdir3=/projects/0/einf2700/oblongl/polygen_comp/PRSice-output_NeuroIMAGE-ADHD-ICs-ricopili-bgn
#outdir4=/projects/0/einf2700/oblongl/polygen_comp/PRSice-output_NeuroIMAGE-subaffADHD-ICs-ricopili-bgn
IC_names=/projects/0/einf2700/oblongl/polygen_comp-NeuroIMAGE/polygen_comp/IC25_names.txt
#IC_names=/projects/0/einf2700/oblongl/polygen_comp-NeuroIMAGE/polygen_comp/ICall_names.txt
#PC_names=/projects/0/einf2700/oblongl/polygen_comp-NeuroIMAGE/polygen_comp/PCall_names.txt
PC_names=/projects/0/einf2700/oblongl/polygen_comp-NeuroIMAGE/polygen_comp/PC10_names.txt

### If the desired component directories don't exist, they are created heare. Note: create outdirs beforehand, otherwise this will not work ###


Rscript /projects/0/einf2700/oblongl/scripts/PRSice.R --dir . \
                               --prsice /projects/0/einf2700/oblongl/scripts/PRSice_linux \
                               --base /projects/0/einf2700/oblongl/polygen_comp/ICA_ic_2.assoc \
                               --target ${target}/adh_amd1_eur_es-postpc2.hg19.ch.fl.bgn \
                               --beta T \
                               --binary-target T \
                               --no-clump  \
                               --lower 0.01 \
                               --quantile 100 \
                               --quant-break 1,5,10,20,40,60,80,90,95,99,100 \
                               --quant-ref 60 \
                               --out  /projects/0/einf2700/oblongl/polygen_comp-NeuroIMAGE_PRSice-output_single_cov/PRSice_ICA_ic_2_SUBAFF_EVER \
                               --pheno ${target}/traits/Mphen_binary_Cswap.csv \
                               --pheno-col SUBAFF_EVER \
                               --cov ${covdir}/Mcovall_4prsice_FID_tab.txt \
                               --cov-col @C[1-4] \
#                               --perm 10000 \
                               --thread max \
                               --upper 1

