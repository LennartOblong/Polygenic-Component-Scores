#!/bin/bash

module load 2021
module load R/4.1.0-foss-2021a
module load PRSice/2.3.3-GCCcore-10.3.0

### set directories for PRSice ###

### set directories for PRSice base and target data. Create PC and IC lists as desired to cat over the list quickly ###
### Always check these!!! ###

indir=/projects/0/einf2700/oblongl/polygen_comp
target=/projects/0/einf2700/oblongl/NeuroIMAGE_genetics_imputedcopy/ricopili_imp_SELECTIVE_UNTAR_4tim/cobg_dir_genome_wide
covdir=/projects/0/einf2700/oblongl/NeuroIMAGE_genetics_imputedcopy
outdir=/projects/0/einf2700/oblongl/test_dir/PRS-PCA_test

comp_names=/projects/0/einf2700/oblongl/polygen_comp/IC_PC_names.txt

# PRSice main loop

for names in ${comp_names} ; do

        for pheno_col in SUBAFF ADHD ; do

                for comp in $(cat ${names}) ; do

                        echo "variables used now:"
                        echo "${names} ${pheno_col} ${comp}"

                       Rscript /projects/0/einf2700/oblongl/scripts/PRSice.R --dir . \
                               --prsice /projects/0/einf2700/oblongl/scripts/PRSice_linux \
                               --base ${indir}/${comp}.assoc \
                               --target ${target}/adh_amd1_eur_es-postpc2.hg19.ch.fl.bgn \
                               --beta T \
                               --binary-target T \
                               --lower 0.01 \
			       --all-score \
                               --no-clump  \
			       --quantile 100 \
			       --quant-break 1,5,10,20,40,60,80,90,95,99,100 \
			       --quant-ref 60 \
                               --out  ${outdir}_${comp}_${pheno_col} \
                               --pheno ${target}/traits/Mphen_binary_Cswap.csv \
                               --pheno-col ${pheno_col}_EVER \
                               --cov ${covdir}/Mcovall_4prsice_FID_tab.txt \
                               --cov-col @C[1-4] \
#                               --perm 10000 \
                               --thread max \
                               --upper 1

                done

        done

done



