#!/bin/bash

### Script to run PRSice2 on a cluster node for ICs, PCs, and two phenotypic variants ###
### Lennart Oblong 31.03.22 ###

#SBATCH -N 1
#SBATCH -p thin
#SBATCH -t 118:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128

### note that the more specifications a jobscript has, the fewer the number of potential nodes able to compute the job. Job may fail because of that so think carefully about batch details ###

### load  modules necessary to perform submitted applications ###

module load 2021
module load R/4.1.0-foss-2021a
module load PRSice/2.3.3-GCCcore-10.3.0

### set directories for PRSice base and target data. Create PC and IC lists as desired to cat over the list quickly ###
### Always check these!!! ###

indir=/projects/0/einf2700/oblongl/polygen_comp
target=/projects/0/einf2700/oblongl/NeuroIMAGE_genetics_imputedcopy/ricopili_imp_SELECTIVE_UNTAR_4tim/cobg_dir_genome_wide
covdir=/projects/0/einf2700/oblongl/NeuroIMAGE_genetics_imputedcopy

comp_names=/projects/0/einf2700/oblongl/polygen_comp/IC_PC_names.txt
#comp_names=/projects/0/einf2700/oblongl/polygen_comp/sig_IC_PC_names.txt

### move to the scratch directory for more efficient computation ###

cd $TMPDIR

mkdir results

### compute the PRSice for desired components and two phenotypes per component. The below loop runs the analysis on ADHD phenotypes in the NeuroIMAGE cohort ###
### Always check indir, targetdir and comp_lists before running to avoid overwriting data ###
### The below loop can run in two ways, depending on the commenting out of flags. If all-score is enabled, all sbj pvals will be extracted to one file across all pval thresholds. ###
### This is needed for PRS-PCA, and makes permutation testing unenecessary. If permutation testing is desired, uncomment that flag but do not forget to comment out all-score, as the file will be huge. ###

        for pheno_col in SUBAFF ADHD ; do

                for comp in $(cat ${comp_names}) ; do

                        echo "variables used now:"
                        echo "${pheno_col} ${comp}"


		Rscript /projects/0/einf2700/oblongl/scripts/PRSice.R --dir . \
                               --prsice /projects/0/einf2700/oblongl/scripts/PRSice_linux \
                               --base ${indir}/${comp}.assoc \
                               --target ${target}/adh_amd1_eur_es-postpc2.hg19.ch.fl.bgn \
                               --beta T \
                               --binary-target T \
                               --no-clump  \
                               --lower 0.01 \
                               --quantile 100 \
                               --quant-break 1,5,10,20,40,60,80,90,95,99,100 \
                               --quant-ref 60 \
                               --out  results/PRSice_${pheno_col}_EVER_${comp} \
                               --pheno ${target}/traits/Mphen_binary_Cswap.csv \
                               --pheno-col ${pheno_col}_EVER \
                               --cov ${covdir}/Mcovall_4prsice_FID_tab.txt \
                               --cov-col @C[1-4] \
                               --perm 10000 \
                               --thread max \
                               --upper 1

                done

        done

### Copy the processed data from scratch back to desired out directory ###
### Always check this !!! ###

#cp results/* /projects/0/einf2700/oblongl/polygen_comp-UKB_PRSice-output/´
#cp results/* /projects/0/einf2700/oblongl/polygen_comp-NeuroIMAGE_PRSice-output/
cp results/* /projects/0/einf2700/oblongl/polygen_comp-NeuroIMAGE_PRSice-output_cov_perm/




