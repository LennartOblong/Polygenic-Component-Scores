#!/bin/bash
#$ -cwd

#after downloading and converting betas/pvals into 3D NIFTI volumes
#This file merges them into a 4D nifti
#The list of all IDPs are read from the all-pheno.txt file
#This step is memory-hungry with thousands of IDPs, so need to be run on the big memory node

source /home/sousoh/.bashrc


fslmerge -t all-33k-beta $(cat all-pheno.txt |while read s;do echo -en "33k/big_33k_beta_${s} ";done)
fslmerge -t all-33k-logp $(cat all-pheno.txt |while read s;do echo -en "33k/big_33k_logp_${s} ";done)


fslmerge -t all-22k-beta $(cat all-pheno.txt |while read s;do echo -en "22k/big_22k_beta_${s} ";done)
fslmerge -t all-11k-beta $(cat all-pheno.txt |while read s;do echo -en "11k/big_11k_beta_${s} ";done)

#This is the SNP-wise min-pvalue across all IDPs, so a 3D volume
#Can be used e.g. for clumping, masking, etc.
fslmaths all-33k-logp -Tmax all-33k-logp-max
