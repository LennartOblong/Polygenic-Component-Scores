### Lennart Oblong 09.03.2022 ###

#!/bin/bash

module load 2021

### set up directories to pull input files from and output dir --- always check these!###

comp_names=/projects/0/einf2700/oblongl/polygen_comp/IC_PC_names.txt
#indir1=/projects/0/einf2700/sourena/33k-no-fold-oics/
indir1=/projects/0/einf2700/sourena/33k-1000GP-oics/ics
indir2=/projects/0/einf2700/sourena/33k-1000GP-oics/ics/ics-pval
outdir=/projects/0/einf2700/oblongl/polygen_comp
#outdir=/projects/0/einf2700/oblongl/polygen_comp_toy

### create input files for PRSice (mirorring GWAS summary stat files (.assoc)) to run. Data is pulled from previous work by Sourena ###
### this will require information from multiple files to be concatenated into one .assoc file ###
# include information on character, SNP ID, base pair,  effective alleles (A1&A2), SNP loadings within comp, and associated pvals. Add headers accordingly#

# get the post-clumping SNPs of the components and isolate them into a new file. Delete column 6 in the process (would mess up risk score analysis if left in).

#awk '$6==1' /projects/0/einf2700/oblongl/scripts/variants-17103079-and-triple-mask.txt | cut -d" " -f1-5 > /projects/0/einf2700/oblongl/polygen_comp/variants_165364-post-clump.txt

for comp in $(cat ${comp_names}) ; do

        paste \
                /projects/0/einf2700/sourena/33k-1000GP-oics/almanac/variants-17103079-1000gp-base.txt \
                ${indir1}/${comp} \
                ${indir2}/${comp}_pval | awk 'BEGIN { print "SNP CHR BP A1 A2 BETA P" } FNR==NR {print}' \
        > ${outdir}/${comp}.assoc

done

