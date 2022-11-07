#!/bin/bash

### Script to call the PRS-PCA function over a job_node to ensure that enough working memory is present to finish calculation ###

#SBATCH -N 1
#SBATCH -p thin
#SBATCH -t 72:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128

### note that the more specifications a jobscript has, the fewer the number of potential nodes able to compute the job. Job may fail because of that so think carefully about batch details ###

### load  modules necessary to perform submitted applications ###

module load 2021
module load R/4.1.0-foss-2021a

cd $TMPDIR

mkdir results_PRS-PCA

### The component for analysis is hard-coded in the script called below - Always check this ###

Rscript /projects/0/einf2700/oblongl/scripts/PRS-PCA.R

cp results_PRS-PCA/* /projects/0/einf2700/oblongl/PRS-PCA_output/
