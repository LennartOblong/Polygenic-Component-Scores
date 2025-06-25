#!/bin/bash
# recommended submission script of non-neuroimaging testing on all phenotypes simultaneously. Running all of this is very memory intensive due to the large sample size (>330,000 participants),
# and large number of phenotypes.

# NOTE - output will be generated per split and will need to be merged after running successfully.

#loop for the split and desired components in binary/binarized phenotypes
# 22 splits where chosen due to large phenotype number after binarization (>1800). NOTE the memory and time requirements within the loop below (50 GB of memory per split)
for i in {1..22}; do
  for prefix in IC PC; do
    infile="/path/to/your/input/data/splits/cleandata_${prefix}nonimage_bin_split_${i}.txt"
    outfile="/path/to/your/desired/ouput/${prefix}nonimage_output_bin_split${i}.txt"

    sbatch <<EOF
#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=1
#SBATCH -p genoa
#SBATCH -t 00:20:00
#SBATCH --mem=50GB

module load 2023
module load R/4.3.2-gfbf-2023a

Rscript /path/to/your/scripts/step3b-modelling-non-neuroimaging-data.R ${infile} ${outfile}
EOF

  done
done

# loop for the split and desired components in continuous phenotypes
for i in {1..3}; do
  for prefix in IC PC; do
    infile="/path/to/your/input/data/splits/cleandata_${prefix}nonimage_cont_split_${i}.txt"
    outfile="/path/to/your/desired/ouput//${prefix}nonimage_output_cont_split${i}.txt"

    sbatch <<EOF
#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=1
#SBATCH -p genoa
#SBATCH -t 00:20:00
#SBATCH --mem=50GB

module load 2023
module load R/4.3.2-gfbf-2023a

Rscript /path/to/your/scripts/multivarreg-prspca-nonimage.R.pairwise ${infile} ${outfile}
EOF

  done
done
