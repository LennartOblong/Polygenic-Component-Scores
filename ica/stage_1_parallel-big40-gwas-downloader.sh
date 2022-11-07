#Sourena Soheili-Nezhad, 2019
#!/bin/bash
#SBATCH -N 1
#SBATCH -p normal
#SBATCH -t 00:20:00

#retreives UKB summary statistics from https://open.win.ox.ac.uk/ukbiobank/big40
#and converts to nifti-1

#put the current working directory below 
#cd /home/ssoheili/genetic-data/genica/big40

#modify the FSLDIR
FSLDIR=/home/ssoheili/software/fsl
. ${FSLDIR}/etc/fslconf/fsl.sh
PATH=${FSLDIR}/bin:${PATH}
export FSLDIR PATH

# declare parallel, modify the path for your cluster
parallel=/lustre2/0/multifac/bin/parallel

#The first argument (below) is a text file with a list of all phenotype codes (IDPs) in rows to be downloaded and processed.
#Optimally this script is submitted multiple times, each time with a different list of e.g. 24 IDPs,
#so the process is parallelized as much as possible
#Alternatively, just put all IDP codes (thousands) in a single text file and pass as an argument. But this will only be run in a single node

pheno_list_24=$1
# wrap steps into a function so all cores are kept occupied
function ukbb_downloader {
        pheno=$1

        wget -q https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/${pheno}.txt.gz -O ${pheno}-33k.txt.gz
        wget -q https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats/${pheno}.txt.gz -O ${pheno}-22k.txt.gz
        wget -q https://open.win.ox.ac.uk/ukbiobank/big40/release2/repro/${pheno}.txt.gz -O ${pheno}-11k.txt.gz

        gunzip ${pheno}-33k.txt.gz
        gunzip ${pheno}-22k.txt.gz
        gunzip ${pheno}-11k.txt.gz

        cat ${pheno}-33k.txt |awk 'NR>1 {print $6}' > temp_${pheno}_beta_33k.txt
        cat ${pheno}-33k.txt |awk 'NR>1 {print $8}' > temp_${pheno}_logp_33k.txt

        cat ${pheno}-22k.txt |awk 'NR>1 {print $6}' > temp_${pheno}_beta_22k.txt
#        cat ${pheno}-22k.txt |awk 'NR>1 {print $8}' > temp_${pheno}_logp_22k.txt

        cat ${pheno}-11k.txt |awk 'NR>1 {print $6}' > temp_${pheno}_beta_11k.txt
#        cat ${pheno}-11k.txt |awk 'NR>1 {print $8}' > temp_${pheno}_logp_11k.txt

#Integer factorization: 17103079 SNPs total = 10223*239*7
#This is a fast and arbitrary workaround to encode the 1D SNP vector in a 3D volumetric file, 
#since NIFTI-1 doesn't support any dimension with more than 32k voxels
#CIFTI and NIFTI-2 can handle that, but MELODIC predates them (and I'm old-fashioned)

	fslascii2img  temp_${pheno}_beta_33k.txt 10223 239 7 1 1 1 1 1  nifti/big_33k_beta_${pheno}
        fslascii2img  temp_${pheno}_logp_33k.txt 10223 239 7 1 1 1 1 1  nifti/big_33k_logp_${pheno}

        fslascii2img  temp_${pheno}_beta_22k.txt 10223 239 7 1 1 1 1 1  nifti/big_22k_beta_${pheno}
#        fslascii2img  temp_${pheno}_logp_22k.txt 10223 239 7 1 1 1 1 1  nifti/big_22k_logp_${pheno}

        fslascii2img  temp_${pheno}_beta_11k.txt 10223 239 7 1 1 1 1 1  nifti/big_11k_beta_${pheno}
#        fslascii2img  temp_${pheno}_logp_11k.txt 10223 239 7 1 1 1 1 1  nifti/big_11k_logp_${pheno}

        rm temp_${pheno}_beta_33k.txt temp_${pheno}_logp_33k.txt ${pheno}-33k.txt
        rm temp_${pheno}_beta_22k.txt ${pheno}-22k.txt #temp_${pheno}_logp_22k.txt
        rm temp_${pheno}_beta_11k.txt ${pheno}-11k.txt #temp_${pheno}_logp_11k.txt
        }

# export the function to the environment so that parallel can access it outside the bash script
export -f ukbb_downloader


# run the function in parallel. Change number of parallel jobs with --jobs N depending on how many cores/threads the CPU has/can handle
# "-a file " takes a text file as input and every line is a argument
# parallel will keep a number of jobs (--jobs N) running and start new ones when the older ones are finished
cat $pheno_list_24 | $parallel --jobs 12 ukbb_downloader {1}

