#!/bin/bash
#$ -cwd

source /home/sousoh/.bashrc
#Check the absolute path or config FSL based on cluster	
melodic=/data/clusterfs/admin/local/apps/fsl/bin/melodic


#for masking based on top 5 percentile in variance, we ditched this approach in BIG-40 data
#fslmaths all-33k-std -thr $(fslstats all-33k-std -p 95) -bin mask-percent-95

for i in 11k 22k 33k
do
for dim in 7 15 20 75 100 150 200 
do
#old version: variance-based masking: fslmaths all-${i}-beta -Tstd all-${i}-std 
$melodic -i all-${i}-beta --no_mm --mask=/data/clusterfs/lag/projects/lg-genica/big40/nifti/clump/mask-clump-all-three --update_mask --Oall -o ica-clumped-no-dim-${dim}-${i} --keep_meanvol --rescale_nht --vn  -dim $dim
done
done
