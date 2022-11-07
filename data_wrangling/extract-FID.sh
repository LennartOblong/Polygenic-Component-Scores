#!/bin/bash

# Lennart Oblong - script to extract middle number out of subject ID. 

ID=/projects/0/einf2700/oblongl/PRSice_stats/Subject.txt
outdir=/projects/0/einf2700/oblongl/PRSice_stats

 for FID in $(cat ${ID}) ; do

	echo "$FID"| cut -d '-' -f 2 

 done  > $outdir/test.txt

