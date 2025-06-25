# Description
Running polygenic score analysis on GenomICA output - code repository for reproducibility. 
Contains the code to reproduce the results of our work to explore GenomICA output out-of-sample, published as a preprint on MedRXiv (https://www.medrxiv.org/content/10.1101/2025.06.06.25329112v1)
This work is a follow up to our previous work using the FSL MELODIC algorithm v3.15 (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/MELODIC) on genetic data.
Please cite our previous work when implementing genomICA in your own work. (https://onlinelibrary.wiley.com/doi/10.1111/gbb.12876)

A new and improved iteration of genomic components with visualizations and associated traits is available under http://genomica.info/ ! Should you wish to gain access to components for your own research please contact lennart.oblong@donders.ru.nl.
The output of the current analysis is also presented on this website, showing the vast trait associations we uncovered using component-based polygenic component scores.

This repository contains the step-by-step analysis outlined in our publication (https://www.medrxiv.org/content/10.1101/2025.06.06.25329112v1).
The code is divided into 3 stages. Following the stages from 1 through 3 demonstrates how we generated the polygenic component scores from genomic indepenedent and principal components.
It also includes the UK Biobank data QC that we performed and the binarization of categorical and ordinal phenotypes. Lastly, general linear and logistic regression models are run. 
To accurately reproduce our results read the comments in the scripts carefully and pay attention to the "important notes" section below.

    STAGE 1
      Calculate polygenic component scores per component using the PRS-PCA method
    STAGE 2
      Data wrangling and data QC on large UK Biobank data
    STAGE 3a & 3b
      Run linear and logistic regression models to compute the phenotypic variance explained by components. Final files contain the contributions of individual components

# Important Notes
This repository does not contain plug-and-play code. The code provided here is meant for checking the reproducibility of our analysis, and shows the implementation of polygenic score analyses on the output of GenomICA components.
All the analyses were run on Snellius (https://www.surf.nl/en/dutch-national-supercomputer-snellius), using the Slurm workload management system (https://slurm.schedmd.com/documentation.html).
Some of the steps are memory hungry and require a high-performance computing environment.
