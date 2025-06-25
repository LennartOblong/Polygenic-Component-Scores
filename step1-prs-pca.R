## Script containing the code to generate polygenic component scores from genomICA components.

##########################################################################################################
## Prior to creating PGS-PCs the initial PGS need to be calculated. We used the PRSice2 toolbox by Choi et al. (https://choishingwan.github.io/PRSice/ | https://doi.org/10.1093/gigascience/giz082.) 
## with the following options.

#componentfile=$1
#
#Rscript /path/to/your/bin/PRSice.R \
#    --prsice /path/to/your/PRSice_linux \
#    --base "$componentfile" \
#        --snp SNP --chr CHR --bp BP --A1 A1 --A2 A2 --stat BETA --pvalue P \ #BETA represents the component loadings
#    --target /path/to/target/data/ukb_imp_QCed_plink2_data \ #plink2 data format
#    --binary-target F \
#    --no-regress \	#no regression at this stage
#    --fastscore \	
#    --all-score \	#output PRS for all p-thresholds
#    --clump-r2 0.1 \	#clump at 0.1
#    --extract /path/to/PRSice.valid \ #if not obtained before run PRSice blind to have this created. Contains all matching SNPs across base and target
#    --out "/path/to/your/desired/output/${componentfile##*/}"

## Default options for p-value thresholds where used (0.001, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1)
##########################################################################################################

########### NOTE ###########
# PGS where calculated in a non-overlapping sample of UK Biobank not included in the inital decomposition. From the resulting data only those subjects with neuroimaging phenotypes where 
# extracted and PRS-PCA was calculated on that sample seperately to derive PCS for association testing.
# The full dataset of participants was used to calculate PRS-PCA to derive PCS for association testing in non-neuroimaging phenotypes.
############################


# The below supplemental R Code for PRS-PCA published in Coombes et al.
# https://onlinelibrary.wiley.com/doi/full/10.1002/gepi.22339
# "A principal component approach to improve association testing with polygenic risk scores"

# Altered to generate polygenic component scores of genomICA derived principal and independent components of brain structure and function.

## inputs needed for prs-pca function:
## dat = n * (p+1) dataframe of PRSs under different settings with the first two columns as IDs (they are ignored during PC computation)
## x = label for PRS (typically named for phenotype (i.e. MDD))
## OUTPUTS:
## list of
##  - data = dataframe with cols (ID , PRS-PCA1 , PRS-PCA2)
##  - r2 = variance explained by each PC of the PRS matrix
##  - loadings = the PC-loadings used to create PRS-PCA1

####### function
prs.pc <- function(dat,x){
  xo <- scale(as.matrix(dat[,-c(1,2)]))  ## scale cols of matrix of only PRSs (remove ID columns)
  g <- prcomp(xo)   ## perform principal components
  pca.r2 <- g$sdev^2/sum(g$sdev^2)    ## calculate variance explained by each PC
  pc1.loadings <- g$rotation[,1];     ## loadings for PC1
  pc2.loadings <- g$rotation[,2]      ## loadings for PC2
  ## flip direction of PCs to keep direction of association
  ## (sign of loadings for PC1 is arbitrary so we want to keep same direction)
  if (mean(pc1.loadings>0)==0){
    pc1.loadings <- pc1.loadings*(-1)
    pc2.loadings <- pc2.loadings*(-1)
  }
  ## calculate PRS-PCA (outputs PC1 and PC2 even though PC1 sufficient)
  pc1 <- xo %*% pc1.loadings
  pc2 <- xo %*% pc2.loadings
  dat[,paste0(x,"_PC1")] <- scale(pc1)   ## rescales PRS-PCA1
  dat[,paste0(x,"_PC2")] <- scale(pc2)  ## rescales PRS-PCA2
  return(list(data=dat,r2=pca.r2,loadings=pc1.loadings))
}

####### directories
indir <- "path/to/ukb/data" #phenotype files, covariate files, keys
outdir <- "path/to/desired/output"
img_pgs_dir <- "path/to/imaging/component/PGS" #component PGS at desired p-val thresholds for neuroimaging phenotypes
nonimg_pgs_dir <- "path/to/nonimaging/component/PGS" #component PGS at desired p-val thresholds for non-neuroimaging phenotypes

####### library
library(data.table)
library(dplyr)
library(tidyverse)
library(stringr)
library(rms)
library(DescTools)
library(reshape2)
library(foreach)
library(doParallel)

#######calculate prs-PCs for all components that have neuroimaging data - store in list object which is used in later scripts. Individual files are also saved for individual examination.
#######adjust filenames and p-thresholds to fit your data naming convention and p-thresholds applied
cat("Commencing PRS-PCA calculation in participants with neuroimaging phenotypes ...", "\n")
output_list_img <- list()
for (i in 1:16) { #16 components - adjust for number of components
          for (j in c("oIC", "pca")) { #calculate prs-pca for ICA and PCA outputs
        filename <- paste0("data_", j, "_", i, ".txt.all_score")
        data <- fread(file.path(img_pgs_dir, filename))
        prs_pca_output <- prs.pc(dat = data, x = "PRS")
        cat("PRS PC calculation completed for", filename, "\n")
        output_df <- prs_pca_output$dataPt_0.001 Pt_0.05 Pt_0.1 Pt_0.2 Pt_0.3 Pt_0.4 Pt_0.5 Pt_1
        output_filename <- paste0("prs-pca_output_all_img", i, "_", j, ".txt")
        #store all output in a list
        output_list_img[[paste0("output_", i, "_", j)]] <- prs_pca_output   ##write output to list. to access the outputs call the elements of the list e.g. output_list[["output_i_oIC"]]$data
        cat ("Stored in output_list_img as element", paste0("output_", i, "_", j), "\n")
        write.table(output_df, file.path(outdir, output_filename), sep="\t", row.names = FALSE) ##seperate file is saved in data_dir
        }
}
cat("Done", "\n")

#examine R2 for prs-pcs from genomic PCs and ICs to confirm salience of first component
cat("Retrieving PRS-PC variance explained ...", "\n")
output_matrix_oIC <- matrix(nrow = 16, ncol = 8)
output_matrix_pca <- matrix(nrow = 16, ncol = 8)

for (i in 1:16) {
    r2_i_oIC <- output_list_img[[paste0("output_", i, "_oIC")]][["r2"]]
    r2_i_pca <- output_list_img[[paste0("output_", i, "_pca")]][["r2"]]
    output_matrix_oIC[i, ] <- unlist(r2_i_oIC)
    output_matrix_pca[i, ] <- unlist(r2_i_pca)
}

output_df_oIC <- as.data.frame(output_matrix_oIC)
output_df_pca <- as.data.frame(output_matrix_pca)

row.names(output_df_oIC) <- paste0("genomic-IC", 1:16)
colnames(output_df_oIC) <- paste0("PRS-PC", 1:8)
row.names(output_df_pca) <- paste0("genomic-PC", 1:16)
colnames(output_df_pca) <- paste0("PRS-PC", 1:8)

write.csv(output_df_oIC, file = paste0(outdir, "/R2_oIC_pcs_img.csv"), row.names = TRUE)
write.csv(output_df_pca, file = paste0(outdir, "/R2_pca_pcs_img.csv"), row.names = TRUE)
cat("Done", "\n")

#######calculate prs-pca for all components that have non-neuroimaging data - store in list object which is used in later scripts. Individual files are also saved for individual examination
cat("Commencing PRS-PCA calculation in participants with non-neuroimaging phenotypes ...", "\n")
output_list_nonimg <- list()
for (i in 1:16) {
  for (j in c("IC", "PC")) {
    filename <- paste0(j,"_", i,".txt.all_score")  ##get PRSs calculated seperately in plink1.9
    data <- fread(file.path(prs_dir, filename))
    prs_pca_output <- prs.pc(dat = data, x = "PRS") ##run prs-pca
    cat("PRS PC calculation completed for", filename, "\n")
    output_df <- prs_pca_output$data ##create output dataframe and filename
    output_filename <- paste0("prs-pca_out_all_nonimg", i, "_", j, ".txt")
    #store all output in a list for easy access in this script
    output_list_nonimg[[paste0("output_", i, "_", j)]] <- prs_pca_output   ##to access the outputs call the elements of the list e.g. output_list[["output_i_j"]]$data
    cat ("Stored in output_list as element", paste0("output_", i, "_", j), "\n")
    write.table(output_df, file.path(outdir, output_filename), sep="\t", row.names = FALSE) ##write each components prs-pca to file for future reference
  }
}

cat("Done", "\n")

cat("Retrieving PRS-PC variance explained ...", "\n")
output_matrix_oIC <- matrix(nrow = 16, ncol = 8)
output_matrix_pca <- matrix(nrow = 16, ncol = 8)

for (i in 1:16) {
  r2_i_oIC <- output_list_nonimg[[paste0("output_", i, "_IC")]][["r2"]]
  r2_i_pca <- output_list_nonimg[[paste0("output_", i, "_PC")]][["r2"]]
  output_matrix_oIC[i, ] <- unlist(r2_i_oIC)
  output_matrix_pca[i, ] <- unlist(r2_i_pca)
}

output_df_oIC <- as.data.frame(output_matrix_oIC)
output_df_pca <- as.data.frame(output_matrix_pca)
row.names(output_df_oIC) <- paste0("genomic-IC", 1:16)
colnames(output_df_oIC) <- paste0("PRS-PC", 1:8)
row.names(output_df_pca) <- paste0("genomic-PC", 1:16)
colnames(output_df_pca) <- paste0("PRS-PC", 1:8)

write.csv(output_df_oIC, file = paste0(outdir, "/R2_oIC_pcs_nonimg.csv"), row.names = TRUE)
write.csv(output_df_pca, file = paste0(outdir, "/R2_pca_pcs_nonimg.csv"), row.names = TRUE)
cat("Done, step 1 completed.", "\n")
