#!/bin/bash
# Lennart Oblong 23.05.2022 #

module load 2021
module load R/4.1.0-foss-2021a
module load PRSice/2.3.3-GCCcore-10.3.0

# INPUTS:
# dat = n x (p+1) dataframe of PRSs under different settings with first column as ID
# x = label for PRS (typically named for phenotype (i.e. MDD))
# OUTPUTS:
# list of 
#  - data = dataframe with cols (ID , PRS-PCA1 , PRS-PCA2)
#  - r2 = variance explained by each PC of the PRS matrix
#  - loadings = the PC-loadings used to create PRS-PCA1

Rscript

dat = sink("/projects/0/einf2700/oblongl/polygen_comp-NeuroIMAGE_PRSice-output_all-score/PRSice_ICA_ic_25_SUBAFF_EVER.all.score")
x = SUBAFF_EVER
prs.pc <- function(dat,x){
  xo <- scale(as.matrix(dat[,-1]))  ## scale cols of matrix of only PRSs (remove ID)
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
  dat[,paste0(x,".prs.pc")] <- scale(pc1)   ## rescales PRS-PCA1 
  dat[,paste0(x,".prs.pc2")] <- scale(pc2)  ## rescales PRS-PCA2
  return(list(data=df,r2=pca.r2,loadings=pc1.loadings))
}

