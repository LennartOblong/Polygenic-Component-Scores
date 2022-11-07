###### FUNCTION TO PERFORM PRS-PCA ##########################################
# INPUTS:
# dat = n x (p+1) dataframe of PRSs under different settings with first column as ID
# x = label for PRS (typically named for phenotype (i.e. MDD))
# OUTPUTS:
# list of 
#  - data = dataframe with cols (ID , PRS-PCA1 , PRS-PCA2)
#  - r2 = variance explained by each PC of the PRS matrix
#  - loadings = the PC-loadings used to create PRS-PCA1

### the TMPDIR out-directory is used when the script is passed as a job to the cluster to increase computational efficiency ###

#outdir <- ("/projects/0/einf2700/oblongl/PRS-PCA_output")
outdir <- ("/$TMPDIR/results_PRS-PCA")

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

dat = read.table("/projects/0/einf2700/oblongl/polygen_comp-NeuroIMAGE_PRSice-output_all-score-copy/test.all.score")
x = "SUBAFF_EVER"

prs.pc(dat, x)

PRS-PCA_SUB-IC2 <- as.data.frame(data)
PRS-PCA_SUB-IC2_R2 <- as.data.frame(pca.r2)
PRS-PCA_SUB-IC2_pc1loadings <- as.data.frame(pc1.loadings)

PRS-PCA_SUB-IC2_file <- paste0(outdir, "/" ,"PRS-PCA_SUB-IC2.csv")
write.csv(PRS-PCA_SUB-IC2, PRS-PCA_SUB-IC2_file, row.names = T)

PRS-PCA_SUB-IC2_R2_file <- paste0(outdir, "/", "PRS-PCA_SUB-IC2_R2.csv")
write.csv(PRS-PCA_SUB-IC2_R2, PRS-PCA_SUB-IC2_R2_file, row.names = T)

PRS-PCA_SUB-IC2_pc1loadings_file <- paste0(outdir, "/", "PRS-PCA_SUB-IC2_pc1loadings.csv")
write.csv(PRS-PCA_SUB-IC2_pc1loadings, PRS-PCA_SUB-IC2_pc1loadings_file, row.names = T)



#save(pc1, file = "PC1_SUBAFF_IC2")
#save(pc2, file = "PC2_SUBAFF_IC2")


#dat2 = list.files(path="/projects/0/einf2700/oblongl/polygen_comp-NeuroIMAGE_PRSice-output_all-score-copy/", pattern="*_ICA_*_SUBAFF_EVER.all.score", full.names=TRUE, recursive=FALSE)
#x2 = "SUBAFF_EVER-ICA"


#dat3 = list.files(path="/projects/0/einf2700/oblongl/polygen_comp-NeuroIMAGE_PRSice-output_all-score-copy/", pattern="*_PCA_*_ADHD_EVER.all.score", full.names=TRUE, recursive=FALSE)
#x3 = "ADHD_EVER-PCA"


#dat4 = list.files(path="/projects/0/einf2700/oblongl/polygen_comp-NeuroIMAGE_PRSice-output_all-score-copy/", pattern="*_PCA_*_SUBAFF_EVER.all.score", full.names=TRUE, recursive=FALSE)
#x4 = "SUBAFF_EVER-PCA"

