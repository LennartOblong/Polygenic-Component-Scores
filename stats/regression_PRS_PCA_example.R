# Supplemental R Code for PRS-PCA
# https://onlinelibrary.wiley.com/doi/full/10.1002/gepi.22339
# A principal component approach to improve association testing with polygenic risk scores

###################
# 09-06-2022
# example code for computing PC1 from PRSice all_scores and conduct regression analysis onto 3 continuous factor outcomes: SES,Hea,Soc
# PRS bases: ADHD, ASD
# Yingjie

###### FUNCTION TO PERFORM PRS-PCA ##########################################
# INPUTS:
# dat = n * (p+1) dataframe of PRSs under different settings with first column as ID
# x = label for PRS (typically named for phenotype (i.e. MDD))
# OUTPUTS:
# list of
#  - data = dataframe with cols (ID , PRS-PCA1 , PRS-PCA2)
#  - r2 = variance explained by each PC of the PRS matrix
#  - loadings = the PC-loadings used to create PRS-PCA1
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
  dat[,paste0(x,"_PC1")] <- scale(pc1)   ## rescales PRS-PCA1
  dat[,paste0(x,"_PC2")] <- scale(pc2)  ## rescales PRS-PCA2
  return(list(data=dat,r2=pca.r2,loadings=pc1.loadings))
}

library(here)

# load data
base_list <- c('ADHD','ASD')

dat_ADHD <- read.table(here('all_score','ADHD_3factor.all_score'), sep = ' ', header = T)[,-1]
dat_SUBAFF <- read.table(here('all_score','ADHD_3factor.all_score'), sep = ' ', header = T)[,-1]

output_ADHD <- prs.pc(dat_ADHD,'ADHD')
PC1_ADHD <- output_ADHD$data[,c(1,10)]

#output_ASD <- prs.pc(dat_ASD,'ASD')
#PC1_ASD <- output_ASD$data[,c(1,10)]

# merge all three phenotypes,and PRS_PC for two bases
pheno_cov <- read.table('H:\\GR Theme groups\\02 PI Group Barbara Franke\\UKB\\CFA_outcome\\pheno_cov_valid.txt', sep = ' ', header = T)
pheno_cov[,c(5,7,8)] <- sapply(pheno_cov[,c(5,7,8)], as.factor)

PRS_seven_PC <- data.frame(IID = pheno_cov$IID,SES = pheno_cov$SES,Hea = pheno_cov$Hea,Soc = pheno_cov$Soc)

PRS_seven_PC <- merge(merge(merge(PRS_seven_PC,PC1_ADHD, by = 'IID'),PC1_ASD, by = 'IID'),pheno_cov[,c(1,5:18)], by = 'IID')

PRS_seven_PC <- PRS_seven_PC[complete.cases(PRS_seven_PC),]

# save df
write.table(PRS_seven_PC, here('df_PCA','PRS_seven_PC.csv'), row.names = F, col.names = T, quote = F)

############
# SES
SES.ADHD.single.PC <- lm(SES ~ ADHD_PC1 + site + age + batch + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = PRS_seven_PC)
SES.ASD.single.PC <- lm(SES ~ ASD_PC1 + site + age + batch + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = PRS_seven_PC)

SES.Null.PC <- lm(SES ~ site + age + batch + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = PRS_seven_PC)

SES.Full.R2.PC <- c(summary(SES.ADHD.single.PC)$r.squared,
                 summary(SES.ASD.single.PC)$r.squared)

SES.Null.R2.PC <- rep(summary(SES.Null.PC)$r.squared,2)
SES.R2.PC.df <- data.frame(SES.Full.R2.PC,SES.Null.R2.PC)

SES.R2.PC.df$PRS.R2 <- 1 - (1 - SES.Full.R2.PC) / (1 - SES.Null.R2.PC) #if you are using PRS from the best threshold, you should be getting the exact same R2 from the PRSice .summary file

###########
# Hea
Hea.ADHD.single.PC <- lm(Hea ~ ADHD_PC1 + site + age + batch + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = PRS_seven_PC)
Hea.ASD.single.PC <- lm(Hea ~ ASD_PC1 + site + age + batch + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = PRS_seven_PC)

Hea.Null.PC <- lm(Hea ~ site + age + batch + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = PRS_seven_PC)

Hea.Full.R2.PC <- c(summary(Hea.ADHD.single.PC)$r.squared,
                    summary(Hea.ASD.single.PC)$r.squared)

Hea.Null.R2.PC <- rep(summary(Hea.Null.PC)$r.squared,7)
Hea.R2.PC.df <- data.frame(Hea.Full.R2.PC,Hea.Null.R2.PC)

Hea.R2.PC.df$PRS.R2 <- 1 - (1 - Hea.Full.R2.PC) / (1 - Hea.Null.R2.PC)

#############
# Soc
Soc.ADHD.single.PC <- lm(Soc ~ ADHD_PC1 + site + age + batch + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = PRS_seven_PC)
Soc.ASD.single.PC <- lm(Soc ~ ASD_PC1 + site + age + batch + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = PRS_seven_PC)

Soc.Null.PC <- lm(Soc ~ site + age + batch + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = PRS_seven_PC)

Soc.Full.R2.PC <- c(summary(Soc.ADHD.single.PC)$r.squared,
                    summary(Soc.ASD.single.PC)$r.squared)

Soc.Null.R2.PC <- rep(summary(Soc.Null.PC)$r.squared,7)
Soc.R2.PC.df <- data.frame(Soc.Full.R2.PC,Soc.Null.R2.PC)

Soc.R2.PC.df$PRS.R2 <- 1 - (1 - Soc.Full.R2.PC) / (1 - Soc.Null.R2.PC)

############ r2 matrix 3*7
PRS.R2 <- data.frame(SES.R2 = SES.R2.PC.df$PRS.R2,
                     Hea.R2 = Hea.R2.PC.df$PRS.R2,
                     Soc.R2 = Soc.R2.PC.df$PRS.R2)
rownames(PRS.R2) <- base_list

write.table(PRS.R2, here('df_PCA','PRS.R2.csv'), row.names = T, col.names = T, quote = F, sep=',')
R2_PCA <- read.table(here('df_PCA','PRS.R2.csv'), sep = ',', header = T)

