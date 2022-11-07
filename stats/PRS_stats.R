# Lennart Oblong 01.06.22 #
#################Script to perform statistical analysis on PRS output - concatenate emp-pvals and then compute adjusted pvals with FDR. Check association with multiple metrics ##############


#p.adjust.methods
# c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
#   "fdr", "none")

###############get packages

require(readr)  
require(dplyr)  
require(tidyr)  
require(purrr)
require(nlme) 
require(lme4)

#############define files and outdir

outdir <- ("/projects/0/einf2700/oblongl/PRSice_stats/")
comp_names <- read.csv("/projects/0/einf2700/oblongl/polygen_comp/IC_PC_names.txt", header=FALSE)
filenames_all <- list.files(path ='/projects/0/einf2700/oblongl/polygen_comp-NeuroIMAGE_PRSice-output_cov_perm', pattern = '.summary', full.names = TRUE)
filenames1 <- list.files(path ='/projects/0/einf2700/oblongl/polygen_comp-NeuroIMAGE_PRSice-output_cov_perm', pattern = 'PRSice_ADHD_EVER_\\s*(.*?)\\s*.summary', full.names = TRUE)
filenames2 <- list.files(path ='/projects/0/einf2700/oblongl/polygen_comp-NeuroIMAGE_PRSice-output_cov_perm', pattern = 'PRSice_SUBAFF_EVER_\\s*(.*?)\\s*.summary', full.names = TRUE)

summary_files_all <- do.call(rbind,lapply(filenames_all,read.csv,sep=""))

summary_files_ADHD <- as.data.frame(lapply(filenames1,function(x){
  read.csv(x, sep="\t", header=TRUE)
}))
summary_files_SUBAFF <- as.data.frame(lapply(filenames2,function(x){
  read.csv(x, sep="\t", header=TRUE)
}))

############write list of unadjusted pvalues - this is hard coded so always check this!!
############ Update 12.07.22: Found a much better way to do this then below, only do at your own risk. Will be updated maybe. Still works though ###########

emp_pvals_ADHD <- t(as.data.frame(summary_files_ADHD[,c(12, 24, 36, 48, 60, 72, 84, 96, 108, 120, 132, 144, 156, 168, 180, 192, 204, 216, 228, 240, 252, 264, 276, 288, 300,
		 312, 324, 336, 348, 360, 372, 384, 396, 408, 420, 432, 444, 456, 468, 480, 492, 504, 516, 528, 540, 552, 564, 576, 588, 600)], header=FALSE))
emp_pvals_ADHD <- cbind(comp_names, emp_pvals_ADHD)

emp_pvals_SUBAFF <- t(as.data.frame(summary_files_SUBAFF[,c(12, 24, 36, 48, 60, 72, 84, 96, 108, 120, 132, 144, 156, 168, 180, 192, 204, 216, 228, 240, 252, 264, 276, 288, 300,
                 312, 324, 336, 348, 360, 372, 384, 396, 408, 420, 432, 444, 456, 468, 480, 492, 504, 516, 528, 540, 552, 564, 576, 588, 600)], header=FALSE))
emp_pvals_SUBAFF <- cbind(comp_names, emp_pvals_SUBAFF)

emp_pvals <- merge(emp_pvals_ADHD, emp_pvals_SUBAFF)

write.csv(emp_pvals, paste0(outdir, 'emp_pvals_unadj.csv'), row.names = FALSE)

##########adjust for false discovery rate and save to new file

emp_pvals_stacked <- cbind(emp_pvals[1], stack(emp_pvals[2:3]))
colnames(emp_pvals_stacked) <- c('Comp','Pval','Phen')

emp_pvals_adj <- p.adjust(emp_pvals_stacked$Pval, method = "fdr", n = length(emp_pvals_stacked$Pval))
emp_pvals_adj <- cbind(emp_pvals_adj,emp_pvals_stacked$Comp,emp_pvals_stacked$Phen)
colnames(emp_pvals_adj) <- c('Pval', 'Comp', 'Phen')

write.csv(emp_pvals_adj, paste0(outdir, 'emp_pvals_adj.csv'), row.names = FALSE)


############## Create table of all PRSice results, including adjusted pvals


results_all <- cbind(summary_files_all, emp_pvals_adj)
names(results_all)[names(results_all) == 'Pval'] <- 'FDR_adj_Empirical.P'

write.csv(results_all, paste0(outdir, 'all_results.csv'), row.names = FALSE)

##########perform association testing with potentially confounding variables
## create data files

dat <- read.table("/projects/0/einf2700/oblongl/PRSice_stats/Mcovall_Dx.txt", header = FALSE, sep = "", dec = ".")
colnames(dat) <- dat[1,]
dat <- dat[-1,]

##fetch PRS of components of interest, here IC2 and IC21. Could be any, always check this!! ##

filenames_scores_IC2 <- list.files(path ='/projects/0/einf2700/oblongl/polygen_comp-NeuroIMAGE_PRSice-output_cov_perm', pattern = 'PRSice_ADHD_EVER_ICA_ic_1.best', full.names = TRUE)
filenames_scores_IC21 <- list.files(path ='/projects/0/einf2700/oblongl/polygen_comp-NeuroIMAGE_PRSice-output_cov_perm', pattern = 'PRSice_SUBAFF_EVER_ICA_ic_21.best', full.names = TRUE)
filenames_scores_PC11 <- list.files(path ='/projects/0/einf2700/oblongl/polygen_comp-NeuroIMAGE_PRSice-output_cov_perm', pattern = 'PRSice_SUBAFF_EVER_PCA_pc_11.best', full.names = TRUE)

scores.IC2 <- do.call(cbind, lapply(filenames_scores_IC2,read.csv,sep=""))
scores.IC21 <- do.call(cbind,lapply(filenames_scores_IC21,read.csv,sep=""))
scores.PC11 <- do.call(cbind,lapply(filenames_scores_PC11,read.csv,sep=""))

colnames(scores.IC2)[4] <- "PRS.IC2"
colnames(scores.IC21)[4] <- "PRS.IC21"
colnames(scores.PC11)[4] <- "PRS.PC11"

#write.csv(scores.IC2, paste0(outdir, 'ADHD_PRS.csv'), row.names = FALSE)
#write.csv(scores.IC21, paste0(outdir, 'SUBAFF_PRS.csv'), row.names = FALSE)

FID <- dat[,1:2]
data1 <- sapply(dat[,3:24], as.numeric,omit.NA=TRUE)

data <- cbind.data.frame(FID, data1)
colnames(data)[1:2] <- c("FamID","FID")

data <- merge(data, scores.IC2, by = c("FID" = "FID"))
data <- merge(data, scores.IC21, by = c("FID" = "FID"))
data <- merge(data, scores.PC11, by = c("FID" = "FID"))

write.csv(data, paste0(outdir, 'regression_file.csv'), row.names=FALSE)

data$FamID <- as.factor(data$FamID)

###if you want to scale the polygenic scores run the code below, if not leave it out

data$PRS.IC2 <- data$PRS.IC2 * 1000
data$PRS.IC21 <- data$PRS.IC21 * 1000
data$PRS.PC11 <- data$PRS.PC11 * 1000
#perform regression

regIC2 <- glmer(ADHD_EVER ~ PRS.IC2 + Sex + Site + Batch + C1 + C2 + C3 + C4 + (1 + PRS.IC2 | FamID), family = binomial(link = "logit"), data = data, na.action = na.omit)
regIC21 <- glmer(SUBAFF_EVER ~ PRS.IC21 + Sex + Site + Batch + C1 + C2 + C3 + C4 + (1 + PRS.IC21 | FamID), family = binomial(link = "logit"), data = data, na.action = na.omit)
regPC11 <- glmer(SUBAFF_EVER ~ PRS.PC11 + Sex + Site + Batch + C1 + C2 + C3 + C4 + (1 + PRS.PC11 | FamID), family = binomial(link = "logit"), data = data, na.action = na.omit)


