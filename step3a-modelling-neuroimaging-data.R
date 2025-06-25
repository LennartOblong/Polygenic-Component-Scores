# Script to fit models for neuroimaging derived phenotypes using polygenic component scores as predictors

####### directories
indir <- "path/to/data/previously/created" #full phenotype, covariate and PCS data, and a set of thresholded phenotype data
outdir <- "path/to/desired/output"

####### library
library(data.table)
library(dplyr)
library(tidyverse)
library(stringr)
library(rms)
library(DescTools)
library(reshape2)
library(foreach)  #####this script can be easily parallelized if desired
library(doParallel)

####### get data created in previous scripts
cat("Loading data ...","\n")

data1 <- read.table(paste0(indir, "/full_imaging_data_ICthresh2.txt"), fill = TRUE, header = TRUE)
data2 <- read.table(paste0(indir, "/full_imaging_data_PCthresh2.txt"), fill = TRUE, header = TRUE)

phenofileica2 <- read.table(paste0(indir, "/fullpheno_ICthresh2.txt"), fill = TRUE, header = TRUE)
phenofilepca2 <- read.table(paste0(indir, "/fullpheno_PCthresh2.txt"), fill = TRUE, header = TRUE)

####### association testing
t1 <- "t1_"

cat("Commencing association testing with genomic independent component derived PCSs ...","\n")

results1 <- list()
for (phenotype in colnames(phenofileica2)) {
  if (!all(is.na(phenofileica2[[phenotype]]))) {  #check if phenotype is all NA
    if (grepl("t1_", phenotype)) {
        cat("Processing", phenotype, "\n")
        covar <- c("Assesment_center.0.0", "Age_at_assesment.0.0", "Genotype_batch.0.0", "BMI.0.0", "Sex.0.0", paste("Principal_Components.0.", 1:40, sep = ""),
                    paste(t1 ,"25756", sep = ""), paste(t1, "25757", sep = ""),  paste(t1, "25758", sep = ""), paste(t1, "25759", sep = ""),
                    paste(t1 ,"25000", sep = ""), paste(t1 ,"25741", sep = ""), paste(t1 ,"25742", sep = ""))
        model_prspc1 <- lm(as.formula(paste0(phenotype, " ~ ", paste0("genIC", 1:16, "_PRS_PC1", collapse = " + "), " + ", paste(covar, collapse = " + "))), data = data1, na.action = na.omit)
        null_model <- lm(as.formula(paste0(phenotype, " ~ ", paste(covar, collapse = " + "))), data = data1, na.action = na.omit)
        model_prspc1_r2 <- summary(model_prspc1)$r.squared
        null_model_r2 <- summary(null_model)$r.squared

        rsquared_prs_pc1 <- data.frame(model_prspc1_r2, null_model_r2)
        rsquared_prs_pc1$ID <- phenotype
        rsquared_prs_pc1$PRS_R2 <- (rsquared_prs_pc1$model_prspc1_r2 - rsquared_prs_pc1$null_model_r2)
        rsquared_prs_pc1$PRS_R2_adj <- (rsquared_prs_pc1$model_prspc1_r2 - rsquared_prs_pc1$null_model_r2) / (1 - rsquared_prs_pc1$null_model_r2)

        p_val <- tryCatch({
            p <- anova(model_prspc1, null_model)[2,6]
        }, error = function(e) NA)

        rsquared_prs_pc1$pval_unadj <- p_val

        comp_pvals <- summary(model_prspc1)$coefficients[paste0("genIC", 1:16, "_PRS_PC1"), "Pr(>|t|)"]
        comp_pvals_df <- as.data.frame(t(comp_pvals))
        comp_beta <- summary(model_prspc1)$coefficients[paste0("genIC", 1:16, "_PRS_PC1"), "Estimate"]
        comp_beta_df <- as.data.frame(t(comp_beta))
        colnames(comp_pvals_df) <- paste0("pval_genIC", 1:16)
        colnames(comp_beta_df) <- paste0("beta_genIC", 1:16)
        results_part <- cbind(rsquared_prs_pc1, comp_pvals_df)
        results_full <- cbind(results_part, comp_beta_df)
        #write results to output
        results1[[length(results1) + 1]] <- results_full

    } else {
    }
  } else {
    cat(phenotype, " is all NAs. Skipping.\n")
  }
}
if (length(results1) > 0) {
  final_results1 <- do.call(rbind, results1)
  output_file <- paste0(outdir, "/oIC_thresh2_imaging_results.txt")
  if (file.exists(output_file)) {
    write.table(final_results1, file = output_file, sep = "\t", row.names = FALSE, col.names = TRUE, append = FALSE, quote = FALSE)
  } else {
    write.table(final_results1, file = output_file, sep = "\t", row.names = FALSE, col.names = TRUE, append = TRUE, quote = FALSE)
  }
} else {
  cat("No results to write.\n")
}

cat("Done...","\n")

#### association testing using genomic principal component PCSs

cat("Commencing association testing with genomic principal component derived PCSs ...","\n")

results2 <- list()
for (phenotype in colnames(phenofilepca2)) {
  if (!all(is.na(phenofilepca2[[phenotype]]))) {  #check if phenotype is all NA
    if (grepl("t1_", phenotype)) {

        covar <- c("Assesment_center.0.0", "Age_at_assesment.0.0", "Genotype_batch.0.0", "BMI.0.0", "Sex.0.0", paste("Principal_Components.0.", 1:40, sep = ""),
                    paste(t1 ,"25756", sep = ""), paste(t1, "25757", sep = ""),  paste(t1, "25758", sep = ""), paste(t1, "25759", sep = ""),
                    paste(t1 ,"25000", sep = ""), paste(t1 ,"25741", sep = ""), paste(t1 ,"25742", sep = ""))
        model_prspc1 <- lm(as.formula(paste0(phenotype, " ~ ", paste0("genPC", 1:16, "_PRS_PC1", collapse = " + "), " + ", paste(covar, collapse = " + "))), data = data2, na.action = na.omit)
        null_model <- lm(as.formula(paste0(phenotype, " ~ ", paste(covar, collapse = " + "))), data = data2, na.action = na.omit)
        model_prspc1_r2 <- summary(model_prspc1)$r.squared
        null_model_r2 <- summary(null_model)$r.squared

        rsquared_prs_pc1 <- data.frame(model_prspc1_r2, null_model_r2)
        rsquared_prs_pc1$ID <- phenotype
        rsquared_prs_pc1$PRS_R2 <- (rsquared_prs_pc1$model_prspc1_r2 - rsquared_prs_pc1$null_model_r2)
        rsquared_prs_pc1$PRS_R2_adj <- (rsquared_prs_pc1$model_prspc1_r2 - rsquared_prs_pc1$null_model_r2) / (1 - rsquared_prs_pc1$null_model_r2)

        p_val <- tryCatch({
            p <- anova(model_prspc1, null_model)[2,6]
        }, error = function(e) NA)

        rsquared_prs_pc1$pval_unadj <- p_val

        comp_pvals <- summary(model_prspc1)$coefficients[paste0("genPC", 1:16, "_PRS_PC1"), "Pr(>|t|)"]
        comp_pvals_df <- as.data.frame(t(comp_pvals))
        comp_beta <- summary(model_prspc1)$coefficients[paste0("genPC", 1:16, "_PRS_PC1"), "Estimate"]
        comp_beta_df <- as.data.frame(t(comp_beta))
        colnames(comp_pvals_df) <- paste0("pval_genPC", 1:16)
        colnames(comp_beta_df) <- paste0("beta_genPC", 1:16)
        results_part <- cbind(rsquared_prs_pc1, comp_pvals_df)
        results_full <- cbind(results_part, comp_beta_df)

        #write results to output
        results2[[length(results2) + 1]] <- results_full

    } else {
    }
  } else {
    cat(phenotype, " is all NAs. Skipping.\n")
  }
}
if (length(results2) > 0) {
  final_results2 <- do.call(rbind, results2)
  output_file <- paste0(outdir, "/pca_thresh2_imaging_results.txt")
  if (file.exists(output_file)) {
    write.table(final_results2, file = output_file, sep = "\t", row.names = FALSE, col.names = TRUE, append = FALSE, quote = FALSE)
  } else {
    write.table(final_results2, file = output_file, sep = "\t", row.names = FALSE, col.names = TRUE, append = TRUE, quote = FALSE)
  }
} else {
  cat("No results to write.\n")
}

cat("All done. Results written to output directory.","\n")
cat("Have a nice day.")
