# Script to run pairwise-regression on the binary/binarized and continuous non-neuroimaging phenotypes.
# The following code needs the datafiles created by the previous scripts. This code is used to take the prepared datafiles as input and run a parallelized analysis of the results.

# This code is set up to take a split up version


#######arguments
args <- commandArgs(trailingOnly = TRUE)

inputfile <- args[1]
outfile <- args[2]

#######dir
#Note that the output of the analysis will be stored in the outfile passed to this script, defined in the wrapper script.

outdir <- "/path/to/desired/outdir" ###always check!!

#######library

library(data.table)
library(dplyr)
library(tidyverse)
library(stringr)
library(foreach)
library(doParallel)
library(rms)
library(MASS)
library(nnet)
library(DescTools)
library(fastDummies)

######set up backend
######NOTE: Setting up the backend will require the estimation of the cores used. This, however, poses some memory issues:
######1. Parallelizing means copying the required data. Here this means that for each regression, a large base datafile must be copied x times (number of cores). Thus, parallelizing is only
######   efficient up to the point where the job won't get immediately oom killed because the data is being copied to many times and therefore exceeds memory limit.
######2. The Snellius cluster will assign x cores for y memory requested proportionally, meaning that when parallelizing it becomes impossible to use all cores assigned due to the above mentioned
######   memory limitation. While this poses no issues for smaller tasks it becomes an issue when working with large datasets like the UKBB.

###### These issues are solved by splitting the output files from the previous code into chunks, and running this parallel script on each chunk seperately.

# set up backend
cl <- makeCluster(12) #make feasable cluster out of available cores
registerDoParallel(cl)

cat("Loading datasets...","\n")

#get input from args
input_filename <- basename(inputfile)

#load data based on filename and store in respective dataframe
if (grepl("^cleandata_ICnonimage_bin_split_", input_filename)) {
    cat("Loading data1_bin...\n")
    data1_bin <- fread(inputfile, fill = TRUE, header = TRUE)
    data1_bin[, (grep("^t1_", names(data1_bin), value = TRUE)) := lapply(.SD, as.factor), .SDcols = grep("^t1_", names(data1_bin), value = TRUE)]
} else if (grepl("^cleandata_PCnonimage_bin_split_", input_filename)) {
    cat("Loading data2_bin...\n")
    data2_bin <- fread(inputfile, fill = TRUE, header = TRUE)
    data2_bin[, (grep("^t1_", names(data2_bin), value = TRUE)) := lapply(.SD, as.factor), .SDcols = grep("^t1_", names(data2_bin), value = TRUE)]
} else if (grepl("^cleandata_ICnonimage_cont_split_", input_filename)) {
    cat("Loading data1_cont...\n")
    data1_cont <- fread(inputfile, fill = TRUE, header = TRUE)
    data1_cont[, (grep("^t1_", names(data1_cont), value = TRUE)) := lapply(.SD, as.numeric), .SDcols = grep("^t1_", names(data1_cont), value = TRUE)]
} else if (grepl("^cleandata_PCnonimage_cont_split_", input_filename)) {
    cat("Loading data2_cont...\n")
    data2_cont <- fread(inputfile, fill = TRUE, header = TRUE)
    data2_cont[, (grep("^t1_", names(data2_cont), value = TRUE)) := lapply(.SD, as.numeric), .SDcols = grep("^t1_", names(data2_cont), value = TRUE)]
} else {
    stop("Unrecognized file pattern. Please check the input file name.")
}

#columns to filter. The models do not run correctly or have not been binarized.
alter_columns <- c(
  paste0("t1_Cancer_register_-_Type_of_cancer.0"),
  paste0("t1_Diagnoses_-_main_ICD10.", 0:10),
  paste0("t1_Diagnoses_-_secondary_ICD10.", 0:24),
  paste0("t1_Number_of_self-reported_cancers.0")
)

#for now we remove the specified columns (could be used in seperate analysis, but binarization is going to need a lot of memory)
if (exists("data1_bin")) {
  toremove <- intersect(alter_columns, names(data1_bin))
  if (length(toremove) > 0) {
    data1_bin[, (toremove) := NULL]
  }
}

if (exists("data2_bin")) {
  toremove <- intersect(alter_columns, names(data2_bin))
  if (length(toremove) > 0) {
    data2_bin[, (toremove) := NULL]
  }
}


########## re-binarize all binary columns such that they take the values "0" and "1". Also enforces consistency across phenotypes. save original value mapping.
########## ensures that phenotypes already coded 0 1 remain this way.
original_value_mapping <- list()

if (exists("data1_bin")) {
  for (col_name in grep("^t1_", names(data1_bin), value = TRUE)) {
    unique_values <- sort(unique(na.omit(data1_bin[[col_name]])))  #sort unique values

    if (length(unique_values) == 2) {
      #if column is already 0 and 1 keep it unchanged
      if (all(unique_values == c(0, 1))) {
        original_value_mapping[[col_name]] <- setNames(unique_values, c("0", "1"))
      } else {
        #else map smallest value to 0 and the largest to 1 so that no flip occurs
        original_value_mapping[[col_name]] <- setNames(unique_values, c("0", "1"))
        data1_bin[[col_name]] <- factor(ifelse(data1_bin[[col_name]] == unique_values[1], "0", "1"), levels = c("0", "1"))
      }
    }
  }
}

if (exists("data2_bin")) {
  for (col_name in grep("^t1_", names(data2_bin), value = TRUE)) {
    unique_values <- sort(unique(na.omit(data2_bin[[col_name]])))  #sort unique values

    if (length(unique_values) == 2) {
      #if column is already 0 and 1 keep it unchanged
      if (all(unique_values == c(0, 1))) {
        original_value_mapping[[col_name]] <- setNames(unique_values, c("0", "1"))
      } else {
        #else map smallest value to 0 and the largest to 1 so that no flip occurs
        original_value_mapping[[col_name]] <- setNames(unique_values, c("0", "1"))
        data2_bin[[col_name]] <- factor(ifelse(data2_bin[[col_name]] == unique_values[1], "0", "1"), levels = c("0", "1"))
      }
    }
  }
}

#save a file that records the mapping for each binary phenotype
mapping_file_path <- file.path(outdir, paste0(basename(inputfile), "_value_mapping.rds"))
saveRDS(original_value_mapping, file = mapping_file_path)
cat("Original value mappings saved to:", mapping_file_path, "\n")

#define phenotypes that are to be tested based on input file to avoid testing covariates

valid_phenotypes <- NULL

if (exists("data1_bin")) {
  valid_phenotypes <- grep("^t1_", colnames(data1_bin), value = TRUE)
} else if (exists("data2_bin")) {
  valid_phenotypes <- grep("^t1_", colnames(data2_bin), value = TRUE)
} else if (exists("data1_cont")) {
  valid_phenotypes <- grep("^t1_", colnames(data1_cont), value = TRUE)
} else if (exists("data2_cont")) {
  valid_phenotypes <- grep("^t1_", colnames(data2_cont), value = TRUE)
}

#initialize output
output_file <- file.path(outfile)
write.table(NULL, file = output_file, sep="\t", row.names = FALSE, quote = FALSE)

cat("Commencing association testing of genomic independent component PCSs ...", "\n")
########NOTE: The PseudoR2() function of Desctools is used for the calculation of Mcfadden's and Nagelkerke's pseudo-R2

#### IC PCSs in binary phenotypes
results_bin1 <- list()
results_bin1 <- foreach(phenotype = valid_phenotypes, .packages = c("foreach","fastDummies", "DescTools", "stats", "stringr"), .combine = 'rbind') %dopar% {
#  if (phenotype %in% c("eid")) { next } #skip eid in every case
  if (exists("data1_bin")) {  #check if datafile exists
    original_phenotype <- phenotype
    phenotype <- str_replace_all(phenotype, c("\\+" = "plus", "[()]" = "", "/" = "")) #remove parentheses and dashes from column names. Replace + with 'plus'. Needed to create pairwise columns accurately
    unique_values <- unique(na.omit(data1_bin[[original_phenotype]])) ##determine unique values omitting NAs
    if (length(unique_values) > 2) {  #check if phenotype is not binary. create dummy cols if not binary
      data1_bin_mod <- data1_bin
      colnames(data1_bin_mod) <- str_replace_all(colnames(data1_bin_mod), c("\\+" = "plus", "[()]" = "", "/" = "")) #remove parentheses and dashes from column names. Replace + with 'plus'.
      data1_bin_mod <- fastDummies::dummy_cols(data1_bin_mod, select_columns = phenotype, remove_first_dummy = FALSE)
      dummy_columns <- grep(paste0("^", phenotype, "_"), names(data1_bin_mod), value = TRUE)
      dummy_columns <- dummy_columns[!grepl("_NA$", dummy_columns)]
      pairwise_columns <- c()
      for (i in 1:(length(dummy_columns) - 1)) {
        for (j in (i + 1):length(dummy_columns)) {

        category_j_name <- sub(".*_", "", dummy_columns[j])

        new_col_name <- paste0(dummy_columns[i], "_vs_", category_j_name)

        data1_bin_mod[[new_col_name]] <- ifelse(
        data1_bin_mod[[dummy_columns[i]]] == 1, 1,
        ifelse(data1_bin_mod[[dummy_columns[j]]] == 1, 0, NA)
          )
        pairwise_columns <- c(pairwise_columns, new_col_name)
        }
      }

      foreach(binary_col = pairwise_columns, .combine = 'rbind', .packages = c("foreach", "fastDummies", "DescTools", "stats")) %do% {
        if (any(data1_bin_mod[[binary_col]] < 0 | data1_bin_mod[[binary_col]] > 1, na.rm = TRUE)) { #check if the dummy cols contain values outside 0 to 1 range
                    stop(paste("Invalid values found in", binary_col))
                }
        full_phenotype_name <- binary_col
        if (grepl("t1_", phenotype) && full_phenotype_name %in% names(data1_bin_mod)) {

          covar <- c("Assesment_center.0.0", "Age_at_assesment.0.0", "Genotype_batch.0.0", "BMI.0.0", "Sex.0.0", paste("Principal_Components.0.", 1:40, sep = ""))

          model_prspc1_logit <- glm(as.formula(paste0("`", full_phenotype_name, "` ~ ", paste0("genIC", 1:16, "_PRS_PC1", collapse = " + "), " + ", paste(covar, collapse = " + "))),
                                        family = binomial(link="logit"), data = data1_bin_mod, na.action = na.omit)

          null_model_logit <- glm(as.formula(paste0("`", full_phenotype_name, "` ~ ", paste(covar, collapse = " + "))), family = binomial(link = "logit"), data = data1_bin_mod, na.action = na.omit)

          rsquared_prs_pc1 <- data.frame(
            model_prspc1_logit_McFr2 = tryCatch(PseudoR2(model_prspc1_logit, "McFadden"), error = function(e) NA),
            null_model_logit_McFr2 = tryCatch(PseudoR2(null_model_logit, "McFadden"), error = function(e) NA),
            model_prspc1_logit_Nagr2 = tryCatch(PseudoR2(model_prspc1_logit, "Nagel"), error = function(e) NA),
            null_model_logit_Nagr2 = tryCatch(PseudoR2(null_model_logit, "Nagel"), error = function(e) NA)
          )
          rsquared_prs_pc1$ID <- full_phenotype_name
          rsquared_prs_pc1$original_phenotype <- original_phenotype
          rsquared_prs_pc1$PRS_McF_pseudoR2 <- rsquared_prs_pc1$model_prspc1_logit_McFr2 - rsquared_prs_pc1$null_model_logit_McFr2
          rsquared_prs_pc1$PRS_McF_pseudoR2_adj <- rsquared_prs_pc1$PRS_McF_pseudoR2 / (1 - rsquared_prs_pc1$null_model_logit_McFr2)
          rsquared_prs_pc1$PRS_Nag_pseudoR2 <- rsquared_prs_pc1$model_prspc1_logit_Nagr2 - rsquared_prs_pc1$null_model_logit_Nagr2
          rsquared_prs_pc1$PRS_Nag_pseudoR2_adj <- rsquared_prs_pc1$PRS_Nag_pseudoR2 / (1 - rsquared_prs_pc1$null_model_logit_Nagr2)

          #calculate p-value
          p_val <- tryCatch({
            p <- anova(model_prspc1_logit, null_model_logit, test = "Chisq")[2,5]
          }, error = function(e) NA)

          rsquared_prs_pc1$pval_unadj <- p_val
          #fetch pvalues and log-odds for individual components
          comp_pvals <- summary(model_prspc1_logit)$coefficients[paste0("genIC", 1:16, "_PRS_PC1"), "Pr(>|z|)"]
          comp_pvals_df <- as.data.frame(t(comp_pvals))
          comp_logodds <- summary(model_prspc1_logit)$coefficients[paste0("genIC", 1:16, "_PRS_PC1"), "Estimate"]
          comp_logodds_df <- as.data.frame(t(comp_logodds))
          colnames(comp_pvals_df) <- paste0("pval_genIC", 1:16)
          colnames(comp_logodds_df) <- paste0("log-odds_genIC", 1:16)
          #bind all into one df
          results_part <- cbind(rsquared_prs_pc1, comp_pvals_df)
          results_full <- cbind(results_part, comp_logodds_df)
          #write result directly to output file
          return(results_full)
        } else {
          cat("Skipping timepoint until decision is made for", full_phenotype_name, ".\n")
        }
      }
    } else {  #phenotype is already binary
      if (grepl("t1_", original_phenotype) && original_phenotype %in% names(data1_bin)) {
        covar <- c("Assesment_center.0.0", "Age_at_assesment.0.0", "Genotype_batch.0.0", "BMI.0.0", "Sex.0.0", paste("Principal_Components.0.", 1:40, sep = ""))

        model_prspc1_logit <- glm(as.formula(paste0("`", original_phenotype, "` ~ ", paste0("genIC", 1:16, "_PRS_PC1", collapse = " + "), " + ", paste(covar, collapse = " + "))),
                                        family = binomial(link="logit"), data = data1_bin, na.action = na.omit)

        null_model_logit <- glm(as.formula(paste0("`", original_phenotype, "` ~ ", paste(covar, collapse = " + "))), family = binomial(link = "logit"), data = data1_bin, na.action = na.omit)

        rsquared_prs_pc1 <- data.frame(
          model_prspc1_logit_McFr2 = tryCatch(PseudoR2(model_prspc1_logit, "McFadden"), error = function(e) NA),
          null_model_logit_McFr2 = tryCatch(PseudoR2(null_model_logit, "McFadden"), error = function(e) NA),
          model_prspc1_logit_Nagr2 = tryCatch(PseudoR2(model_prspc1_logit, "Nagel"), error = function(e) NA),
          null_model_logit_Nagr2 = tryCatch(PseudoR2(null_model_logit, "Nagel"), error = function(e) NA)
        )

        rsquared_prs_pc1$ID <- phenotype
        rsquared_prs_pc1$original_phenotype <- original_phenotype
        rsquared_prs_pc1$PRS_McF_pseudoR2 <- rsquared_prs_pc1$model_prspc1_logit_McFr2 - rsquared_prs_pc1$null_model_logit_McFr2
        rsquared_prs_pc1$PRS_McF_pseudoR2_adj <- rsquared_prs_pc1$PRS_McF_pseudoR2 / (1 - rsquared_prs_pc1$null_model_logit_McFr2)
        rsquared_prs_pc1$PRS_Nag_pseudoR2 <- rsquared_prs_pc1$model_prspc1_logit_Nagr2 - rsquared_prs_pc1$null_model_logit_Nagr2
        rsquared_prs_pc1$PRS_Nag_pseudoR2_adj <- rsquared_prs_pc1$PRS_Nag_pseudoR2 / (1 - rsquared_prs_pc1$null_model_logit_Nag

        #calculate p-value
        p_val <- tryCatch({
          p <- anova(model_prspc1_logit, null_model_logit, test = "Chisq")[2,5]
        }, error = function(e) NA)

        rsquared_prs_pc1$pval_unadj <- p_val
          comp_pvals <- summary(model_prspc1_logit)$coefficients[paste0("genIC", 1:16, "_PRS_PC1"), "Pr(>|z|)"]
          comp_pvals_df <- as.data.frame(t(comp_pvals))
          comp_logodds <- summary(model_prspc1_logit)$coefficients[paste0("genIC", 1:16, "_PRS_PC1"), "Estimate"]
          comp_logodds_df <- as.data.frame(t(comp_logodds))
          colnames(comp_pvals_df) <- paste0("pval_genIC", 1:16)
          colnames(comp_logodds_df) <- paste0("log-odds_genIC", 1:16)
          results_part <- cbind(rsquared_prs_pc1, comp_pvals_df)
          results_full <- cbind(results_part, comp_logodds_df)
          #write result directly to output file
          return(results_full)

      } else {
        cat("Skipping timepoint in favor of larger sample.\n")
      }
    }
  } else {
    cat("data1_bin does not exist. Skipping the phenotype processing.\n")
  }
# }, error = function(e) {
# message("Error in processing phenotype: ", phenotype, " Error message: ", e$message)
# return(NULL)
# })
}
if (!is.null(results_bin1) && nrow(results_bin1) > 0) {
  write.table(results_bin1, file = output_file, sep = "\t", row.names = FALSE, col.names = TRUE, append = TRUE, quote = FALSE)
}

##### PC PCSs in binary phenotypes

results_bin2 <- list()
results_bin2 <- foreach(phenotype = valid_phenotypes, .packages = c("foreach","fastDummies", "DescTools", "stats", "stringr"), .combine = 'rbind') %dopar% {
# tryCatch({
#  if (phenotype %in% c("eid")) { next } #skip eid in every case
  if (exists("data2_bin")) {  #check if datafile exists
    original_phenotype <- phenotype
    phenotype <- str_replace_all(phenotype, c("\\+" = "plus", "[()]" = "", "/" = ""))#remove parentheses and dashes from column names. Replace + with 'plus'. Needed to create pairwise columns accurately
    unique_values <- unique(na.omit(data2_bin[[original_phenotype]]))
    if (length(unique_values) > 2) {  #check if phenotype is not binary. create dummy cols if not binary
      data2_bin_mod <- data2_bin
      colnames(data2_bin_mod) <- str_replace_all(colnames(data2_bin_mod), c("\\+" = "plus", "[()]" = "", "/" = ""))
      data2_bin_mod <- fastDummies::dummy_cols(data2_bin_mod, select_columns = phenotype, remove_first_dummy = FALSE)
      dummy_columns <- grep(paste0("^", phenotype, "_"), names(data2_bin_mod), value = TRUE)
      dummy_columns <- dummy_columns[!grepl("_NA$", dummy_columns)]
      pairwise_columns <- c()
      for (i in 1:(length(dummy_columns) - 1)) {
        for (j in (i + 1):length(dummy_columns)) {

        category_j_name <- sub(".*_", "", dummy_columns[j])

        new_col_name <- paste0(dummy_columns[i], "_vs_", category_j_name)

        data2_bin_mod[[new_col_name]] <- ifelse(
        data2_bin_mod[[dummy_columns[i]]] == 1, 1,
        ifelse(data2_bin_mod[[dummy_columns[j]]] == 1, 0, NA)
          )
        pairwise_columns <- c(pairwise_columns, new_col_name)
        }
      }
      foreach(binary_col = pairwise_columns, .combine = 'rbind', .packages = c("foreach", "fastDummies", "DescTools", "stats")) %do% {
        if (any(data2_bin_mod[[binary_col]] < 0 | data2_bin_mod[[binary_col]] > 1, na.rm = TRUE)) { #check if the dummy cols contain values outside 0 to 1 range
                    stop(paste("Invalid values found in", binary_col))
                }
        full_phenotype_name <- binary_col
        if (grepl("t1_", phenotype) && full_phenotype_name %in% names(data2_bin_mod)) {

          covar <- c("Assesment_center.0.0", "Age_at_assesment.0.0", "Genotype_batch.0.0", "BMI.0.0", "Sex.0.0", paste("Principal_Components.0.", 1:40, sep = ""))

          model_prspc1_logit <- glm(as.formula(paste0("`", full_phenotype_name, "` ~ ", paste0("genPC", 1:16, "_PRS_PC1", collapse = " + "), " + ", paste(covar, collapse = " + "))),
                                         family = binomial(link="logit"), data = data2_bin_mod, na.action = na.omit)

          null_model_logit <- glm(as.formula(paste0("`", full_phenotype_name, "` ~ ", paste(covar, collapse = " + "))), family = binomial(link = "logit"), data = data2_bin_mod, na.action = na.omit)

          rsquared_prs_pc1 <- data.frame(
            model_prspc1_logit_McFr2 = tryCatch(PseudoR2(model_prspc1_logit, "McFadden"), error = function(e) NA),
            null_model_logit_McFr2 = tryCatch(PseudoR2(null_model_logit, "McFadden"), error = function(e) NA),
            model_prspc1_logit_Nagr2 = tryCatch(PseudoR2(model_prspc1_logit, "Nagel"), error = function(e) NA),
            null_model_logit_Nagr2 = tryCatch(PseudoR2(null_model_logit, "Nagel"), error = function(e) NA)
          )

          rsquared_prs_pc1$ID <- full_phenotype_name
          rsquared_prs_pc1$original_phenotype <- original_phenotype
          rsquared_prs_pc1$PRS_McF_pseudoR2 <- rsquared_prs_pc1$model_prspc1_logit_McFr2 - rsquared_prs_pc1$null_model_logit_McFr2
          rsquared_prs_pc1$PRS_McF_pseudoR2_adj <- rsquared_prs_pc1$PRS_McF_pseudoR2 / (1 - rsquared_prs_pc1$null_model_logit_McFr2)
          rsquared_prs_pc1$PRS_Nag_pseudoR2 <- rsquared_prs_pc1$model_prspc1_logit_Nagr2 - rsquared_prs_pc1$null_model_logit_Nagr2
          rsquared_prs_pc1$PRS_Nag_pseudoR2_adj <- rsquared_prs_pc1$PRS_Nag_pseudoR2 / (1 - rsquared_prs_pc1$null_model_logit_Nagr2)

          #calculate p-value
          p_val <- tryCatch({
            p <- anova(model_prspc1_logit, null_model_logit, test = "Chisq")[2,5]
          }, error = function(e) NA)

          rsquared_prs_pc1$pval_unadj <- p_val

          comp_pvals <- summary(model_prspc1_logit)$coefficients[paste0("genPC", 1:16, "_PRS_PC1"), "Pr(>|z|)"]
          comp_pvals_df <- as.data.frame(t(comp_pvals))
          comp_logodds <- summary(model_prspc1_logit)$coefficients[paste0("genPC", 1:16, "_PRS_PC1"), "Estimate"]
          comp_logodds_df <- as.data.frame(t(comp_logodds))
          colnames(comp_pvals_df) <- paste0("pval_genPC", 1:16)
          colnames(comp_logodds_df) <- paste0("log-odds_genPC", 1:16)
          results_part <- cbind(rsquared_prs_pc1, comp_pvals_df)
          results_full <- cbind(results_part, comp_logodds_df)
          #write result directly to output file
          return(results_full)

        } else {
          cat("Skipping timepoint until decision is made for", full_phenotype_name, ".\n")
        }
      }
    } else {  #phenotype is already binary
      if (grepl("t1_", original_phenotype) && original_phenotype %in% names(data2_bin)) {

        covar <- c("Assesment_center.0.0", "Age_at_assesment.0.0", "Genotype_batch.0.0", "BMI.0.0", "Sex.0.0", paste("Principal_Components.0.", 1:40, sep = ""))

        model_prspc1_logit <- glm(as.formula(paste0("`", original_phenotype, "` ~ ", paste0("genPC", 1:16, "_PRS_PC1", collapse = " + "), " + ", paste(covar, collapse = " + "))),
                                        family = binomial(link="logit"), data = data2_bin, na.action = na.omit)

        null_model_logit <- glm(as.formula(paste0("`", original_phenotype, "` ~ ", paste(covar, collapse = " + "))), family = binomial(link = "logit"), data = data2_bin, na.action = na.omit)

        rsquared_prs_pc1 <- data.frame(
          model_prspc1_logit_McFr2 = tryCatch(PseudoR2(model_prspc1_logit, "McFadden"), error = function(e) NA),
          null_model_logit_McFr2 = tryCatch(PseudoR2(null_model_logit, "McFadden"), error = function(e) NA),
          model_prspc1_logit_Nagr2 = tryCatch(PseudoR2(model_prspc1_logit, "Nagel"), error = function(e) NA),
          null_model_logit_Nagr2 = tryCatch(PseudoR2(null_model_logit, "Nagel"), error = function(e) NA)
        )

        rsquared_prs_pc1$ID <- phenotype
        rsquared_prs_pc1$original_phenotype <- original_phenotype
        rsquared_prs_pc1$PRS_McF_pseudoR2 <- rsquared_prs_pc1$model_prspc1_logit_McFr2 - rsquared_prs_pc1$null_model_logit_McFr2
        rsquared_prs_pc1$PRS_McF_pseudoR2_adj <- rsquared_prs_pc1$PRS_McF_pseudoR2 / (1 - rsquared_prs_pc1$null_model_logit_McFr2)
        rsquared_prs_pc1$PRS_Nag_pseudoR2 <- rsquared_prs_pc1$model_prspc1_logit_Nagr2 - rsquared_prs_pc1$null_model_logit_Nagr2
        rsquared_prs_pc1$PRS_Nag_pseudoR2_adj <- rsquared_prs_pc1$PRS_Nag_pseudoR2 / (1 - rsquared_prs_pc1$null_model_logit_Nagr2)

        # Calculate p-value
        p_val <- tryCatch({
          p <- anova(model_prspc1_logit, null_model_logit, test = "Chisq")[2,5]
        }, error = function(e) NA)


        rsquared_prs_pc1$pval_unadj <- p_val
          comp_pvals <- summary(model_prspc1_logit)$coefficients[paste0("genPC", 1:16, "_PRS_PC1"), "Pr(>|z|)"]
          comp_pvals_df <- as.data.frame(t(comp_pvals))
          comp_logodds <- summary(model_prspc1_logit)$coefficients[paste0("genPC", 1:16, "_PRS_PC1"), "Estimate"]
          comp_logodds_df <- as.data.frame(t(comp_logodds))
          colnames(comp_pvals_df) <- paste0("pval_genPC", 1:16)
          colnames(comp_logodds_df) <- paste0("log-odds_genPC", 1:16)
          results_part <- cbind(rsquared_prs_pc1, comp_pvals_df)
          results_full <- cbind(results_part, comp_logodds_df)
          #write result directly to output file
          return(results_full)
      } else {
        cat("Skipping timepoint in favor of larger sample.\n")
      }
    }
  } else {
    cat("data2_bin does not exist. Skipping the phenotype processing.\n")
  }
}
if (!is.null(results_bin2) && nrow(results_bin2) > 0) {
write.table(results_bin2, file = output_file, sep="\t", row.names = FALSE, col.names = TRUE, append = TRUE, quote = FALSE)
}

# IC PCSs in continuous phenotypes

results_cont1 <- list()
results_cont1 <- foreach(phenotype = valid_phenotypes, .packages = c("foreach","fastDummies", "DescTools", "stats", "stringr"), .combine = 'rbind') %dopar% {
  if (exists("data1_cont")) {  #check if dataframe exists
    if (grepl("t1_", phenotype) && phenotype %in% names(data1_cont)) {

        covar <- c("Assesment_center.0.0", "Age_at_assesment.0.0", "Genotype_batch.0.0", "BMI.0.0", "Sex.0.0", paste("Principal_Components.0.", 1:40, sep = ""))
        model_prspc1 <- lm(as.formula(paste0("`", phenotype, "` ~ ", paste0("genIC", 1:16, "_PRS_PC1", collapse = " + "), " + ", paste(covar, collapse = " + "))), data = data1_cont, na.action = na.omit)
        null_model <- lm(as.formula(paste0("`", phenotype, "` ~ ", paste(covar, collapse = " + "))), data = data1_cont, na.action = na.omit)
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
        return(results_full)
    } else {
        cat("Skipping timepoint until decision is made.\n")
    }
  } else {
    cat("data1_cont does not exist. Skipping.\n")
  }
}
if (!is.null(results_cont1) && nrow(results_cont1) > 0) {
write.table(results_cont1, file = output_file, sep="\t", row.names = FALSE, col.names = TRUE, append = TRUE, quote = FALSE)
}

#PC PCSs in continuous phenotypes

results_cont2 <- list()
results_cont2 <- foreach(phenotype = valid_phenotypes, .packages = c("foreach","fastDummies", "DescTools", "stats", "stringr"), .combine = 'rbind') %dopar% {
  if (exists("data2_cont")) {  #check if dataframe exists
    if (grepl("t1_", phenotype) && phenotype %in% names(data2_cont)) {

        covar <- c("Assesment_center.0.0", "Age_at_assesment.0.0", "Genotype_batch.0.0", "BMI.0.0", "Sex.0.0", paste("Principal_Components.0.", 1:40, sep = ""))
        model_prspc1 <- lm(as.formula(paste0("`", phenotype, "` ~ ", paste0("genPC", 1:16, "_PRS_PC1", collapse = " + "), " + ", paste(covar, collapse = " + "))), data = data2_cont, na.action = na.omit)
        null_model <- lm(as.formula(paste0("`", phenotype, "` ~ ", paste(covar, collapse = " + "))), data = data2_cont, na.action = na.omit)
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
        return(results_full)
    } else {
        cat("Skipping timepoint until decision is made.\n")
    }
  } else {
    cat("data2_cont does not exist. Skipping.\n")
  }
}
if (!is.null(results_cont2) && nrow(results_cont2) > 0) {
write.table(results_cont2, file = output_file, sep="\t", row.names = FALSE, col.names = TRUE, append = TRUE, quote = FALSE)
}

#stop cluster
stopCluster(cl)

cat("All done.", "\n")
cat("Have a nice day.","\n")
