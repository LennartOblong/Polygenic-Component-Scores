# Script to create full datasets to run regression analysis.
# This script is subset into two sections: 1. Creating necessary files for neuroimaging phenotype analysis. 2. Creating necessary files for non-neuroimaging phenotype analysis.
# Both sections are based on two different PRS-PCA outputs, depending on the sample used.
# Throughout this script neuroimaging phenotype data is often referred to as "IDP", meaning "imaging derived phenotype".

###### NOTE #######
# This script is memory hungry particularly in the second section, due to the binarization and wrangling of large data.
# ~200GB memory required, ~3 hours runtime.

# source prs-pca script. This will run prs-pca on provided polygenic score p-thresholds (calculate prior to running step1 using the command in the description of step1), and store the output for all components in a list object.
# this will also define all necessary directories which contain keys, covariate data, phenotype data, etc.
source("path/to/scripts/step1-prs-pca.R")

####### get data for later scripts. R does not like number headers and puts an X in front of each number. This is removed to enable match with phenotype key.
####### this assumes the UK Biobank data naming convention e.g. 22222.1.1 (datafieldID.instance.array)
####### change the names of the files below to fit your naming convention.

cat("Loading neuroimaging data and covariates...", "\n")
covarfile <- read.table(paste0(indir, "/covariate_data.txt"), fill = TRUE, header = TRUE) #file containing covariate data
names(covarfile) <- gsub("X", "", names(covarfile))

covarkey <- read.table(paste0(indir, "/covar_identifiers.txt"), sep="\t", header=FALSE, colClasses=c("numeric", "character"), stringsAsFactors=FALSE) #file containing the covariate UKBB catalogue IDs next to the names of the covariates
covarkey[, 2] <- gsub(" ", "_", covarkey[, 2])

phenotype_key <- read.table(paste0(indir, "/imaging-ID-to-UKBB-ID.txt"), fill = TRUE, header = TRUE) #This file relates the IDs used in Smith et al., 2021 to the UKBB catalogue IDs

imgid <- read.table(paste0(indir, "/imaging_data_eids_2023.txt"), sep="\t", fill = TRUE, header = TRUE) #participant IDs that have imaging data

imgcovarfile <- read.table(paste0(indir, "/mri_specific_covariate_data.txt"), fill = TRUE, header = TRUE) #imaging only covariate data (e.g. headmotion in scanner)
names(imgcovarfile) <- gsub("X", "", names(imgcovarfile))

imgphenofile <- read.table(paste0(indir, "/mri_phenotype_data.txt"), fill = TRUE, header = TRUE) #neuroimaging phenotype data
names(phenofile) <- gsub("X", "", names(phenofile))

cat("Loading non-neuroimaging  data...", "\n")
pheno <- read.table(paste0(indir, "/non_mri_phenotype_data.csv"), fill = TRUE, header = TRUE) #non-neuroimaigng phenotype data
names(pheno) <- gsub("X", "", names(pheno))

phenoident <- read.table(paste0(indir, "/nonimage_key.csv"), sep="\t", fill = TRUE, header = TRUE) #file matching the non-neuroimaging phenotype UKBB IDs with the names of the phenotypes
phenoident[, 2] <- gsub(" ", "_", phenoident[, 2])


#Wrangle lists to create dataframe that contains the first PRS-PC for each genomic component (i.e. PCS). Prepare merging with phenotype files. (Requires list objects from previous script!)
#neuroimaging data PCSs. Note that 16 components where generated in our study, therefore it is hard coded.
cat("Getting salient polygenic component score per component ...", "\n")
prs_pca_oIC_img_output <- NULL

for (i in 1:16) {
  prs_pc1 <- output_list_img[[paste0("output_", i, "_oIC")]]$data$PRS_PC1

  if (is.null(prs_pca_oIC_img_output)) {
    prs_pca_oIC_img_output <- data.frame(IID = output_list_img[[paste0("output_", i, "_oIC")]]$data$IID, prs_pc1)
  } else {
    if (nrow(prs_pca_oIC_img_output) == length(prs_pc1)) {
      prs_pca_oIC_img_output <- cbind(prs_pca_oIC_img_output, prs_pc1)
    } else {
      warning(paste("Number of rows doesn't match for PRS_PC1 data in output_", i, "_oIC"))
    }
  }
}
colnames(prs_pca_oIC_img_output) <- c("eid", paste0("genIC", 1:16, "_PRS_PC1"))

prs_pca_pca_img_output <- NULL

for (i in 1:16) {
  prs_pc1 <- output_list_img[[paste0("output_", i, "_pca")]]$data$PRS_PC1

  if (is.null(prs_pca_pca_img_output)) {
    prs_pca_pca_img_output <- data.frame(IID = output_list_img[[paste0("output_", i, "_pca")]]$data$IID, prs_pc1)
  } else {
    if (nrow(prs_pca_pca_img_output) == length(prs_pc1)) {
      prs_pca_pca_img_output <- cbind(prs_pca_pca_img_output, prs_pc1)
    } else {
      warning(paste("Number of rows doesn't match for PRS_PC1 data in output_", i, "_pca"))
    }
  }
}
colnames(prs_pca_pca_img_output) <- c("eid", paste0("genPC", 16:1, "_PRS_PC1"))  ##Note, the naming is reversed to account for the reverse PC output of MELODIC during genPC generation

#non-neuroimaging data PCSs
prs_pca_oIC_nonimg_output <- NULL
for (i in 1:16) {
  prs_pc1 <- output_list_nonimg[[paste0("output_", i, "_IC")]]$data$PRS_PC1

  if (is.null(prs_pca_oIC_nonimg_output)) {
    prs_pca_oIC_nonimg_output <- data.frame(IID = output_list_nonimg[[paste0("output_", i, "_IC")]]$data$IID, prs_pc1)
  } else {
    if (nrow(prs_pca_oIC_nonimg_output) == length(prs_pc1)) {
      prs_pca_oIC_nonimg_output <- cbind(prs_pca_oIC_nonimg_output, prs_pc1)
    } else {
      warning(paste("Number of rows doesn't match for PRS_PC1 data in output_", i, "_IC"))
    }
  }
}
colnames(prs_pca_oIC_nonimg_output) <- c("eid", paste0("genIC", 1:16, "_PRS_PC1"))

prs_pca_pca_nonimg_output <- NULL
for (i in 1:16) {
  prs_pc1 <- output_list_nonimg[[paste0("output_", i, "_PC")]]$data$PRS_PC1

  if (is.null(prs_pca_pca_nonimg_output)) {
    prs_pca_pca_nonimg_output <- data.frame(IID = output_list_nonimg[[paste0("output_", i, "_PC")]]$data$IID, prs_pc1)
  } else {
    if (nrow(prs_pca_pca_nonimg_output) == length(prs_pc1)) {
      prs_pca_pca_nonimg_output <- cbind(prs_pca_pca_nonimg_output, prs_pc1)
    } else {
      warning(paste("Number of rows doesn't match for PRS_PC1 data in output_", i, "_PC"))
    }
  }
}
colnames(prs_pca_pca_nonimg_output) <- c("eid", paste0("genPC", 16:1, "_PRS_PC1"))  ##Note, the naming is reversed to account for the reverse PC output of MELODIC during genPC generation

cat("Done.", "\n")
############################# Data wrangling, changing headers and creating final dataframe for imaging data

cat("Data Wrangling neuroimaging data ...", "\n")

# function to match datafield IDs with corresponding description using a key filer
# this function again assumes the standard UK Biobank naming convention in the form of datafieldID.instance.array
replacenumbers <- function(colnames, lookup_table, pheno_filename, lookup_filename) {
  tracking_file <- paste0(ukb_dir, "/replacement_tracking2_", pheno_filename, "_", lookup_filename, ".csv")
  changes <- list()

  for (i in seq_along(colnames)) {
    original_name <- colnames[i]
    match <- sub("^(t[0-9]+_)([0-9]+)\\..*", "\\1\\2", colnames[i])

    if (grepl("^t[0-9]+_[0-9]+", match)) {
      prefix <- sub("(_[0-9]+).*", "\\1", match)
      lookup_number <- sub("^t[0-9]+_", "", match)

      for (j in seq_len(nrow(lookup_table))) {
        if (lookup_number == lookup_table[j, 1]) {
          if (grepl(":", lookup_table[j, 2])) {
            replace <- sub(":.*", "", lookup_table[j, 2])
          } else {
            replace <- lookup_table[j, 2]
          }
          colnames[i] <- sub(paste0("\\b", lookup_number, "\\b"), replace, colnames[i])
          changes[[length(changes) + 1]] <- data.frame(
            Input_Colname = original_name,
            Matched_String = lookup_table[j, 2],
            Output_Colname = colnames[i],
            stringsAsFactors = FALSE
          )
          break
        }
      }
    }
  }

  tracking_df <- do.call(rbind, changes)
  write.csv(tracking_df, file = tracking_file, row.names = FALSE)

  return(colnames)
}

#use function to name covariates
colnames(covarfile) <- replacenumbers(colnames(covarfile), lookuptable, "covariates", "key")

#prepare phenofile and imgcovarfile for function. Replaces the instance number with a prefix. In imaging the array indicator is always 0 and can thus be ignored and removed.
colnames(imgphenofile) <- ifelse(grepl("\\.2\\.0$", colnames(imgphenofile)), paste0("t1_", sub("\\.2\\.0$", "", colnames(imgphenofile))), colnames(imgphenofile))
colnames(imgphenofile) <- ifelse(grepl("\\.3\\.0$", colnames(imgphenofile)), paste0("t2_", sub("\\.3\\.0$", "", colnames(imgphenofile))), colnames(imgphenofile))

colnames(imgcovarfile) <- ifelse(grepl("\\.2\\.0$", colnames(imgcovarfile)), paste0("t1_", sub("\\.2\\.0$", "", colnames(imgcovarfile))), colnames(imgcovarfile))
colnames(imgcovarfile) <- ifelse(grepl("\\.3\\.0$", colnames(imgcovarfile)), paste0("t2_", sub("\\.3\\.0$", "", colnames(imgcovarfile))), colnames(imgcovarfile))

#bind headmotion phenotypes (QC) to covarfile to prevent them from being filtered out by IDP thresholding later
covar <- cbind(covarfile, phenofile[match(covarfile$eid, phenofile$eid), c("t1_25741", "t1_25742")])

#get identifier data used for thresholding IDPs, and thresholed IDP lists
ica2 <- read.table(paste0(data_dir, "/ICA_IDPloadings2.txt"), fill = TRUE, header = FALSE)
pca2 <- read.table(paste0(data_dir, "/PCA_IDPloadings2.txt"), fill = TRUE, header = FALSE)

ica2ident <- merge(ica2, phenotype_key, by.x = "V1", by.y = "Pheno")
pca2ident <- merge(pca2, phenotype_key, by.x = "V1", by.y = "Pheno")

#filter phenofile for the different IDP ica/pca threshold loaded above. The IDP threshold files need to be created prior to this (filter component IDP-loadings for desired threshold)
#IDP-loadings
phenofile_columns <- gsub("^t1_", "", colnames(imgphenofile))

ica2ident_ids <- ica2ident$UKB_ID
matching_ids_ica2 <- intersect(phenofile_columns, ica2ident_ids)
if (length(matching_ids_ica2) > 0) {
  matching_column_names <- c("eid", paste0("t1_", matching_ids_ica2))
  phenofileica2 <- imgphenofile[, matching_column_names, drop = FALSE]
} else {
  phenofileica2 <- NULL
}
phenofileica2 <- phenofileica2[, !colnames(phenofileica2) %in% c("t1_25755"), drop = FALSE]
write.table(phenofileica2, file.path(paste0(outdir, "/fullpheno_ICthresh2.txt")), sep="\t", row.names = FALSE)

pca2ident_ids <- pca2ident$UKB_ID
matching_ids_pca2 <- intersect(phenofile_columns, pca2ident_ids)
if (length(matching_ids_pca2) > 0) {
  matching_column_names <- c("eid", paste0("t1_", matching_ids_pca2))
  phenofilepca2 <- imgphenofile[, matching_column_names, drop = FALSE]
} else {
  phenofilepca2 <- NULL
}
phenofilepca2 <- phenofilepca2[, !colnames(phenofilepca2) %in% c("t1_25755"), drop = FALSE]
write.table(phenofilepca2, file.path(paste0(outdir, "/fullpheno_PCthresh2.txt")), sep="\t", row.names = FALSE)

cat("Neuroimaging data merging into complete dataset...", "\n")
data1 <- prs_pca_oIC_img_output %>%
  inner_join(phenofileica2, by = "eid") %>%
  inner_join(covar, by = "eid") %>%
  inner_join(imgcovarfile, by = "eid")
data1[is.na(data1) | data1=="Inf"] = NA
data2 <- prs_pca_pca_img_output %>%
  inner_join(phenofilepca2, by = "eid") %>%
  inner_join(covar, by = "eid") %>%
  inner_join(imgcovarfile, by = "eid")
data2[is.na(data2) | data2=="Inf"] = NA

#remove invalid phenotypes
data1 <- data1[, !colnames(data1) %in% c("t1_25755", "t2_25755"), drop = FALSE]
data2 <- data2[, !colnames(data2) %in% c("t1_25755", "t2_25755"), drop = FALSE]

cat("Saving final neuroimaging datafiles for later reference...", "\n")
write.table(data1, file.path(paste0(outdir, "/full_imaging_data_ICthresh2.txt")), sep="\t", row.names = FALSE)
write.table(data2, file.path(paste0(outdir, "/full_imaging_data_PCthresh2.txt")), sep="\t", row.names = FALSE)

##########################################################################################################

cat("Data Wrangling, quality control and decoding datafield names in non-neuroimaging data ...", "\n")

##########################################################################################################################################################
#Replace values to remove data coding oddities. Reorder values to make interpretation more consistent. Set values to NA to remove from later analysis.
#These data cleaning are consistent with the recommendations in Watanabe et al., 2019

cat("Data cleaning...", "\n")

column_prefixes <- sub("\\..*", "", names(pheno))

#replace -10 with 0 in columns where applicable
cat("Recoding datafields...", "\n")
pheno[column_prefixes == "699"] <- lapply(
  pheno[column_prefixes == "699"],
  function(x) ifelse(x == -10, 0, x)
)

#replace -7 with 0 in columns where applicable
pheno[column_prefixes %in% c("680", "670", "6164", "6162", "6160", "6159", "6157", "6152", "6150", "6149", "6146", "6145", "6143", "6142", "6139")] <- lapply(
  pheno[column_prefixes %in% c("680", "670", "6164", "6162", "6160", "6159", "6157", "6152", "6150", "6149", "6146", "6145", "6143", "6142", "6139")],
  function(x) ifelse(x == -7, 0, x)
)

#replace -5 with NA in columns where applicable and remove 2 and 3 from other fields.
pheno[column_prefixes == "3591"] <- lapply(
  pheno[column_prefixes == "3591"],
  function(x) ifelse(x == -5, NA, x)
)

pheno[column_prefixes == "2834"] <- lapply(
  pheno[column_prefixes == "2834"],
  function(x) ifelse(x == -5, NA, x)
)

pheno[column_prefixes == "2724"] <- lapply(
  pheno[column_prefixes == "2724"],
  function(x) ifelse(x %in% c(2, 3), NA, x)
)

#replace 5 with NA in columns where applicable
pheno[column_prefixes == "2267"] <- lapply(
  pheno[column_prefixes == "2267"],
  function(x) ifelse(x == 5, NA, x)
)

#field 20536: exclude phenotype code 3 and swap 2 to 0
pheno[column_prefixes == "20536"] <- lapply(
  pheno[column_prefixes == "20536"],
  function(x) {
    x[x == 3] <- NA
    x[x == 2] <- 0
    x
  }
)

#fields 20117 and 20116: replace -7 with NA
pheno[column_prefixes %in% c("20117", "20116")] <- lapply(
  pheno[column_prefixes %in% c("20117", "20116")],
  function(x) {
    x[x == -7] <- NA
    x
  }
)

#fields 20111, 20110 and 20107: replace -17 and -27 with 0
pheno[column_prefixes %in% c("20111", "20110", "20107")] <- lapply(
  pheno[column_prefixes %in% c("20111", "20110", "20107")],
  function(x) {
    x[x %in% c(-17, -27)] <- 0
    x
  }
)

#fields 20002: replace 99999 with 0
pheno[column_prefixes %in% c("20002")] <- lapply(
  pheno[column_prefixes %in% c("20002")],
  function(x) {
    x[x == 99999] <- 0
    x
  }
)

#field 1757: swap 2 (older) with 3 (average)
pheno[column_prefixes == "1757"] <- lapply(
  pheno[column_prefixes == "1757"],
  function(x) {
    x[x == 2] <- 99  #temp value
    x[x == 3] <- 2
    x[x == 99] <- 3  #restore temp value to 3
    x
  }
)

#field 1707: swap 2 (left) with 3 (both)
pheno[column_prefixes == "1707"] <- lapply(
  pheno[column_prefixes == "1707"],
  function(x) {
    x[x == 2] <- 99  #temp value
    x[x == 3] <- 2
    x[x == 99] <- 3  #restore temp value to 3
    x
  }
)

#field 1697: swap 2 (taller) with 3 (about average)
pheno[column_prefixes == "1697"] <- lapply(
  pheno[column_prefixes == "1697"],
  function(x) {
    x[x == 2] <- 99  #temp value
    x[x == 3] <- 2
    x[x == 99] <- 3  #restore temp value to 3
    x
  }
)

#field 1687: swap 2 (plumper) with 3 (about average)
pheno[column_prefixes == "1687"] <- lapply(
  pheno[column_prefixes == "1687"],
  function(x) {
    x[x == 2] <- 99  #temp value
    x[x == 3] <- 2
    x[x == 99] <- 3  #restore temp value to 3
    x
  }
)

#field 1618: remove -6 (it varies).
pheno[column_prefixes == "1618"] <- lapply(
  pheno[column_prefixes == "1618"],
  function(x) {
    x[x == -6] <- NA
    x
  }
)

#field 1518: remove -2 (do not drink hot drinks) and reverse values
pheno[column_prefixes == "1518"] <- lapply(
  pheno[column_prefixes == "1518"],
  function(x) {
    x[x == -2] <- NA
    x <- max(x, na.rm = TRUE) - x  #reverse values
    x
  }
)

#fields 1508: replace 4 with 0
pheno[column_prefixes %in% c("1508")] <- lapply(
  pheno[column_prefixes %in% c("1508")],
  function(x) {
    x[x == 4] <- 0
    x
  }
)

#field 1468: replace 5 with 0
pheno[column_prefixes == "1468"] <- lapply(
  pheno[column_prefixes == "1468"],
  function(x) {
    x[x == 5] <- 0
    x
  }
)

#fields 1448: replace 4 with 0
pheno[column_prefixes %in% c("1448")] <- lapply(
  pheno[column_prefixes %in% c("1448")],
  function(x) {
    x[x == 4] <- 0
    x
  }
)

#field 1428: replace 3 with 4 and remove 0
pheno[column_prefixes == "1428"] <- lapply(
  pheno[column_prefixes == "1428"],
  function(x) {
    x[x == 3] <- 4
    x[x == 0] <- NA
    x
  }
)

#fields 1418: replace 5 with 0 and remove 6
pheno[column_prefixes %in% c("1418")] <- lapply(
  pheno[column_prefixes %in% c("1418")],
  function(x) {
    x[x == 5] <- 0
    x[x == 6] <- NA
    x
  }
)

#field 1249: reverse phenotype codes
pheno[column_prefixes == "1249"] <- lapply(
  pheno[column_prefixes == "1249"],
  function(x) {
    x <- max(x, na.rm = TRUE) - x
    x
  }
)

#field 2178: reverse phenotype codes
pheno[column_prefixes == "2178"] <- lapply(
  pheno[column_prefixes == "2178"],
  function(x) {
    x <- max(x, na.rm = TRUE) - x
    x
  }
)


#field 1239: swap phenotype codes 1 and 2
pheno[column_prefixes == "1239"] <- lapply(
  pheno[column_prefixes == "1239"],
  function(x) {
    x[x == 1] <- 99  #temp value
    x[x == 2] <- 1
    x[x == 99] <- 2  #restore temp value to 2
    x
  }
)

#field 1150: swap 2 (right) and 3 (both) and reverse
pheno[column_prefixes == "1150"] <- lapply(
  pheno[column_prefixes == "1150"],
  function(x) {
    x[x == 2] <- 99  #temp value
    x[x == 3] <- 2
    x[x == 99] <- 3  #restore temp value to 3
    x <- max(x, na.rm = TRUE) - x  #reverse values
    x
  }
)

#field 1140: remove 3 and swap 0 and 1
pheno[column_prefixes == "1140"] <- lapply(
  pheno[column_prefixes == "1140"],
  function(x) {
    x[x == 3] <- NA
    x[x == 0] <- 99  #temp value
    x[x == 1] <- 0
    x[x == 99] <- 1  #restore temp value to 1
    x
  }
)

#####all field replacements for inconclusive data. Replace -3, -818, -1, -4, -11, -13, -21,and -23 with NA (removes "do not know", "do not recall","prefer not to answer" replies).
#####Also removes -10 values as they are typically few ("less then one of x a week").

rep <- c(-3, -818, -1, -4, -11, -13, -21, -23, -121, -10)

pheno[] <- lapply(pheno, function(column) {
  if (is.numeric(column)) {
    column[column %in% rep] <- NA
  }
  return(column)
})

############################################################################################################################################################

#replace pheno IDs with phenotype name. replace UKB timepoint coding (i.e. ".0", ".1", etc. with string prefixes (i.e. "t1_", "t2_", etc.)
#colnames(pheno) <- replacenumbers(colnames(pheno), phenoident, "phenotypes", "key")
cat("Adjusting datafield names in file...", "\n")
colnames(pheno) <- ifelse(grepl("\\.0\\.\\d+$", colnames(pheno)),
                          paste0("t1_", sub("\\.0\\.(\\d+)$", ".\\1", colnames(pheno))),
                          colnames(pheno))
colnames(pheno) <- ifelse(grepl("\\.1\\.\\d+$", colnames(pheno)),
                          paste0("t2_", sub("\\.1\\.(\\d+)$", ".\\1", colnames(pheno))),
                          colnames(pheno))
colnames(pheno) <- ifelse(grepl("\\.2\\.\\d+$", colnames(pheno)),
                          paste0("t3_", sub("\\.2\\.(\\d+)$", ".\\1", colnames(pheno))),
                          colnames(pheno))
colnames(pheno) <- ifelse(grepl("\\.3\\.\\d+$", colnames(pheno)),
                          paste0("t4_", sub("\\.3\\.(\\d+)$", ".\\1", colnames(pheno))),
                          colnames(pheno))

cat("Binarizing multiple answer phenotypes...", "\n")

#remove ICD10 diagnoses and cancer diagnoses because of high volume of <10,000 phenotypes per datafield - makes binarization unfeasable
patterns <- c("41202", "41204", "40006")
to_remove <- paste0(".*(", paste(patterns, collapse = "|"), ").*")
pheno <- pheno[, !grepl(to_remove, colnames(pheno))]

col_names <- colnames(pheno)
col_names <- col_names[col_names != "eid"]

#get unique phenotype names and suffixes
phenotype_unique <- str_extract(col_names, "^[^\\.]+")
sfx <- str_extract(col_names, "\\.[0-9]+$")

#filter for phenotypes that have more than the ".0" suffix, thereby indicating array structure
#saved in list structure for easy access
phenotype_dfs <- list()
for (phenotype in unique(phenotype_unique)) {
  phenotype_cols <- col_names[phenotype_unique == phenotype]
  extra_suffixes <- sum(!is.na(sfx[phenotype_unique == phenotype])) > 1
  if (extra_suffixes) {
    phenotype_dfs[[phenotype]] <- pheno %>%
      select(eid, all_of(phenotype_cols))  #add back eid to each df
  }
}

#define continuous datafields in the analysis. For these, the mean is derived across arrays instead of binarizing them.
cont_var <- c("4079", "4080", "20006", "20008", "102", "20009", "20132","20133", "20230")

###### calculate the mean across all cont columns for each participant. discard the original columns - THIS IS THE MEMORY HUNGRY PORTION OF THE CODE
###### if not continuous, convert the dataframe into wide format to binarize along answers given/diagnoses obtained. Retain NA columns carefully.
###### if parallel approach is used, use the appropriate version of the loop engaging foreach() and closing the cluster after the loop.
#phenotype_dfs <- foreach(name = names(phenotype_dfs), .combine = 'c', .packages = c("dplyr", "reshape2")) %dopar% {
 for (name in names(phenotype_dfs)) {
  cat("Processing dataframe:", name, "\n")
  #check if dataframe is continuous or categorcial
  if (any(grepl(paste(cont_var, collapse = "|"), name))) {
    cat("Identified as a continuous variable dataframe.\n")
    #if continuous derive mean and save column containing the mean, discard other columns
    mean_colname <- paste0(name, ".mean")
    phenotype_dfs[[name]] <- phenotype_dfs[[name]] %>%
      mutate(!!mean_colname := rowMeans(select(., -eid), na.rm = TRUE)) %>%
      select(eid, all_of(mean_colname))

  } else {
    cat("Identified as a categorical variable dataframe.\n")
    #if categorcial pull out participants that have only NAs, saved for adding back later
    na_data <- phenotype_dfs[[name]] %>%
     filter(rowSums(is.na(select(., -eid))) == ncol(select(., -eid))) %>%
     select(eid)
    #process non NA data
    non_na_data <- phenotype_dfs[[name]] %>%
      filter(rowSums(is.na(select(., -eid))) < ncol(select(., -eid)))
    #convert into long format by column and answer value. filter out NAs. Binarize answers
    long_data <- non_na_data %>%
      gather(key = "answer_column", value = "answer_value", -eid) %>%
      filter(!is.na(answer_value)) %>%
      mutate(answer_column = paste0(answer_column, ".", answer_value),
             answer_value = 1)
    #turn into wide format according to answer columns. Fill missing values with 0s to populate full df and indicate the negative.
    wide_data <- long_data %>%
      dcast(eid ~ answer_column, value.var = "answer_value", fill = 0)
    #create NA dummy data to add back NA participants, ensures that they have the same data structure and same N. Necessary to retain consistent row count across dataframes.
    na_dummy <- wide_data[0, ]
    na_dummy[1, ] <- NA
    na_dummy$eid <- NULL
    na_data <- na_data %>% bind_cols(na_dummy[rep(1, nrow(na_data)), ])
    #finally, bind wide and NA data
    phenotype_dfs[[name]] <- bind_rows(wide_data, na_data)

    cat("Finalized dataframe for:", name, "\n")
  }
  list(name = df)
  cat("Completed processing for dataframe:", name, "\n\n")
}
#stopCluster(cl)
##after stopping parallel cluster bind results back into one list
#phenotype_dfs <- setNames(phenotype_dfs, names(phenotype_dfs))

###### above output contains the binarized answers but specific to each array column. To save the data more efficiently, and make it interpretable, we collapse the dataframes into answer-columns across arrays
for (name in names(phenotype_dfs)) {
  #skip continuous variables, they only contain the mean
  if (!any(grepl(paste(cont_var, collapse = "|"), name))) {
    cat("Aggregating columns for categorical dataframe:", name, "\n")

    df <- phenotype_dfs[[name]]
    #extract possible responses from dataframe column names
    response_pattern <- paste0(name, "\\.[0-9]+\\.(-?[0-9]+)$")
    unique_responses <- unique(gsub(response_pattern, "\\1", colnames(df[-1])))
    unique_responses <- gsub("[^0-9-]", "", unique_responses)
    #create aggregated columns for each unique response
    aggregated_columns <- lapply(unique_responses, function(response) {
      #select columns matching the current response
      response_columns <- grep(paste0("\\.", response, "$"), colnames(df), value = TRUE)
      #aggregate columns, preserve NAs. Ensure that the only possible outcomes are 1, 0 or NA.
      aggregated_col <- apply(df[, response_columns, drop = FALSE], 1, function(row) {
        if (all(is.na(row))) {
          return(NA) #all values are NA
        } else if (any(row == 1, na.rm = TRUE)) {
          return(1)  #any value is 1
        } else {
          return(0)  #otherwise, 0
        }
      })
    return(aggregated_col)
    })

    #combine aggregated columns into new dataframe
    aggregated_df <- as.data.frame(aggregated_columns)
    colnames(aggregated_df) <- paste0(name, ".", unique_responses)

    #bind with eid column and overwrite old dataframe
    phenotype_dfs[[name]] <- cbind(eid = df$eid, aggregated_df)

    cat("Completed aggregation for:", name, "\n\n")
  }
}

#unwrap the df list into a single dataframe
pheno_unwrapped <- Reduce(function(x, y) merge(x, y, by = "eid", all = TRUE), phenotype_dfs)

#### bind the new binarized phenotypes back to the original dataframe
#remove original arrays from pheno
phenotype_remove <- unlist(lapply(names(phenotype_dfs), function(name) {
  grep(paste0("^", name, "\\."), colnames(pheno), value = TRUE)
}))
pheno <- pheno[, !colnames(pheno) %in% phenotype_remove]

#merge remaining pheno with unwrapped, binarized columns
pheno <- merge(pheno, pheno_unwrapped, by = "eid", all.x = TRUE)

#replace phenotype IDs with names using key
colnames(pheno) <- replacenumbers(colnames(pheno), phenoident, "phenotypes", "key")

##################################################################################################################################################
cat("Filtering phenotypes for sample size and excluding participants with neuroimaging data", "\n")
####remove all participants with imaging data from phenofile (they where potentially used in the creation of components and are therefore need to be removed)

imgid <- as.data.table(imgid)
pheno <- as.data.table(pheno)
pheno <- pheno[!eid %in% imgid$eid]

#filter pheno for GWAS sample sizes (50000 participants minimum per phenotype)
count <- pheno[, lapply(.SD, function(x) sum(!is.na(x)))]
cols <- names(count)[count > 50000]
pheno <- pheno[, ..cols]

#save for further processing
cat("Saving full non-neuroimaging phenotype data", "\n")
write.table(pheno, file.path(paste0(data_dir, "/cleanpheno_nonimage.txt")), sep="\t", row.names = FALSE, quote = FALSE)

# data subsetting
cat("Splitting data into categorical and continuous phenotypes for downstream analysis...", "\n")

#categorize pheno dataframe
upd <- data.frame(
  column_name = character(),
  unique_value_count = integer(),
  categorization = character(),
  stringsAsFactors = FALSE
)

for (col in names(pheno)) {
  #get unique count
  if (col == "eid") next
  unique_values <- unique(na.omit(pheno[[col]]))
  count_unique <- length(unique_values)

  #categorize
  if (col == "t1_Numeric_memory_test_-_Maximum_digits_remembered_correctly.0") { #exception made for this phenoype as it is continuous
    categorization <- "continuous"
  } else {
    categorization <- if (count_unique == 2) {
      "binary"
    } else if (count_unique >= 3 && count_unique <= 12) {
      "categorical"
    } else if (count_unique > 12) {
      "continuous"
    } else {
      NA
    }
  }

  #add to df
  upd <- rbind(
    upd,
    data.frame(
      column_name = col,
      unique_values = count_unique,
      categorization = categorization,
      stringsAsFactors = FALSE
    )
  )
}

#sort by categorization column
bin_cols <- upd %>%
  filter(categorization %in% c("binary", "categorical")) %>%
  pull(column_name)
cont_cols <- upd %>%
  filter(categorization %in% c("continuous", "interval")) %>%
  pull(column_name)

#ensure inclusion of eid
bin_cols <- c("eid", bin_cols)
cont_cols <- c("eid", cont_cols)

#subset pheno and convert to factors except continuous (NOTE: data.table)
pheno_bin <- pheno[, ..bin_cols]
pheno_bin[, (setdiff(names(pheno_bin), "eid")) := lapply(.SD, as.factor), .SDcols = setdiff(names(pheno_bin), "eid")]
pheno_cont <- pheno[, ..cont_cols]
pheno_cont[, (setdiff(names(pheno_cont), "eid")) := lapply(.SD, as.numeric), .SDcols = setdiff(names(pheno_cont), "eid")]

#save phenotype files
write.table(pheno_bin, file.path(paste0(outdir, "/cleanpheno_nonimage_bin.txt")), sep="\t", row.names = FALSE, quote = FALSE)
write.table(pheno_cont, file.path(paste0(outdir, "/cleanpheno_nonimage_cont.txt")), sep="\t", row.names = FALSE, quote = FALSE)

######################################################################
#######bind data into single df, may take a while due to data size. NOTE: This requires the prs.pca output and covar file from the previous script!
cat("Merging and saving phenotypes with component PCSs to create final datasets...", "\n")

##IC data
data1_bin <- prs_pca_oIC_nonimg_output %>%
  inner_join(covar, by = "eid") %>%
  inner_join(pheno_bin, by = "eid")
data1_cont <- prs_pca_oIC_nonimg_output %>%
  inner_join(covar, by = "eid") %>%
  inner_join(pheno_cont, by = "eid")

data1_bin <- as.data.table(data1_bin)
data1_bin[, (grep("^t1_", names(data1_bin), value = TRUE)) := lapply(.SD, as.factor), .SDcols = grep("^t1_", names(data1_bin), value = TRUE)]
data1_cont <- as.data.table(data1_cont)
data1_cont[, (grep("^t1_", names(data1_cont), value = TRUE)) := lapply(.SD, as.numeric), .SDcols = grep("^t1_", names(data1_cont), value = TRUE)]

write.table(data1_bin, file.path(paste0(outdir, "/cleandata_ICnonimage_bin.txt")), sep="\t", row.names = FALSE, quote = FALSE)
write.table(data1_cont, file.path(paste0(outdir, "/cleandata_ICnonimage_cont.txt")), sep="\t", row.names = FALSE, quote = FALSE)

##PC data

data2_bin <- prs_pca_pca_nonimg_output %>%
  inner_join(covar, by = "eid") %>%
  inner_join(pheno_bin, by = "eid")
data2_cont <- prs_pca_pca_nonimg_output %>%
  inner_join(covar, by = "eid") %>%
  inner_join(pheno_cont, by = "eid")

data2_bin <- as.data.table(data2_bin)
data2_bin[, (grep("^t1_", names(data2_bin), value = TRUE)) := lapply(.SD, as.factor), .SDcols = grep("^t1_", names(data2_bin), value = TRUE)]
data2_cont <- as.data.table(data2_cont)
data2_cont[, (grep("^t1_", names(data2_cont), value = TRUE)) := lapply(.SD, as.numeric), .SDcols = grep("^t1_", names(data2_cont), value = TRUE)]

write.table(data2_bin, file.path(paste0(outdir, "/cleandata_PCnonimage_bin.txt")), sep="\t", row.names = FALSE, quote = FALSE)
write.table(data2_cont, file.path(paste0(outdir, "/cleandata_PCnonimage_cont.txt")), sep="\t", row.names = FALSE, quote = FALSE)

cat("Data files created and saved.","\n")
cat("All done. Have a nice day.", "\n")
