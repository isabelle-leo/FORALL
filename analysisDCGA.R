###################################################################################################
##################----- calls differential correlation DGCA pairwise pipeline -----################
###################################################################################################
# RNA protein correlation
rm(list = ls())
library(Biobase)
library(dplyr)
library(DGCA)
library(parallel)
compare_name = "RNA_protein_cor"
n_cores <- 4 # number of cores to use for parallelization (8 cores crashed , 4 cores took 6 hours)

protein_eSet <- readRDS("~/projects/ALL/luay_analysis/manuscript_analysis/Data_Analysis_2021/Data_Objects/proteins.RDS")
colnames(protein_eSet) <- protein_eSet$Cell.Line.R

RNA_eSet <- readRDS("~/projects/ALL/luay_analysis/manuscript_analysis/Data_Analysis_2021/Data_Objects/RNA.RDS")

# remove EBV
protein_eSet <- protein_eSet[, grep("EBV",  protein_eSet$Type, invert = T)]


# take common cell lines
common_samples <-  intersect(colnames(RNA_eSet), colnames(protein_eSet))
RNA_eSet <- RNA_eSet[, common_samples]
protein_eSet <- protein_eSet[, common_samples]

compare <- c("B", "T") 
Type <- ifelse(RNA_eSet$Type == "T-ALL", "T", "B") 

not_na_zero_sd <- function(mat  = "", two_groups_vec = Type){
  groups <- unique(two_groups_vec)
  apply(X = mat, MARGIN = 1, FUN = function(x){
    no_na <- sum(is.na(x)) == 0
    sd_1 <- sd(x[two_groups_vec == groups[1]], na.rm = T) != 0
    sd_2 <- sd(x[two_groups_vec == groups[2]], na.rm = T) != 0
    return(all(no_na, sd_1, sd_2))})
}


RNA_eSet <- RNA_eSet[not_na_zero_sd(mat = exprs(RNA_eSet), two_groups_vec = Type), ]

protein_eSet <- protein_eSet[not_na_zero_sd(mat = exprs(protein_eSet), two_groups_vec = Type), ]


design <-
  model.matrix(object = ~ Type + 0)   %>%
  `rownames<-` (colnames(RNA_eSet)) %>%
  `colnames<-`(compare)

n_features <- nrow(protein_eSet)
split_indices <- split(x = 1:n_features, f =  ceiling(seq_along(1:n_features)/ceiling(n_features/n_cores)))


for ( i in  c("pearson", "spearman")){
  ddcorAll_output <- mclapply(X = split_indices, 
                              FUN = function(x){
                                temp_ddcorAll_output <-  ddcorAll(inputMat = exprs(RNA_eSet)[x, ],
                                                                  design = design,
                                                                  compare = compare,
                                                                  inputMatB = exprs(protein_eSet) ,
                                                                  splitSet = NULL,
                                                                  impute = FALSE,
                                                                  corrType = i,
                                                                  nPairs = "all",
                                                                  sortBy = "zScoreDiff",
                                                                  adjust = "perm",
                                                                  nPerms = 10,
                                                                  classify = TRUE,
                                                                  sigThresh = 1,
                                                                  corSigThresh = 0.05,
                                                                  heatmapPlot = FALSE,
                                                                  color_palette = NULL,
                                                                  verbose = FALSE,
                                                                  plotFdr = FALSE,
                                                                  corr_cutoff = 0.99,
                                                                  signType = "none",
                                                                  getDCorAvg = FALSE,
                                                                  dCorAvgType = "gene_average",
                                                                  dCorAvgMethod = "median",
                                                                  oneSidedPVal = FALSE,
                                                                  customize_heatmap = FALSE,
                                                                  heatmapClassic = FALSE,
                                                                  corPower = 2)
                                return(temp_ddcorAll_output)
                              }, mc.cores = n_cores)
  ddcorAll_output <- do.call("rbind", ddcorAll_output)
  
  colnames(ddcorAll_output)[1:2] <- c("RNA", "protein")
  ddcorAll_output[, 3:10] <- signif(ddcorAll_output[, 3:10], 3)
  
  saveRDS(object = ddcorAll_output, file = paste(compare_name, "_B_vs_T_", i, "_ddcorAll_output.RDS", sep = ""))
  
  ddcorAll_output <- ddcorAll_output[ddcorAll_output$pValDiff < 0.05, ]
  write.table(x =  ddcorAll_output, file = paste(compare_name, "_B_vs_T_", i, "_ddcorAll_output_pValDiff_0.05.txt", sep = ""), sep = "\t", row.names = F, quote = F)
}

# RNA sDSS correlation
rm(list = ls())
library(Biobase)
library(dplyr)
library(DGCA)
library(parallel)
compare_name = "RNA_sDSS_cor"

sDSS_eSet <- readRDS("~/projects/ALL/luay_analysis/manuscript_analysis/Data_Analysis_2021/Data_Objects/sDSS.RDS")

RNA_eSet <- readRDS("~/projects/ALL/luay_analysis/manuscript_analysis/Data_Analysis_2021/Data_Objects/RNA.RDS")

# remove EBV and B-ALL
sDSS_eSet <- sDSS_eSet[, grep("EBV|B-ALL",  sDSS_eSet$Type, invert = T)]

# take common cell lines
common_samples <-  intersect(colnames(RNA_eSet), colnames(sDSS_eSet))
RNA_eSet <- RNA_eSet[, common_samples]
sDSS_eSet <- sDSS_eSet[, common_samples]

compare <- c("B", "T") 
Type <- ifelse(RNA_eSet$Type == "T-ALL", "T", "B") 

not_na_zero_sd <- function(mat  = "", two_groups_vec = Type){
  groups <- unique(two_groups_vec)
  apply(X = mat, MARGIN = 1, FUN = function(x){
    no_na <- sum(is.na(x)) == 0
    sd_1 <- sd(x[two_groups_vec == groups[1]], na.rm = T) != 0
    sd_2 <- sd(x[two_groups_vec == groups[2]], na.rm = T) != 0
    return(all(no_na, sd_1, sd_2))})
}


RNA_eSet <- RNA_eSet[not_na_zero_sd(mat = exprs(RNA_eSet), two_groups_vec = Type), ]

sDSS_eSet <- sDSS_eSet[not_na_zero_sd(mat = exprs(sDSS_eSet), two_groups_vec = Type), ]

design <-
  model.matrix(object = ~ Type + 0)   %>%
  `rownames<-` (colnames(RNA_eSet)) %>%
  `colnames<-`(compare)

for ( i in  c("pearson", "spearman")){
  
  ddcorAll_output <-  ddcorAll(inputMat = exprs(RNA_eSet),
                               design = design,
                               compare = compare,
                               inputMatB = exprs(sDSS_eSet) ,
                               splitSet = NULL,
                               impute = FALSE,
                               corrType = i,
                               nPairs = "all",
                               sortBy = "zScoreDiff",
                               adjust = "perm",
                               nPerms = 10,
                               classify = TRUE,
                               sigThresh = 1,
                               corSigThresh = 0.05,
                               heatmapPlot = FALSE,
                               color_palette = NULL,
                               verbose = FALSE,
                               plotFdr = FALSE,
                               corr_cutoff = 0.99,
                               signType = "none",
                               getDCorAvg = FALSE,
                               dCorAvgType = "gene_average",
                               dCorAvgMethod = "median",
                               oneSidedPVal = FALSE,
                               customize_heatmap = FALSE,
                               heatmapClassic = FALSE,
                               corPower = 2)
  
  
  
  colnames(ddcorAll_output)[1:2] <- c("RNA", "drug")
  ddcorAll_output[, 3:10] <- signif(ddcorAll_output[, 3:10], 3)
  
  ddcorAll_output <- cbind(ddcorAll_output, fData(sDSS_eSet)[ddcorAll_output$drug, c("Mechanism.Targets", "Putative.Target.Protein", "Drug_class") ])
  ddcorAll_output$Drug_name <- gsub(".+.\\.\\.\\.", "", ddcorAll_output$drug)
  
  saveRDS(object = ddcorAll_output, file = paste(compare_name, "_B_vs_T_", i, "_ddcorAll_output.RDS", sep = ""))
  
  ddcorAll_output <- ddcorAll_output[ddcorAll_output$pValDiff < 0.05, ]
  write.table(x =  ddcorAll_output, file = paste(compare_name, "_B_vs_T_", i, "_ddcorAll_output_pValDiff_0.05.txt", sep = ""), sep = "\t", row.names = F, quote = F)
}

# protein sDSS correlation
rm(list = ls())
library(Biobase)
library(dplyr)
library(DGCA)
library(parallel)
compare_name = "protein_sDSS_cor"

sDSS_eSet <- readRDS("~/projects/ALL/luay_analysis/manuscript_analysis/Data_Analysis_2021/Data_Objects/sDSS.RDS")

protein_eSet <- readRDS("~/projects/ALL/luay_analysis/manuscript_analysis/Data_Analysis_2021/Data_Objects/proteins.RDS")
colnames(protein_eSet) <- protein_eSet$Cell.Line.R
# remove EBV and B-ALL
sDSS_eSet <- sDSS_eSet[, grep("EBV|B-ALL",  sDSS_eSet$Type, invert = T)]

# take common cell lines
common_samples <-  intersect(colnames(protein_eSet), colnames(sDSS_eSet))
protein_eSet <- protein_eSet[, common_samples]
sDSS_eSet <- sDSS_eSet[, common_samples]

compare <- c("B", "T") 
Type <- ifelse(protein_eSet$Type == "T-ALL", "T", "B") 

not_na_zero_sd <- function(mat  = "", two_groups_vec = Type){
  groups <- unique(two_groups_vec)
  apply(X = mat, MARGIN = 1, FUN = function(x){
    no_na <- sum(is.na(x)) == 0
    sd_1 <- sd(x[two_groups_vec == groups[1]], na.rm = T) != 0
    sd_2 <- sd(x[two_groups_vec == groups[2]], na.rm = T) != 0
    return(all(no_na, sd_1, sd_2))})
}


protein_eSet <- protein_eSet[not_na_zero_sd(mat = exprs(protein_eSet), two_groups_vec = Type), ]

sDSS_eSet <- sDSS_eSet[not_na_zero_sd(mat = exprs(sDSS_eSet), two_groups_vec = Type), ]

design <-
  model.matrix(object = ~ Type + 0)   %>%
  `rownames<-` (colnames(protein_eSet)) %>%
  `colnames<-`(compare)

for ( i in  c("pearson", "spearman")){
  
  ddcorAll_output <-  ddcorAll(inputMat = exprs(protein_eSet),
                               design = design,
                               compare = compare,
                               inputMatB = exprs(sDSS_eSet) ,
                               splitSet = NULL,
                               impute = FALSE,
                               corrType = i,
                               nPairs = "all",
                               sortBy = "zScoreDiff",
                               adjust = "perm",
                               nPerms = 10,
                               classify = TRUE,
                               sigThresh = 1,
                               corSigThresh = 0.05,
                               heatmapPlot = FALSE,
                               color_palette = NULL,
                               verbose = FALSE,
                               plotFdr = FALSE,
                               corr_cutoff = 0.99,
                               signType = "none",
                               getDCorAvg = FALSE,
                               dCorAvgType = "gene_average",
                               dCorAvgMethod = "median",
                               oneSidedPVal = FALSE,
                               customize_heatmap = FALSE,
                               heatmapClassic = FALSE,
                               corPower = 2)
  
  
  colnames(ddcorAll_output)[1:2] <- c("protein", "drug")
  ddcorAll_output[, 3:10] <- signif(ddcorAll_output[, 3:10], 3)
  
  ddcorAll_output <- cbind(ddcorAll_output, fData(sDSS_eSet)[ddcorAll_output$drug,  c("Mechanism.Targets", "Putative.Target.Protein", "Drug_class")])
  ddcorAll_output$Drug_name <- gsub(".+.\\.\\.\\.", "", ddcorAll_output$drug)
  
  saveRDS(object = ddcorAll_output, file = paste(compare_name, "_B_vs_T_", i, "_ddcorAll_output.RDS", sep = ""))
  
  ddcorAll_output <- ddcorAll_output[ddcorAll_output$pValDiff < 0.05, ]
  write.table(x =  ddcorAll_output, file = paste(compare_name, "_B_vs_T_", i, "_ddcorAll_output_pValDiff_0.05.txt", sep = ""), sep = "\t", row.names = F, quote = F)
}

###################################################################################################
###################################################################################################
###################################################################################################
# extract corresponding RNA-protein pairs
rm(list = ls())
setwd("~/projects/ALL/luay_analysis/manuscript_analysis/Data_Analysis_2021/DGCA_analysis/")
for( i in c("pearson", "spearman")){
  RNA_protein <- readRDS(paste("RNA_protein_cor_B_vs_T", i, "ddcorAll_output.RDS", sep = "_"))
  keep <- RNA_protein$RNA == RNA_protein$protein
  RNA_protein <-  RNA_protein[keep, ]
  RNA_protein <- RNA_protein[order(RNA_protein$pValDiff), ]
  write.table(RNA_protein, file = paste("matched_RNA_protein_cor_B_vs_T", i , "ddcorAll_output.txt", sep = "\t"), sep = "\t", row.names = F, quote = F)
}
###################################################################################################
###################################################################################################
