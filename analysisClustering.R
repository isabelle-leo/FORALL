####################################################################################
####################################################################################
# consensus leukemic cluster detection
####################################################################################
####################################################################################

# dependencies
library(Biobase)
library(doMC)
library(mixtools)
library(ConsensusClusterPlus)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(cluster)

sourceDir <- function(path, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    source(file.path(path, nm), ...)
  }
}
sourceDir(path = "functions")

# constants
name         <- "[Consensus leukemic clustering] "
useVersion   <- "publication"
outputFolder <- file.path("output", useVersion, "consensus_leukemic_clustering")
colorScheme  <- getColorScheme()

if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder, recursive = TRUE)


####################################################################################
message(name, "Load data")
####################################################################################

proteins <- readRDS(file = file.path("output", useVersion, "proteins.RDS"))

# rename cluster columns, if already present
colnames(pData(proteins)) <- gsub("^hierarCluster", "previous_hierarCluster", colnames(pData(proteins)))


####################################################################################
message(name, "Read already processed data")
####################################################################################

# comment these lines out if you want to start from the beginning
# proteins <- readRDS(file = file.path("output", useVersion, "consensus_leukemic_clustering", "all", "proteins.RDS"))
# proteins.B <- readRDS(file = file.path("output", useVersion, "consensus_leukemic_clustering", "B", "proteins.B.RDS"))
# proteins.T <- readRDS(file = file.path("output", useVersion, "consensus_leukemic_clustering", "T", "proteins.T.RDS"))


####################################################################################
message(name, "Make sample subsets")
####################################################################################

proteins.B <- proteins[, pData(proteins)$Type_Paper %in% c("BCP-ALL", "B-ALL")]
proteins.T <- proteins[, pData(proteins)$Type_Paper %in% c("T-ALL")]


####################################################################################
message(name, "Identify high-variance proteins and sample outliers")
####################################################################################

# filter for high-variance proteins
proteins <- detectHighVarianceExpressions(eSet = proteins,
                                          identifier = "gene_symbol",
                                          nCores = 2,
                                          maxit = 2000,
                                          outputFolder = file.path(outputFolder, "all", "high_variance"))

proteins.B <- detectHighVarianceExpressions(eSet = proteins.B,
                                            identifier = "gene_symbol",
                                            nCores = 2,
                                            maxit = 2000,
                                            outputFolder = file.path(outputFolder, "B", "high_variance"))

proteins.T <- detectHighVarianceExpressions(eSet = proteins.T,
                                            identifier = "gene_symbol",
                                            nCores = 2,
                                            maxit = 2000,
                                            outputFolder = file.path(outputFolder, "T", "high_variance"))

# get correlated samples (p >= 0.5)
proteins <- detectSampleOutliers(eSet = proteins,
                                 pearsonCutOff = .5,
                                 outputFolder = file.path(outputFolder, "all", "sample_non_outliers"))

proteins.B <- detectSampleOutliers(eSet = proteins.B,
                                   pearsonCutOff = .5,
                                   outputFolder = file.path(outputFolder, "B", "sample_non_outliers"))

proteins.T <- detectSampleOutliers(eSet = proteins.T,
                                   pearsonCutOff = .5,
                                   outputFolder = file.path(outputFolder, "T", "sample_non_outliers"))


####################################################################################
message(name, "Sample consensus clustering")
####################################################################################

consClust <- consensusClusterSamples(eSet = proteins[fData(proteins)$hasHighVariance,
                                                     !pData(proteins)$isReplicate & !pData(proteins)$isOutlier],
                                     maxK = 15,
                                     outputFolder = file.path(outputFolder, "all", "consensus_clustering_samples"))

consClust.B <- consensusClusterSamples(eSet = proteins.B[fData(proteins.B)$hasHighVariance,
                                                         !pData(proteins.B)$isReplicate & !pData(proteins.B)$isOutlier],
                                       maxK = 12,
                                       outputFolder = file.path(outputFolder, "B", "consensus_clustering_samples"))

consClust.T <- consensusClusterSamples(eSet = proteins.T[fData(proteins.T)$hasHighVariance,
                                                         !pData(proteins.T)$isReplicate & !pData(proteins.T)$isOutlier],
                                       maxK = 12,
                                       outputFolder = file.path(outputFolder, "T", "consensus_clustering_samples"))

# with sample outliers
consClust <- consensusClusterSamples(eSet = proteins[fData(proteins)$hasHighVariance,
                                                     !pData(proteins)$isReplicate],
                                     maxK = 15,
                                     outputFolder = file.path(outputFolder, "all", "consensus_clustering_samples_with_sample_outliers"))

consClust.B <- consensusClusterSamples(eSet = proteins.B[fData(proteins.B)$hasHighVariance,
                                                         !pData(proteins.B)$isReplicate],
                                       maxK = 12,
                                       outputFolder = file.path(outputFolder, "B", "consensus_clustering_samples_with_sample_outliers"))

consClust.T <- consensusClusterSamples(eSet = proteins.T[fData(proteins.T)$hasHighVariance,
                                                         !pData(proteins.T)$isReplicate],
                                       maxK = 12,
                                       outputFolder = file.path(outputFolder, "T", "consensus_clustering_samples_with_sample_outliers"))

# Silhouette plots
analyzeSilhouettes(eSet = proteins[fData(proteins)$hasHighVariance,
                                   !pData(proteins)$isReplicate & !pData(proteins)$isOutlier],
                   consClust = consClust,
                   outputFolder = file.path(outputFolder, "all", "silhouettes"))

analyzeSilhouettes(eSet = proteins.B[fData(proteins.B)$hasHighVariance,
                                     !pData(proteins.B)$isReplicate & !pData(proteins.B)$isOutlier],
                   consClust = consClust.B,
                   outputFolder = file.path(outputFolder, "B", "silhouettes"))

analyzeSilhouettes(eSet = proteins.T[fData(proteins.T)$hasHighVariance,
                                     !pData(proteins.T)$isReplicate & !pData(proteins.T)$isOutlier],
                   consClust = consClust.T,
                   outputFolder = file.path(outputFolder, "T", "silhouettes"))

# the following number of clusters are fine:
numClusters <- list("all" = 7, "B" = 5, "T" = 7)
message("Number of sample clusters used:")
message("all: ", numClusters$all)
message("B: ", numClusters$B)
message("T: ", numClusters$T)


####################################################################################
message(name, "Cluster samples to core leukemic clusters (CLCs)")
####################################################################################

proteins <- clusterSamplesHierarchical(eSet = proteins,
                                       k = numClusters$all,
                                       method = "pearson",
                                       linkage = "ward.D2",
                                       heatmapAnnotation = "cell_line",
                                       outputFolder = file.path(outputFolder, "all", "CLCs"))

proteins.B <- clusterSamplesHierarchical(eSet = proteins.B,
                                         k = numClusters$B,
                                         method = "pearson",
                                         linkage = "ward.D2",
                                         heatmapAnnotation = "cell_line",
                                         outputFolder = file.path(outputFolder, "B", "CLCs"))

proteins.T <- clusterSamplesHierarchical(eSet = proteins.T,
                                         k = numClusters$T,
                                         method = "pearson",
                                         linkage = "ward.D2",
                                         heatmapAnnotation = "cell_line",
                                         outputFolder = file.path(outputFolder, "T", "CLCs"))


####################################################################################
message(name, "Calculate median expressions for CLCs")
####################################################################################

clcString <- "__clc_median_expression_"
proteins <- addMedianExpressions(eSet = proteins,
                                 use = "hierarCluster_pearson_ward.D2",
                                 recognitionString = clcString)

proteins.B <- addMedianExpressions(eSet = proteins.B,
                                   use = "hierarCluster_pearson_ward.D2",
                                   recognitionString = clcString)

proteins.T <- addMedianExpressions(eSet = proteins.T,
                                   use = "hierarCluster_pearson_ward.D2",
                                   recognitionString = clcString)


####################################################################################
message(name, "Add type medians")
####################################################################################

typeString <- "__type_median_expression_"
proteins <- addMedianExpressions(eSet = proteins,
                                 use = "Type_Paper",
                                 recognitionString = typeString)

proteins.B <- addMedianExpressions(eSet = proteins.B,
                                   use = "Type_Paper",
                                   recognitionString = typeString)

proteins.T <- addMedianExpressions(eSet = proteins.T,
                                   use = "Type_Paper",
                                   recognitionString = typeString)


####################################################################################
message(name, "Cluster proteins")
####################################################################################

# proteins <- clusterFeatures(eSet = proteins,
#                             method = "pearson",
#                             linkage = "average")
# 
# proteins.B <- clusterFeatures(eSet = proteins.B,
#                               method = "pearson",
#                               linkage = "average")
# 
# proteins.T <- clusterFeatures(eSet = proteins.T,
#                               method = "pearson",
#                               linkage = "average")


####################################################################################
message(name, "Save (modified) objects")
####################################################################################

saveRDS(object = proteins,
        file = file.path(outputFolder, "all", "proteins.RDS"))

saveRDS(object = proteins.B,
        file = file.path(outputFolder, "B", "proteins.B.RDS"))

saveRDS(object = proteins.T,
        file = file.path(outputFolder, "T", "proteins.T.RDS"))


####################################################################################
message(name, "Save txt tables")
####################################################################################

saveESetAsTxt(eSet = proteins,
              file = file.path(outputFolder, "all", "proteins.txt"),
              colNameColumn = "Cell.Line.R",
              rowNameColumn = "gene_symbol",
              additionalRows = c("Subtype", "Gender", "hierarCluster_pearson_ward.D2"),
              extendedColnames = FALSE)

saveESetAsTxt(eSet = proteins.B,
              file = file.path(outputFolder, "B", "proteins.B.txt"),
              colNameColumn = "Cell.Line.R",
              rowNameColumn = "gene_symbol",
              additionalRows = c("Subtype", "Gender", "hierarCluster_pearson_ward.D2"),
              extendedColnames = FALSE)

saveESetAsTxt(eSet = proteins.T,
              file = file.path(outputFolder, "T", "proteins.T.txt"),
              colNameColumn = "Cell.Line.R",
              rowNameColumn = "gene_symbol",
              additionalRows = c("Subtype", "Gender", "hierarCluster_pearson_ward.D2"),
              extendedColnames = FALSE)


####################################################################################
message(name, "Store the clustering")
####################################################################################

# full dataset
meta <- proteins %>%
  pData() %>%
  dplyr::select(starts_with("hierarCluster"), Cell.Line.R)

# B dataset
meta_B <- proteins.B %>%
  pData() %>%
  dplyr::select(starts_with("hierarCluster"), Cell.Line.R)

# T dataset
meta_T <- proteins.T %>%
  pData() %>%
  dplyr::select(starts_with("hierarCluster"), Cell.Line.R)

# fuse clustering data
metaExport <- meta %>%
  dplyr::select(Cell.Line.R, hierarCluster_pearson_ward.D2_all = hierarCluster_pearson_ward.D2) %>%
  left_join(meta_B %>% dplyr::select(Cell.Line.R, hierarCluster_pearson_ward.D2_B = hierarCluster_pearson_ward.D2),
            by = "Cell.Line.R") %>%
  left_join(meta_T %>% dplyr::select(Cell.Line.R, hierarCluster_pearson_ward.D2_T = hierarCluster_pearson_ward.D2),
            by = "Cell.Line.R")

write_tsv(x = metaExport,
          file = file.path("meta", "cell_line_clustering.txt"))


####################################################################################
message(name, "Clear environment")
####################################################################################

rm(list = ls())