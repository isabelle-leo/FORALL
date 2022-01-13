###################################################################################################
#####################----- finding highly variable protein coding genes -----######################
###################################################################################################
library(doMC)
library(Biobase)
library(mixtools)
library(dplyr)
outputFolder <- file.path("RNA_consensus_clustering")

RNA <- readRDS(file = "RNA.RDS")
# selecting protein coding RNAs
RNA <- RNA[fData(RNA)$gene_type == "protein_coding", ]
# log2 expression matrix
exprs(RNA) <- log2(exprs(RNA) + 0.01)
# remove any previous clustering
pData(RNA) <- pData(RNA)[, grep("hierar", colnames(pData(RNA)), invert = T)]


source("removeNAsFromESet.R")
source("detectHighVarianceExpressions_RNA.R")
types <- c("all", "B", "T")
reg_exp <- c("", "preB|B-ALL", "T-ALL") 

for( i in 1:length(types)){
  temp_eSet <- RNA[, grep(reg_exp[i], pData(RNA)$Type)]
  sd_row <- apply(exprs(temp_eSet), 1, sd)
  temp_eSet <- temp_eSet[sd_row > 0, ]
  
  RNA_out<- detectHighVarianceExpressions_RNA(eSet = temp_eSet,
                                              identifier = "gene_name",
                                              nCores = 2, maxit = 1000,
                                              outputFolder = file.path(outputFolder, types[i]))
  saveRDS(RNA_out, file = file.path(outputFolder,
                                    types[i],
                                    paste0("RNA.", types[i], ".RDS")))
}
###################################################################################################
################################----- RNA consensus clustering -----###############################
###################################################################################################
library(Biobase)
library(ConsensusClusterPlus)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(cluster)

name         <- "[RNA Consensus leukemic clustering] "
outputFolder <- file.path("RNA_consensus_clustering")

if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder, recursive = TRUE)
####################################################################################
message(name, "Load data")
####################################################################################
source("consensusClusterSamples.R")
for( i in c("all", "B", "T")){
  RNA <- readRDS(file = file.path( "RNA_consensus_clustering_MStahl_script", i, paste0("RNA.", i,".RDS")))
  
  # rename cluster columns, if already present
  colnames(pData(RNA)) <- gsub("hierarCluster", "previous_hierarCluster", colnames(pData(RNA)))
  
  ####################################################################################
  message(name, "Sample consensus clustering (without outliers)")
  ####################################################################################
  
  consClust <- consensusClusterSamples(eSet = RNA[fData(RNA)$hasHighVariance,                                                      !pData(RNA)$isReplicate ],
                                       maxK = 15,
                                       outputFolder = file.path(outputFolder, i, "consensus_clustering_samples"))
  saveRDS(consClust, file = file.path(outputFolder, i, paste(i, "consClust.RDS", sep = "_")))
}
###################################################################################################
###############----- using Silhouette method to find optimal number of clusters -----##############
###################################################################################################
outputFolder <- "RNA_consensus_clustering"
source("analyzeSilhouette.R")

for( i in c("all", "B", "T")){
  temp_eSet <- readRDS(paste(outputFolder, i, "/RNA.",i, ".RDS", sep = ""))
  analyzeSilhouettes(eSet = temp_eSet[fData(temp_eSet)$hasHighVariance,
                                      !pData(temp_eSet)$isReplicate ],
                     consClust = readRDS(file.path(outputFolder, i, paste(i, "_consClust.RDS", sep = ""))),
                     outputFolder = file.path(outputFolder, i, "silhouettes"))
}


# the following number of clusters are fine:
numClusters <- list("all" = 5, "B" = 6, "T" = 5) # by Silhouettes
message("Number of sample clusters used:")
message("all: ", numClusters$all)
message("B: ", numClusters$B)
message("T: ", numClusters$T)

method <- "pearson"
linkage <- "ward.D2"
outputFolder = "RNA_consensus_clustering"
# annotate ExpressionSets
library(Biobase)
library(ComplexHeatmap)
library(tidyverse)
library(circlize)
# 
source("clusterSamplesHierarchical_RNA.R")
source("removeNAsFromESet.R")

for( i in c("all", "B", "T")){
  RNA_updated_eSet <- clusterSamplesHierarchical(eSet =  readRDS(paste(outputFolder, i, "/RNA.", i, ".RDS", sep = "")),
                                                 k = numClusters[[i]],
                                                 method = method,
                                                 linkage = linkage,
                                                 heatmapAnnotation = "cell_line",
                                                 outputFolder = file.path(outputFolder, i, "CLCs"))
  saveRDS(RNA_updated_eSet, file = file.path(outputFolder, i, paste("/RNA.",i, ".RDS", sep = "")))
}
###################################################################################################