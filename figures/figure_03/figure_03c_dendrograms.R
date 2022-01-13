# dependencies
library(Biobase)
library(tidyverse)
library(doMC)
library(networkD3)

sourceDir <- function(path, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    source(file.path(path, nm), ...)
  }
}
sourceDir(path = "functions")

# constants
useVersion   <- "publication"
colorScheme  <- getColorScheme()

# Set outputfolder
outputFolder <- file.path("figures", "output", "figure_03")

# create maijor output folder, if it doesn't exist
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder)

# load data
proteins_B <- readRDS(file = file.path("output", useVersion, "consensus_leukemic_clustering", "B", "proteins.B.RDS")) %>%
  .[, pData(.)$Cell.Line.R != "SEM.BR.NOPS"] %>%
  .[, !pData(.)$isReplicate]
proteins_T <- readRDS(file = file.path("output", useVersion, "consensus_leukemic_clustering", "T", "proteins.T.RDS")) %>%
  .[, !pData(.)$isReplicate]


# B dendrograms
dendrogram2D(proteins = proteins_B,
             nodeSize = 600,
             colorScheme = colorScheme,
             outputFolder = file.path(outputFolder, "dendrogram_B"))

# T dendrograms
dendrogram2D(proteins = proteins_T,
             nodeSize = 600,
             colorScheme = colorScheme,
             outputFolder = file.path(outputFolder, "dendrogram_T"))
