# dependencies
library(tidyverse)
library(Biobase)
library(ComplexHeatmap)
library(circlize)

source(file.path("functions", "getColorScheme.R"))
source(file.path("functions", "removeNAsFromESet.R"))

# constants
useVersion   <- "publication"
colorScheme  <- getColorScheme()

# Set outputfolder
outputFolder <- file.path("figures", "output", "figure_SI01")

# create major output folder, if it doesn't exist
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder)

## load data for SEM replicates
sems <- readRDS(file = file.path("output", useVersion, "sems.RDS")) %>%
  removeNAsFromESet()

## get correlations
# sampleOrder <- c("Set1",
#                  "Set2",
#                  "Set3",
#                  "Set4",
#                  "Set4BR2",
#                  "Set5",
#                  "Set5BR2",
#                  "Set6",
#                  "Set7",
#                  "Set10")

sampleCorrelations <- sems %>%
  exprs() %>%
  `colnames<-`(colnames(.) %>%  gsub("_tmt10plex_", "_", .)) %>% 
  `colnames<-`(colnames(.) %>%  gsub("_130C", "_BR2", .)) %>%
  cor(method = "pearson")
  
 
sampleOrder <- order(as.numeric( gsub("Set|_.+.", "",colnames(sampleCorrelations))))
sampleCorrelations <- sampleCorrelations[sampleOrder, sampleOrder]

## heatmap
hm <- Heatmap(
  matrix = sampleCorrelations,
  col = colorRamp2(c(0, 1),
                   c(colorScheme$white,
                     colorScheme$blue)),
  cluster_rows = FALSE,
  cluster_columns = FALSE
)

print(hm)

pdf(file = file.path(outputFolder, "figure_SI01_sem_correlation_heatmap.pdf"), width = 5.7, height = 5)
hm
dev.off()
