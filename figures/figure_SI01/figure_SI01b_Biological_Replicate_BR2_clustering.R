# dependencies
library(Biobase)
library(tidyverse)
library(doMC)
library(ComplexHeatmap)
library(circlize)
library(randomcoloR)
library(networkD3)
library(svglite)


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
outputFolder <- file.path("figures", "output", "figure_SI01")

# create major output folder, if it doesn't exist
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder)

####################################################################################
message(name, "Heatmap with all samples and replicates")
####################################################################################

# load data
proteins <- readRDS(file = file.path("output", useVersion, "consensus_leukemic_clustering", "all", "proteins.RDS"))

# filter out replicates and NAs
# proteins_filt <- proteins[, !pData(proteins)$isReplicate] %>%
#   removeNAsFromESet()

# filter out the SEM.BR.NOPS 
proteins_filt <- proteins[, !pData(proteins)$Cell.Line.R %in% c("SEM.BR.NOPS")]

proteins_filt <- proteins %>%
removeNAsFromESet()

proteins_filt %>% pData()
pData(proteins_filt)

# top annotation
ha_top <- HeatmapAnnotation(df = proteins_filt %>% pData() %>% dplyr::select(Type_Paper, hierarCluster_pearson_ward.D2, Subtype_Paper, Gender, Tissue),
                            col = list(Type_Paper = colorScheme$type,
                                       hierarCluster_pearson_ward.D2 = colorScheme$samplesClusters$hierar,
                                       Subtype_Paper = colorScheme$subtype.alt,
                                       Gender = colorScheme$gender,
                                       Tissue = colorScheme$tissue),
                            show_annotation_name = TRUE,
                            annotation_name_gp = gpar(fontsize = 7),
                            show_legend = TRUE)

# bottom annotation
ha_bottom <- HeatmapAnnotation(
  text = anno_text(x = pData(proteins_filt)$Cell.Line.R,
                   rot = 90,
                   just = "right",
                   gp = gpar(fontsize = 7),
                   location = unit(1.7, "cm")))

# generate data matrix
mat <- exprs(proteins_filt)

# the heatmap
hm <- Heatmap(name = "log2 expression",
              matrix = mat,
              col = colorRamp2(c(-1.5, 0, 1.5),
                               c(colorScheme$blue,
                                 colorScheme$white,
                                 colorScheme$red)),
              show_row_names = FALSE,
              show_column_names = FALSE,
              cluster_rows = TRUE,
              clustering_distance_rows = "pearson",
              clustering_method_rows = "ward.D2",
              cluster_columns = TRUE,
              clustering_distance_columns = "pearson",
              clustering_method_columns = "ward.D2",
              show_heatmap_legend = TRUE,
              # column_split = factor(pData(proteins_filt)$hierarCluster_pearson_ward.D2, levels = c(1, 2, 4, 6, 3, 7, 5)),
              cluster_column_slices = FALSE,
              top_annotation = ha_top,
              bottom_annotation = ha_bottom)

#view plot
hm

pdf(file = file.path(outputFolder, "figure_SI01_Biological_Replicate_BR2_clustering.pdf"), width = 10)
hm
dev.off()

# # export to svg
# library(svglite)
# ggsave(
#   "figure_SI01_Biological_Replicate_BR2_clustering.svg",
#   plot = hm,
#   device = NULL,
#   path = outputFolder,
#   scale = 1,
#   width = 6,
#   height = 6,
#   units = c("in"),
#   dpi = 600,
#   limitsize = FALSE)


