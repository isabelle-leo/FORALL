####################################################################################
####################################################################################
# heatmap selected genes
####################################################################################
####################################################################################

# dependencies
library(Biobase)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(randomcoloR)
library(svglite)

sourceDir <- function(path, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    source(file.path(path, nm), ...)
  }
}
sourceDir(path = "functions")

# constants
name         <- "[Produce figure 1D] "
useVersion   <- "publication"
colorScheme  <- getColorScheme()

# set output folder
outputFolder <- file.path("figures", "output", "figure_01")
# create output folder, if it doesn't exist
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder)


####################################################################################
message(name, "Heatmap with all samples for selected proteins with replicates")
####################################################################################

# load data
proteins <- readRDS(file = file.path("output", useVersion, "proteins.RDS"))

# choose subtypes, cell lines you want to plot/keep from the eset.
proteins_filt <- proteins %>%
  .[, pData(.)$Subtype %in% c("BCR.ABL1",
                              "NUP214.ABL1",
                              "KMT2A.AFF1",
                              "KMT2A.MLLT1",
                              "TCF3.PBX1",
                              "ETV6.RUNX1",
                              "SPFQ.ABL1")] %>%
  .[fData(.)$gene_symbol %in% c("IGF2BP1", "PIK3C3", "RAG1", "HCK", "DSC3", "LMTK3", "SLC29A2", "MYC", "TAGLN2","ANKS1B", 
                                "RORB", "PHACTR3", "FAT1", "WNT16", "MERTK", "NLGN1", "ROR1","PBX1", "HMGA2","ABL1", 
                                "CAV1", "ABI2", "CDKN3", "BLM1", "LIMD1", "PECAM1", "CDH4","FLT3", "HOXA9", "RUNX2", 
                                "MEIS1", "DNTT", "DUSP1", "CD200", "SMAD1", "PROM1", "RYBP", "PBX3", "KLRK1","ABL1", 
                                "STAT3", "ULBP1", "TCEAL9", "PRSS57", "TFR2", "PTEN", "EPCAM", "KIT"), ]# %>%
 # removeNAsFromESet() #Not removed to include other markers of interest

# top annotation
ha_top <- HeatmapAnnotation(df = proteins_filt %>%
                              pData() %>%
                              dplyr::select(Type_Paper, Subtype_Paper),
                            col = list(Type_Paper = colorScheme$type,
                                       Subtype_Paper = colorScheme$subtype.alt),
                            show_annotation_name = TRUE,
                            annotation_name_gp = gpar(fontsize = 7),
                            show_legend = TRUE,
                            annotation_label = c("Type", "Subtype"))

# bottom annotation
ha_bottom <- HeatmapAnnotation(
  text = anno_text(x = pData(proteins_filt)$Cell_Line_Name_Paper,
                   rot = 90,
                   just = "right",
                   gp = gpar(fontsize = 7),
                   location = unit(2.2, "cm")))

# generate data matrix
mat <- exprs(proteins_filt)
mat <- remove_problematic_combs(mat, 1)

# the heatmap
hm <- Heatmap(name = "log2 protein level",
              matrix = mat,
              col = colorRamp2(c(-1.5, 0, 1.5),
                               c(colorScheme$blue,
                                 colorScheme$white,
                                 colorScheme$red)),
              na_col = "gray90",
              show_row_names = TRUE,
              row_names_gp = gpar(fontsize = 7),
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

pdf(file = file.path(outputFolder, "heatmap_selected_proteins.pdf"), width = 10)
hm
dev.off()

