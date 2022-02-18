# dependencies
library(Biobase)
library(tidyverse)
library(ggpubr)
library(ggthemes)

# dependencies heatmap
library(doMC)
library(ComplexHeatmap)
library(circlize)
library(randomcoloR)
library(networkD3)
library(viridis)

sourceDir <- function(path, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    source(file.path(path, nm), ...)
  }
}
sourceDir(path = "functions")

theme_set(theme_tufte())

# constants
useVersion   <- "publication"
colorScheme  <- getColorScheme()

#Choose output folder
outputFolder <- file.path("figures", "output", "figure_SI04")
#outputFolder <- file.path("output", useVersion, "figures", "figure_5")
# create major output folder, if it doesn't exist
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder)

####################### LOAD AND PREPARE DATA #######################
drugsens <- readRDS(file = file.path("data", "dsrt", "drugsens_2021-03-22.RDS"))

# filter away cell lines from eset
# drugsens_filt <- drugsens[, !pData(drugsens)$Cell.Line.R == "MN.60" & !pData(drugsens)$Cell.Line.R == "TANOUE"] %>%
#   removeNAsFromESet()

drugsens_filt <- drugsens

drugsens_filt <- drugsens_filt[fData(drugsens_filt)$FIMM.ID %in% c("FIMM002374",
                                                                   "FIMM000344",
                                                                   "FIMM003635",
                                                                   "FIMM000556",
                                                                   "FIMM136500",
                                                                   "FIMM133807",
                                                                   "FIMM003707",
                                                                   "FIMM133851",
                                                                   "FIMM023808",
                                                                   "FIMM136549",
                                                                   "FIMM003797",
                                                                   "FIMM136528",
                                                                   "FIMM023797",
                                                                   "FIMM003761",
                                                                   "FIMM002337",
                                                                   "FIMM136386",
                                                                   "FIMM023798",
                                                                   "FIMM000491",
                                                                   "FIMM000406",
                                                                   "FIMM003771",
                                                                   "FIMM136469",
                                                                   "FIMM100375",
                                                                   "FIMM003720",
                                                                   "FIMM136453",
                                                                   "FIMM136374"),]

# filter out B-ALL
# drugsens_filt <- drugsens_filt[, !pData(drugsens_filt)$Type == "B-ALL"] %>%
#  removeNAsFromESet()

# # filter out replicates, NAs and EBV
# drugsens_filt <- drugsens[, !pData(drugsens)$isReplicate & !pData(drugsens)$Type == "EBV"] %>%
#   removeNAsFromESet()

# drugsens_filt %>% pData()
# pData(drugsens_filt)

# top annotation
ha_top <- HeatmapAnnotation(df = drugsens_filt %>% pData() %>% dplyr::select(Type_Paper, hierarCluster_pearson_ward.D2_all, Subtype_Paper),
                            col = list(Type_Paper = colorScheme$type,
                                       hierarCluster_pearson_ward.D2_all = colorScheme$samplesClusters$hierar,
                                       Subtype_Paper = colorScheme$subtype.alt),
                            show_annotation_name = TRUE,
                            annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
                            show_legend = TRUE)

# bottom annotation
ha_bottom <- HeatmapAnnotation(
  text = anno_text(x = pData(drugsens_filt)$Cell.Line.R,
                   rot = 90,
                   just = "right",
                   gp = gpar(fontsize = 8, fontface = "bold"),
                   location = unit(1.7, "cm")))

# generate data matrix
mat <- exprs(drugsens_filt)


# the heatmap using viridis
hm <- Heatmap(name = "sDSS",
              matrix = mat,
              col = viridis(20, alpha = 0.1, option = "inferno"),
              show_row_names = TRUE,
              show_column_names = FALSE,
              cluster_rows = TRUE,
              clustering_distance_rows = "pearson",
              clustering_method_rows = "average",
              cluster_columns = TRUE,
              clustering_distance_columns = "pearson",
              clustering_method_columns = "average",
              show_heatmap_legend = TRUE,
              # column_split = factor(pData(proteins_filt)$hierarCluster_pearson_ward.D2, levels = c(1, 2, 4, 6, 3, 7, 5)),
              cluster_column_slices = FALSE,
              top_annotation = ha_top,
              bottom_annotation = ha_bottom)

# view heatmap
hm

# print heatmap to pdf
pdf(file = file.path(outputFolder, "figure_SI04_analysisDrugSens_Heatmap_Selected_Drugs_rev.pdf"), width = 10, height = 6)
hm
dev.off()

# # export to svg doesnt work for grids
# library(svglite)
# ggsave(
#   "figure_SI04_analysisDrugSens_Heatmap_Selected_Drugs_rev.svg",
#   plot = hm,
#   device = NULL,
#   path = outputFolder,
#   scale = 1,
#   width = 6,
#   height = 6,
#   units = c("in"),
#   dpi = 600,
#   limitsize = FALSE)

