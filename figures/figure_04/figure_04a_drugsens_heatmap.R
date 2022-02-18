# dependencies
library(Biobase)
library(tidyverse)
library(ggpubr)
library(ggthemes)
library(svglite)

# dependencies heatmap
library(ComplexHeatmap)
library(circlize)
library(randomcoloR)
library(viridis)

source(file.path("functions", "getColorScheme.R"))
source(file.path("figures", "all_theme.R"))
theme_set(all_theme())

# constants
useVersion   <- "publication"
colorScheme  <- getColorScheme()

#Choose outputfolder
outputFolder <- file.path("figures", "output", "figure_04")

# create major output folder, if it doesn't exist
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder)


####################### LOAD AND PREPARE DATA #######################
drugsens <- readRDS(file = file.path("data", "dsrt", "drugsens_2021-03-22.RDS")) 

# filter the drugs with zero change
drugsens_filt <- drugsens %>%
  .[, !pData(.)$Cell.Line.R %in% c("SEM.BR.NOPS")] %>% #"TANOUE", "TANOUE.BR2", "MN.60", "MN.60.BR2", "SEM.BR.NOPS"
  .[apply(exprs(.), 1, sd) != 0, ]

# drugsens_filt %>% pData()
# pData(drugsens_filt)

# top annotation
ha_top <- HeatmapAnnotation(df = drugsens_filt %>%
                              pData() %>%
                              dplyr::select(Type_Paper, hierarCluster_pearson_ward.D2_all, Subtype_Paper),
                            col = list(Type_Paper = colorScheme$type,
                                       hierarCluster_pearson_ward.D2_all = colorScheme$samplesClusters$hierar,
                                       Subtype_Paper = colorScheme$subtype.alt),
                            show_annotation_name = TRUE,
                            annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
                            show_legend = TRUE)

# bottom annotation
ha_bottom <- HeatmapAnnotation(
  text = anno_text(x = pData(drugsens_filt)$Cell_Line_Name_Paper,
                   rot = 90,
                   just = "right",
                   gp = gpar(fontsize = 10, fontface = "bold"),
                   location = unit(1.7, "cm")))

# generate data matrix
mat <- exprs(drugsens_filt)

# the heatmap using viridis
hm <- Heatmap(name = "sDSS",
              matrix = mat,
              col = viridis(20, alpha = 0.1, option = "inferno"),
              show_row_names = FALSE,
              show_column_names = FALSE,
              cluster_rows = TRUE,
              clustering_distance_rows = "pearson",
              clustering_method_rows = "complete",
              cluster_columns = TRUE,
              clustering_distance_columns = "pearson",
              clustering_method_columns = "complete",
              show_heatmap_legend = TRUE,
              cluster_column_slices = FALSE,
              top_annotation = ha_top,
              bottom_annotation = ha_bottom)

# view heatmap
hm

# print heatmap to pdf (remember to change file name to match fig)
pdf(file = file.path(outputFolder, "figure_04_drugsens_heatmap_rev.pdf"), width = 8, height = 10)
hm
dev.off()
