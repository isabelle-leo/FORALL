# dependencies
library(ggplot2)
library(ggrepel)
library(ggthemes)
library(tidyverse)

# dependencies
sourceDir <- function(path, ...) {
        for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
                source(file.path(path, nm), ...)
        }
}
sourceDir(path = "functions")

source(file.path("functions", "getColorScheme.R"))
source(file.path("figures", "all_theme.R"))
theme_set(all_theme())

# constants
useVersion   <- "publication"
colorScheme  <- getColorScheme()

# Set output folder
outputFolder <- file.path("figures", "output", "figure_SI01")
# create major output folder, if it doesn't exist
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder)

# Set genes to plot/highlight on the scatter plot from file
genes_of_interest <- read_tsv(file = file.path("output", useVersion, "Luay", "StJude", "genes_of_interest_stjude_vs_proteins_ETV6-RUNX1.txt"), col_names = F) %>% .[[1]]

x_DEG_output <- file.path("output", useVersion, "Luay", "StJude", "DEqMS_output_ETV6-RUNX1_vs_BCP-ALL_2021_08_19_10_08_59.txt")
y_DEG_output <- file.path("output", useVersion, "Luay", "StJude", "edgeR_output_ETV6-RUNX1_vs_Control_2021_08_19_10_03_57.txt")
title_lab = "figure_SI01_scatterplot_ETV6-RUNX1_DE_stjude_common_deqms_rev"
FC_scatter_plot(x_DEG_output = x_DEG_output,
                y_DEG_output = y_DEG_output,x_lab = "ETV6-RUNX1 cell lines protein FC", y_lab = "ETV6-RUNX1 clinical RNA FC", title_lab = title_lab,
                filter_by = "pvalue", threshold = 0.01, points_color =  "#CDCC03", fitting_line_color = '#2C3E50', points_size = 2, genes_of_interest = genes_of_interest, labels_size = 4,
                labels_color = "#555555" , labeld_points_color =  "#35648A")