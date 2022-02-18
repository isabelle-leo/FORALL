## Libraries
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

# Set outputfolder
outputFolder <- file.path("figures", "output", "figure_SI06")
# create maijor output folder, if it doesn't exist
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder)

# Set gens to plot/highlight on the scatterplot from file
genes_of_interest <- read_tsv(file = file.path("output", useVersion, "Luay", "StJude", "genes_of_interest_stjude_vs_proteins.txt"), col_names = F) %>% .[[1]]

# Set gens to plot/highlight on the scatterplot (below from volcano plot)
# genes_of_interest <- c("HDAC9", "NFATC1", "LCK", "MEF2D", "MEF2C", "IDH", "NRIP1", "AMER2", "PTPRZ1", "CKM", "CMTM3", 
#                  "ADSS1", "TCF4", "MAGI1", "PLCL2", "CAMK4", "ARHGAP6", "PTK2", "FBXO27", "FRMPD1", "BCL11A", "CD5", "DNTT", "RRAS", "KRAS", "NRAS")

####################################################################################################################################
# Filter genes on plot by pVal ALL MEF2D vs all non-MEF2D (excluding all not identified subtype merged stranded & unstranded)
####################################################################################################################################

# cell_lines_protein_vs_stjude_RNA_FC edgeR_output_MEF2D_vs_Rest_exclude_not_identified_Stranded_unstranded_2021_03_12_12_03_05
x_DEG_output <- file.path("output", useVersion, "Luay", "StJude", "MEF2D.HNRNPUL1-BCPALL_DEqMS_table_common.txt")
y_DEG_output <- file.path("output", useVersion, "Luay", "StJude", "edgeR_output_MEF2D_vs_Rest_exclude_not_identified_Stranded_unstranded_2021_03_12_12_03_05.txt")
title_lab = "Supp_Figure_5A_cell_lines_protein_vs_clinical_RNA_FC_MEF2D"
FC_scatter_plot(x_DEG_output = x_DEG_output,
                y_DEG_output = y_DEG_output,x_lab = "MEF2D-HNRNPUL1 cell lines protein FC", y_lab = "MEF2D-rearranged clinical RNA FC", title_lab = title_lab,
                filter_by = "pvalue", threshold = 0.01, points_color =  "#CDCC03", fitting_line_color = '#2C3E50', points_size = 2, genes_of_interest = genes_of_interest, labels_size = 4,
                labels_color = "#555555" , labeld_points_color =  "#35648A")

####################################################################################################################################
# Filter genes on plot by pVal Gu et al 2016 Nat Comm
####################################################################################################################################
#cell_lines_protein_vs_stjude_RNA_FC MEF2D-BCL9
x_DEG_output <- file.path("output", useVersion, "Luay", "StJude", "MEF2D.HNRNPUL1-BCPALL_DEqMS_table_common.txt")
y_DEG_output <- file.path("output", useVersion, "Luay", "StJude", "stjude_gu_ncomm_output_MEF2D_combined_no_dups.txt")
title_lab <- "Supp_Figure_5A_cell_lines_protein_vs_stjude_Gu_et_al_RNA_FC"
FC_scatter_plot(x_DEG_output = x_DEG_output,
                y_DEG_output = y_DEG_output,x_lab = "MEF2D-HNRNPUL1 cell lines protein FC", y_lab = "MEF2D-rearranged Gu et al. 2016 RNA_FC", title_lab = title_lab,filter_by = "pvalue",
                threshold = 0.01, points_color =  "#CDCC03", fitting_line_color = '#2C3E50', points_size = 2, genes_of_interest = genes_of_interest, labels_size = 4,
                labels_color = "#555555" , labeld_points_color =  "#35648A")
