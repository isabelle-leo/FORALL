# dependencies
library(tidyverse)
library(ggrepel)
library(ggthemes)
library(svglite)

source(file.path("functions", "getColorScheme.R"))
source(file.path("figures", "all_theme.R"))
theme_set(all_theme())

# constants
useVersion   <- "publication"
colorScheme  <- getColorScheme()

#Choose output folder
outputFolder <- file.path("figures", "output", "figure_SI05")

# create major output folder, if it doesn't exist
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder)

# aesthetic constants
lineColor <- "#999999"
labelColor <- "#555555"
segmentColor <- "#AAAAAA"
# xLabel <- expression(paste(log[2], " fold-change MEF2D-HNRNPUL1 vs. other BCP-ALL"))
# yLabel <- expression(paste(-log[10]~italic(p), "-value"))
xLabel <- expression(paste("Corr B - Corr T"))
yLabel <- expression(paste(-log[10]~italic(AdjP), "-value"))
densityColor <- "Red"

# labeling
# genesToLabel <- c("HDAC9", "NFATC1", "LCK", "MEF2D", "MEF2C", "IDH", "NRIP1", "AMER2", "PTPRZ1", "CKM", "CMTM3", 
#                  "ADSS1", "TCF4", "MAGI1", "PLCL2", "CAMK4", "TEAD1", "ARHGAP6", "PTK2", "FBXO27", "FRMPD1",
#                 "BCL11A", "CD5", "DNTT", "DGKH", "DSTYK")


# read DEqMS result
# file <- file.path("output", useVersion, "de", "MEF2D-HNRNPUL1_vs_BCP-ALL", "MEF2D.HNRNPUL1-BCPALL_DEqMS_table.txt")
# 
# data <- read_tsv(file = file) %>%
#   mutate(negLog10.adj.P.Val = -log10(adj.P.Val))

DGCA_output <- readRDS(file = file.path("output", useVersion, "dsrt", "protein_sDSS_B_vs_T_pearson_ddcorAll_output.RDS"))
#DGCA_output <- readRDS(file = file.path("output", useVersion, "dsrt", "protein_sDSS_B_vs_T_spearman_ddcorAll_output.RDS"))
DGCA_output_sign <- DGCA_output[DGCA_output$pValDiff_adj < 0.15, ] #filter here pValDiff_adj
DGCA_output2 <- DGCA_output
DGCA_output2 <- DGCA_output2[DGCA_output2$Gene1 %in% DGCA_output_sign$Gene1, ]
DGCA_output2 <- DGCA_output2[DGCA_output2$Gene2 %in% DGCA_output_sign$Gene2, ]
DGCA_output <- DGCA_output2
rm("DGCA_output_sign", "DGCA_output2")

# prepare plot df
# gene <- paste(DGCA_output$Gene1, gsub("+.+\\.\\.\\.", "",DGCA_output$Gene2), sep = "-")
# data <- data.frame("gene" = gene, "logFC" = DGCA_output$B_cor - DGCA_output$T_cor, "negLog10.adj.P.Val" = -log10(DGCA_output$pValDiff_adj), stringsAsFactors = F, row.names = gene)
# genesToLabel <- data$gene[data$negLog10.adj.P.Val > 1] #Filter here

# prepare plot df and filter based on pval and also correlation.
gene <- paste(DGCA_output$Gene1, gsub("+.+\\.\\.\\.", "",DGCA_output$Gene2), sep = "-")
data <- data.frame("gene" = gene, "logFC" = DGCA_output$B_cor - DGCA_output$T_cor, "negLog10.adj.P.Val" = -log10(DGCA_output$pValDiff_adj), stringsAsFactors = F, "B_cor" = DGCA_output$B_cor, "T_cor" = DGCA_output$T_cor, row.names = gene)
#genesToLabel <- data$gene[(data$negLog10.adj.P.Val > 1) & (abs(data$B_cor) > 0.3) & (abs(data$T_cor) > 0.3)]

# labeling from the txt file
genesToLabel <- read.table(file = file.path("meta", "DCGA.txt"),
                           header = FALSE) %>% .[[1]]
  

# calculate density
# adapted from https://slowkow.com/notes/ggplot2-color-by-density/
getDensity <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

data$density <- getDensity(x = data$logFC,
                           y = data$negLog10.adj.P.Val,
                           h = c(2, 2),
                           n = 200)

# plot
p <- ggplot(mapping = aes(x = logFC, y = negLog10.adj.P.Val, color = density)) +
  geom_hline(yintercept = -log10(0.1), color = lineColor, linetype = "dotted") + #change y axis
  geom_vline(xintercept = 0.58, color = lineColor, linetype = "dotted") +
  geom_vline(xintercept = -0.58, color = lineColor, linetype = "dotted") +
  geom_point(data = data %>% filter(!gene %in% genesToLabel),
             size = 1.5,
             stroke = 0,
             alpha = 0.9) +
  scale_color_gradient_tableau(densityColor) +
  geom_label_repel(data = data %>% filter(gene %in% genesToLabel),
                   aes(label = gene),
                   color = labelColor,
                   size = 2,
                   direction = "both",
                   box.padding = 0.05,
                   point.padding = 0.2,
                   label.size = 0.15,
                   segment.size = 0.6,
                   segment.color = segmentColor,
                   max.overlaps = Inf,
                   min.segment.length = 0.2) +
  geom_point(data = data %>% filter(gene %in% genesToLabel),
             size = 3,
             stroke = 0,
             alpha = 0.7,
             color = colorScheme$blue) +
  xlab(xLabel) +
  ylab(yLabel) +
  theme(legend.position = "none")

#view plot
p

# export to svg
ggsave(
  "volcano_DCGA.svg",
  plot = p,
  device = NULL,
  path = outputFolder,
  scale = 1,
  width = 8,
  height = 8,
  units = c("cm"),
  dpi = 600,
  limitsize = FALSE)
