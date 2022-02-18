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
outputFolder <- file.path("figures", "output", "figure_06")

# create major output folder, if it doesn't exist
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder)

# aesthetic constants
lineColor <- "#999999"
labelColor <- "#555555"
segmentColor <- "#AAAAAA"
xLabel <- expression(paste(log[2], " fold-change MEF2D-HNRNPUL1 vs. other BCP-ALL"))
yLabel <- expression(paste(-log[10]~italic(p), "-value"))
densityColor <- "Red"

# labeling
genesToLabel <- c("HDAC9", "NFATC1", "LCK", "MEF2D", "MEF2C", "IDH", "NRIP1", "AMER2", "PTPRZ1", "CKM", "CMTM3", 
                  "ADSS1", "TCF4", "MAGI1", "PLCL2", "CAMK4", "TEAD1", "ARHGAP6", "PTK2", "FBXO27", "FRMPD1",
                  "BCL11A", "CD5", "DNTT", "DGKH", "DSTYK")

# read DEqMS result
file <- file.path("output", useVersion, "de", "MEF2D_HNRNPUL1_vs_BCP-ALL", "MEF2D.HNRNPUL1-BCPALL_DEqMS_table.txt")

data <- read_tsv(file = file) %>%
  mutate(negLog10.adj.P.Val = -log10(adj.P.Val))


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
  geom_hline(yintercept = -log10(0.05), color = lineColor, linetype = "dotted") +
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
  "volcano_mef2d_bcp_all.svg",
  plot = p,
  device = NULL,
  path = outputFolder,
  scale = 1,
  width = 8,
  height = 8,
  units = c("cm"),
  dpi = 600,
  limitsize = FALSE)
