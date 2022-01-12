library(Biobase)
library(tidyverse)
library(ggpubr)
library(ggthemes)
library(svglite)

sourceDir <- function(path, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    source(file.path(path, nm), ...)
  }
}
sourceDir(path = "functions")

source(file.path("figures", "all_theme.R"))
theme_set(all_theme())

# constants
useVersion   <- "publication"
colorScheme  <- getColorScheme()

# set output folder
outputFolder <- file.path("figures", "output", "figure_01")

# create major output folder, if it doesn't exist
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder)

# load proteins (highly variable)
proteins <- readRDS(file = file.path("output", useVersion, "consensus_leukemic_clustering", "all", "proteins.RDS")) %>%
  .[fData(.)$hasHighVariance, ] %>%
  removeNAsFromESet()

# pca
pca <- prcomp(x = t(exprs(proteins)),
              scale. = TRUE)

variances <- pca %>%
  summary() %>%
  .$importance %>%
  .[1, ] %>%
  .^2

variances_explained <- round(100 * variances / sum(variances), 2)

df <- pca$x %>%
  as.data.frame() %>%
  rownames_to_column("proteomics_id") %>%
  left_join(pData(proteins), by = "proteomics_id") %>%
  dplyr::select(PC1, PC2, Cell_Line_Name_Paper, Type_Paper)

# plot
p <- df %>%
  ggplot(aes(x = PC1, y = PC2, color = Type_Paper)) +
  geom_vline(xintercept = 0, linetype = 2, color = "#AAAAAA") +
  geom_hline(yintercept = 0, linetype = 2, color = "#AAAAAA") +
  geom_point(size = 2.5, stroke = 0, alpha = 0.8) +
  scale_color_manual(values = colorScheme$type) +
  xlab(paste0("PC1 (", variances_explained[["PC1"]], "%)")) +
  ylab(paste0("PC2 (", variances_explained[["PC2"]], "%)"))
p

# export
ggsave(
  "pca.svg",
  plot = p,
  device = NULL,
  path = outputFolder,
  scale = 1,
  width = 9,
  height = 6,
  units = c("cm"),
  dpi = 600,
  limitsize = FALSE)
