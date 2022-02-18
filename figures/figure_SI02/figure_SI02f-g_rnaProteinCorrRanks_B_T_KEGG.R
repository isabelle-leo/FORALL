library(tidyverse)
library(ggrepel)
library(ggpubr)
library(ggthemes)

source(file.path("functions", "getColorScheme.R"))
source(file.path("figures", "all_theme.R"))
theme_set(all_theme())

# constants
useVersion   <- "publication"
colorScheme  <- getColorScheme()

# Set outputfolder
outputFolder <- file.path("figures", "output", "figure_SI02")
# create major output folder, if it doesn't exist
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder)

# load data
data <- read_tsv(file = file.path("data", "mrna_protein", "ALL_corresponding_RNA_protein_correlation_2021_03_02_19_29_30.txt")) %>%
  filter(!is.na(spearman_corr)) %>%
  arrange(desc(spearman_corr)) %>%
  mutate(rank = 1:nrow(.))

# parameters
labelColor <- "#555555"   
segmentColor <- "#AAAAAA"
labelTop <- 10
labelBottom <- 10

# labeling #Genes of Interest from the txt file
genesToLabel <- read_tsv(file = file.path("data", "mrna_protein", "KEGG_B_CELL_RECEPTOR.txt"),
                         col_names = F) %>% .[[1]]

# plot
pB <- ggplot(mapping = aes(x = rank, y = spearman_corr)) +
  geom_point(data = data %>% filter(!gene %in% genesToLabel),
             size = 3,
             stroke = 0,
             alpha = 0.2,
             color = "#CDCC03") + #Yellow:"#CDCC03" Original #453781 
  # scale_color_manual(values = colorScheme$Class.explained) +
  geom_label_repel(data = data %>% filter(gene %in% genesToLabel),
                   aes(label = gene),
                   color = labelColor,
                   size = rel(3), #change gene label size
                   direction = "both",
                   box.padding = 0.05,
                   point.padding = 0.6,
                   label.size = 0.15,
                   segment.size = 0.6,
                   segment.color = segmentColor,
                   min.segment.length = 0.2,
                   force = 5) +
  geom_point(data = data %>% filter(gene %in% genesToLabel),
             size = 5.5,
             stroke = 0,
             alpha = 0.7,
             color = "#35648A") + #Blue:"#35648A" Original #29AF7F
  xlab("Rank") +
  ylab("Spearman mRNA-protein correlation")

#view plot
pB

# export to png
ggsave(
  "supp_figure_2C_Prot_mRNA_Correlation_vs_rank_KEGG_B_CELL_RECEPTOR_batch_corr.png",
  plot = pB,
  device = NULL,
  path = outputFolder,
  scale = 1,
  width = 6,
  height = 6,
  units = c("in"),
  dpi = 600,
  limitsize = FALSE)

# export to svg
library(svglite)
ggsave(
  "supp_figure_2C_Prot_mRNA_Correlation_vs_rank_KEGG_B_CELL_RECEPTOR_batch_corr.svg",
  plot = pB,
  device = NULL,
  path = outputFolder,
  scale = 1,
  width = 6,
  height = 6,
  units = c("in"),
  dpi = 600,
  limitsize = FALSE)

####################################################################################
####################################################################################
# Supp_Figure_2C_Prot_mRNA_Correlation_vs_rank_KEGG_T_CELL_RECEPTOR
####################################################################################
####################################################################################
source(file.path("functions", "getColorScheme.R"))
source(file.path("figures", "all_theme.R"))
theme_set(all_theme())

# constants
name         <- "[Produce Supp Figure 2C] "
useVersion   <- "publication"
colorScheme  <- getColorScheme()

# Set outputfolder
outputFolder <- file.path("figures", "output")
# create major output folder, if it doesn't exist
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder)

# load data
data <- read_tsv(file = file.path("data", "mrna_protein", "ALL_corresponding_RNA_protein_correlation_2021_03_02_19_29_30.txt")) %>%
  filter(!is.na(spearman_corr)) %>%
  arrange(desc(spearman_corr)) %>%
  mutate(rank = 1:nrow(.))

# parameters
labelColor <- "#555555"   
segmentColor <- "#AAAAAA"
labelTop <- 10
labelBottom <- 10

# labeling #Genes of Interest from the txt file
genesToLabel <- read_tsv(file = file.path("data", "mrna_protein", "KEGG_T_CELL_RECEPTOR.txt"),
                         col_names = F) %>% .[[1]]

# plot
pT <- ggplot(mapping = aes(x = rank, y = spearman_corr)) +
  geom_point(data = data %>% filter(!gene %in% genesToLabel),
             size = 3,
             stroke = 0,
             alpha = 0.2,
             color = "#CDCC03") + #Yellow:"#CDCC03" Original #453781 
  # scale_color_manual(values = colorScheme$Class.explained) +
  geom_label_repel(data = data %>% filter(gene %in% genesToLabel),
                   aes(label = gene),
                   color = labelColor,
                   size = rel(3), #change gene label size
                   direction = "both",
                   box.padding = 0.05,
                   point.padding = 0.6,
                   label.size = 0.15,
                   segment.size = 0.6,
                   segment.color = segmentColor,
                   min.segment.length = 0.2,
                   force = 5) +
  geom_point(data = data %>% filter(gene %in% genesToLabel),
             size = 5.5,
             stroke = 0,
             alpha = 0.7,
             color = "#35648A") + #Blue:"#35648A" Original #29AF7F
  xlab("Rank") +
  ylab("Spearman mRNA-protein correlation")

#view plot
pT

# export to png
ggsave(
  "supp_figure_2C_Prot_mRNA_Correlation_vs_rank_KEGG_T_CELL_RECEPTOR_Batch_Corr.png",
  plot = pT,
  device = NULL,
  path = outputFolder,
  scale = 1,
  width = 6,
  height = 6,
  units = c("in"),
  dpi = 600,
  limitsize = FALSE)

# export to svg
library(svglite)
ggsave(
  "supp_figure_2C_Prot_mRNA_Correlation_vs_rank_KEGG_T_CELL_RECEPTOR_Batch_Corr.svg",
  plot = pT,
  device = NULL,
  path = outputFolder,
  scale = 1,
  width = 6,
  height = 6,
  units = c("in"),
  dpi = 600,
  limitsize = FALSE)
