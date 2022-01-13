library(tidyverse)
library(ggrepel)
library(ggpubr)
library(ggthemes)
library(svglite)

source(file.path("functions", "getColorScheme.R"))
source(file.path("figures", "all_theme.R"))
theme_set(all_theme())

# constants
useVersion   <- "publication"
colorScheme  <- getColorScheme()

# Set outputfolder
outputFolder <- file.path("figures", "output", "figure_02")

# create maijor output folder, if it doesn't exist
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder)

# load data
data <- read_tsv(file = file.path("data", "mrna_protein", "ALL_corresponding_RNA_protein_correlation.txt")) %>%
  filter(!is.na(spearman_corr)) %>%
  arrange(desc(spearman_corr)) %>%
  mutate(rank = 1:nrow(.))

# parameters
labelColor <- "#555555"   
segmentColor <- "#AAAAAA"
labelTop <- 10
labelBottom <- 10

# labeling #Genes of Interest from the txt file
genesToLabel <- read.table(file = file.path("meta", "genes_of_interest_rna_protein_corr.txt"),
                           header = FALSE) %>%
  .[[1]]

# plot
p <- ggplot(mapping = aes(x = rank, y = spearman_corr)) +
  geom_point(data = data %>% filter(!gene %in% genesToLabel),
             size = 2,
             stroke = 0,
             alpha = 0.2,
             color = "#CDCC03") + #Yellow:"#CDCC03" Original #453781 
  # scale_color_manual(values = colorScheme$Class.explained) +
  geom_label_repel(data = data %>% filter(gene %in% genesToLabel),
                   aes(label = gene),
                   color = labelColor,
                   size = 4, #change gene label size
                   direction = "both",
                   box.padding = 0.05,
                   point.padding = 0.6,
                   label.size = 0.15,
                   segment.size = 0.6,
                   segment.color = segmentColor,
                   min.segment.length = 0.2,
                   force = 20) +
  geom_point(data = data %>% filter(gene %in% genesToLabel),
             size = 3.5,
             stroke = 0,
             alpha = 0.7,
             color = "#35648A") + #Blue:"#35648A" Original #29AF7F
  xlab("Rank") +
  ylab("Spearman mRNA-protein correlation")

#view the plot
p

# export to svg
ggsave(
  "rna_protein_corr_ranks.svg",
  plot = p,
  device = NULL,
  path = outputFolder,
  scale = 1,
  width = 12,
  height = 8,
  units = c("cm"),
  dpi = 600,
  limitsize = FALSE)
