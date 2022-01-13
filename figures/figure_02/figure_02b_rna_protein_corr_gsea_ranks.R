####################################################################################
####################################################################################
# rnaProteinCorrRanks GSEA KEGG 
####################################################################################
####################################################################################

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

# create major output folder, if it doesn't exist
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder)

# load data
data <- read_tsv(file = file.path("data", "mrna_protein", "fgea_ALL_corresponding_RNA_protein_correlation_8981_genes.txt")) %>%
  filter(!is.na(ES)) %>%
  arrange(desc(ES)) %>%
  mutate(rank = 1:nrow(.))

# parameters
labelColor <- "#555555"
segmentColor <- "#AAAAAA"


# labeling #Genes of Interest from the txt file
# genesToLabel <- read_tsv(file = file.path("data", "mrna_protein", "Genes of Interest ALL_2020-10-12.txt"),
#                          col_names = F) %>% .[[1]]

# labeling specific pathways
genesToLabel <- c("B_CELL_RECEPTOR_SIGNALING_PATHWAY", 
                  "CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
                  "T_CELL_RECEPTOR_SIGNALING_PATHWAY",
                  "P53_SIGNALING_PATHWAY",
                  "CALCIUM_SIGNALING_PATHWAY",
                  "LEUKOCYTE_TRANSENDOTHELIAL_MIGRATION",
                  "FC_EPSILON_RI_SIGNALING_PATHWAY",
                  "PHOSPHATIDYLINOSITOL_SIGNALING_SYSTEM",
                  "MTOR_SIGNALING_PATHWAY",
                  "OXIDATIVE_PHOSPHORYLATION",
                  "RIBOSOME",
                  "PROTEASOME",
                  "SPLICEOSOME",
                  "BASAL_TRANSCRIPTION_FACTORS",
                  "NOTCH_SIGNALING_PATHWAY") 

#Label top and bottom pathways
# labelTop <- 5
# labelBottom <- 5
# genesToLabel <- NULL
# if (is.null(genesToLabel)) {
#   genesToLabel <- c()
#   
#   if (is.numeric(labelTop)) {
#     genesToLabel <- c(genesToLabel, data$pathway[seq_len(labelTop)])
#   }
#   
#   if (is.numeric(labelBottom)) {
#     genesToLabel <- c(genesToLabel, rev(data$pathway)[seq_len(labelBottom)])
#   }
# }

# plot
p <- ggplot(mapping = aes(x = rank, y = ES)) +
  geom_point(data = data %>% filter(!pathway %in% genesToLabel),
             size = 2,
             stroke = 0,
             alpha = 0.7,
             color = "#CDCC03") + #Yellow:"#CDCC03" Original #453781 #sandybrown 
  # scale_color_manual(values = colorScheme$Class.explained) +
  geom_label_repel(data = data %>% filter(pathway %in% genesToLabel),
                   aes(label = pathway),
                   color = labelColor,
                   size = 2,
                   direction = "both",
                   box.padding = 0.05,
                   point.padding = 0.6,
                   label.size = 0.15,
                   segment.size = 0.6,
                   segment.color = segmentColor,
                   min.segment.length = 0.2,
                   force = 5) +
  geom_point(data = data %>% filter(pathway %in% genesToLabel),
             size = 2.5,
             stroke = 0,
             alpha = 0.9,
             color = "#35648A") + #Blue:"#35648A" Original #29AF7F
  xlab("Rank") +
  ylab("Enrichment Score")

#view the plot
p

# export to svg
ggsave(
  "rna_protein_corr_GSEA_ranks.svg",
  plot = p,
  device = NULL,
  path = outputFolder,
  scale = 1,
  width = 6,
  height = 4,
  units = c("cm"),
  dpi = 600,
  limitsize = FALSE)

