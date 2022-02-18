# dependencies
library(Biobase)
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
outputFolder <- file.path("figures", "output", "figure_SI04")
# create major output folder, if it doesn't exist
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder)

# load drug-protein correlation data
#drug_protein_corr <- readRDS(file = file.path("output", useVersion, "dsrt", "drugsens_protein_corr.RDS")) #ENS92 Search Corr.
drug_protein_corr <- readRDS(file = file.path("output", useVersion, "dsrt", "drugsens_protein_corr_ENS99_Filtered_2021-03-22.RDS")) #Latest ENS99 Search Corr.

####################################################################################
# proteinDrugCorrelationByGoi_ABL1
####################################################################################

# parameters
goi <- "ABL1"
labelColor <- "#555555"
segmentColor <- "#AAAAAA"
labelTop <- 5
labelBottom <- 5

# labeling
drugsToLabel <- NULL #Plots top 10 or bottom 10

# prepare data
drug_protein_corr_tidy <- drug_protein_corr %>%
  exprs() %>%
  as.data.frame() %>%
  cbind(drug_protein_corr %>% fData()) %>%
  filter(gene_symbol == goi) %>%
  gather(key = "drug", value = "r", 1:ncol(drug_protein_corr)) %>%
  left_join(drug_protein_corr %>% pData(), by = c("drug" = "FIMM.ID")) %>%
  filter(!is.na(gene_symbol)) %>%
  arrange(desc(r)) %>%
  mutate(rank = seq_len(nrow(.)))

if (is.null(drugsToLabel)) {
  drugsToLabel <- c()
  
  if (is.numeric(labelTop)) {
    drugsToLabel <- c(drugsToLabel, drug_protein_corr_tidy$drug[seq_len(labelTop)])
  }
  
  if (is.numeric(labelBottom)) {
    drugsToLabel <- c(drugsToLabel, rev(drug_protein_corr_tidy$drug)[seq_len(labelBottom)])
  }
}


# plot
ABL1 <- ggplot(mapping = aes(x = rank, y = r)) +
  geom_point(data = drug_protein_corr_tidy %>% filter(!drug %in% drugsToLabel),
             size = 2,
             stroke = 0,
             alpha = 0.7,
             color = colorScheme$blue) +
  #scale_color_manual(values = colorScheme$Class.explained) +
  geom_label_repel(data = drug_protein_corr_tidy %>% filter(drug %in% drugsToLabel),
                   aes(label = DRUG.NAME),
                   color = labelColor,
                   size = rel(2),
                   direction = "both",
                   box.padding = 0.05,
                   point.padding = 0.6,
                   label.size = 0.15,
                   segment.size = 0.6,
                   segment.color = segmentColor,
                   min.segment.length = 0.2) +
  geom_point(data = drug_protein_corr_tidy %>% filter(drug %in% drugsToLabel),
             size = 2.5,
             stroke = 0,
             alpha = 0.9,
             color = colorScheme$red) +
  xlab("Rank") +
  ylab("Pearson correlation") +
  ggtitle(goi)

#view the plot
ABL1

# export to svg
library(svglite)
ggsave(
  "figure_SI05_proteinDrugCorrelationByGoi_ABL1_rev.svg",
  plot = ABL1,
  device = NULL,
  path = outputFolder,
  scale = 1,
  width = 6,
  height = 4,
  units = c("cm"),
  dpi = 600,
  limitsize = FALSE)

####################################################################################
# proteinDrugCorrelationByGoi_FLT3
####################################################################################

# parameters
goi <- "FLT3"
labelColor <- "#555555"
segmentColor <- "#AAAAAA"
labelTop <- 5
labelBottom <- 5

# labeling
drugsToLabel <- NULL #Plots top 10 or bottom 10

# prepare data
drug_protein_corr_tidy <- drug_protein_corr %>%
  exprs() %>%
  as.data.frame() %>%
  cbind(drug_protein_corr %>% fData()) %>%
  filter(gene_symbol == goi) %>%
  gather(key = "drug", value = "r", 1:ncol(drug_protein_corr)) %>%
  left_join(drug_protein_corr %>% pData(), by = c("drug" = "FIMM.ID")) %>%
  filter(!is.na(gene_symbol)) %>%
  arrange(desc(r)) %>%
  mutate(rank = seq_len(nrow(.)))

if (is.null(drugsToLabel)) {
  drugsToLabel <- c()
  
  if (is.numeric(labelTop)) {
    drugsToLabel <- c(drugsToLabel, drug_protein_corr_tidy$drug[seq_len(labelTop)])
  }
  
  if (is.numeric(labelBottom)) {
    drugsToLabel <- c(drugsToLabel, rev(drug_protein_corr_tidy$drug)[seq_len(labelBottom)])
  }
}


# plot
FLT3 <- ggplot(mapping = aes(x = rank, y = r)) +
  geom_point(data = drug_protein_corr_tidy %>% filter(!drug %in% drugsToLabel),
             size = 2,
             stroke = 0,
             alpha = 0.7,
             color = colorScheme$blue) +
  #scale_color_manual(values = colorScheme$Class.explained) +
  geom_label_repel(data = drug_protein_corr_tidy %>% filter(drug %in% drugsToLabel),
                   aes(label = DRUG.NAME),
                   color = labelColor,
                   size = rel(2),
                   direction = "both",
                   box.padding = 0.05,
                   point.padding = 0.6,
                   label.size = 0.15,
                   segment.size = 0.6,
                   segment.color = segmentColor,
                   min.segment.length = 0.2) +
  geom_point(data = drug_protein_corr_tidy %>% filter(drug %in% drugsToLabel),
             size = 2.5,
             stroke = 0,
             alpha = 0.9,
             color = colorScheme$red) +
  xlab("Rank") +
  ylab("Pearson correlation") +
  ggtitle(goi)

#view the plot
FLT3

# export to svg
library(svglite)
ggsave(
  "figure_SI05_proteinDrugCorrelationByGoi_FLT3_rev.svg",
  plot = FLT3,
  device = NULL,
  path = outputFolder,
  scale = 1,
  width = 6,
  height = 4,
  units = c("cm"),
  dpi = 600,
  limitsize = FALSE)
