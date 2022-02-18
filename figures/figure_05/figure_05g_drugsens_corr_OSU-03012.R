library(Biobase)
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

#Choose output folder
outputFolder <- file.path("figures", "output", "figure_05")

# create major output folder, if it doesn't exist
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder)


# load drug-protein correlation data (calculated using BCP-ALL)
drug_protein_corr <- readRDS(file = file.path("output", useVersion, "dsrt", "drugsens_protein_corr_ENS99_Filtered_2020-12-19.RDS"))

# parameters
drugId <- "FIMM133877" #FIMM133877...OSU-03012
labelColor <- "#555555"
segmentColor <- "#AAAAAA"
labelTop <- 10
labelBottom <- 10

# labeling top 10 or bottom 10
genesToLabel <- NULL

# labeling #FIMM133877...OSU-03012
genesToLabel <- c("TCP1", "CCT2", "CCT3", "CCT4", "CCT6A", "HSP90AB1", "PDPK1","AKT1","AKT2", "AKT3","PDPK1","PIK3CA","PIK3CG","PKN2","RPS6KA3","RPS6KB1","SGK1","TSC2") #"MAPK1" "MAPK3"

#genesToLabel <- read.table(".txt",header=F, check.names = FALSE, stringsAsFactors = F,na.strings=c("", "NA", "Inf"), sep="\t")[,1]

# prepare data
drug_protein_corr_tidy <- drug_protein_corr %>%
  exprs() %>%
  as.data.frame() %>%
  cbind(drug_protein_corr %>% fData()) %>%
  gather(key = "drug", value = "r", 1:ncol(drug_protein_corr)) %>%
  left_join(drug_protein_corr %>% pData(), by = c("drug" = "FIMM.ID")) %>%
  filter(!is.na(gene_symbol)) %>%
  filter(drug == drugId) %>%
  arrange(desc(r)) %>%
  mutate(rank = seq_len(nrow(.))) %>%
  mutate(label = if_else(gene_symbol %in% genesToLabel, gene_symbol, ""))

if (is.null(genesToLabel)) {
  genesToLabel <- c()
  
  if (is.numeric(labelTop)) {
    genesToLabel <- c(genesToLabel, drug_protein_corr_tidy$gene_symbol[seq_len(labelTop)])
  }
  
  if (is.numeric(labelBottom)) {
    genesToLabel <- c(genesToLabel, rev(drug_protein_corr_tidy$gene_symbol)[seq_len(labelBottom)])
  }
}


# plot
p <- ggplot(mapping = aes(x = rank, y = r)) +
  geom_point(data = drug_protein_corr_tidy %>% filter(!gene_symbol %in% genesToLabel),
             size = 2,
             stroke = 0,
             alpha = 0.7,
             color = "#A7C6DD") +
  geom_label_repel(data = drug_protein_corr_tidy %>% filter(gene_symbol %in% genesToLabel),
                   aes(label = gene_symbol),
                   color = labelColor,
                   size = 2,
                   direction = "both",
                   box.padding = 0.05,
                   point.padding = 0.6,
                   label.size = 0.15,
                   segment.size = 0.6,
                   segment.color = segmentColor,
                   min.segment.length = 0.2) +
  geom_point(data = drug_protein_corr_tidy %>% filter(gene_symbol %in% genesToLabel),
             size = 2.5,
             stroke = 0,
             alpha = 0.9,
             color = "#b46cf0") +
  xlab("Rank") +
  ylab("Correlation") +
  ggtitle(drug_protein_corr_tidy$DRUG.NAME %>% unique() %>% paste(collapse = ", "))


#view plots
p

# export to svg
ggsave(
  "figure_05_drugsens_corr_OSU-03012_rev.svg",
  plot = p,
  device = NULL,
  path = outputFolder,
  scale = 1,
  width = 6,
  height = 4,
  units = c("cm"),
  dpi = 600,
  limitsize = FALSE)
