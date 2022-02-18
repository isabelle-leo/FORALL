# dependencies
library(Biobase)
library(tidyverse)
library(ggpubr)
library(ggthemes)
library(svglite)

source(file.path("figures", "all_theme.R"))
theme_set(all_theme())

source(file.path("functions", "getColorScheme.R"))

# constants
useVersion   <- "publication"
colorScheme  <- getColorScheme()

#Choose output folder
outputFolder <- file.path("figures", "output", "figure_SI05")

# create major output folder, if it doesn't exist
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder)

# load proteins and drugs
proteins <- readRDS(file = file.path("output", useVersion, "consensus_leukemic_clustering", "all", "proteins.RDS"))
# filter example to remove B-ALL
proteins <- proteins[, !pData(proteins)$Cell.Line.R %in% c("SEM.BR.NOPS")]
# load drug sensitivity for drugs
drugsens <- readRDS(file = file.path("data", "dsrt", "drugsens_2021-03-22.RDS")) 

#pick protein and drugs to plot
gois <- c("PSMC2") #, "NRIP1"
drug <- c("FIMM133887", "FIMM136375", "FIMM136458") #BAY_87-2243 FIMM136458 #Mubritinib FIMM136375 #NMS FIMM136458

xLabel <- expression(paste("relative ", log[2], " expression"))
yLabel <- "sDSS"

# prepare and filter data
proteins_df <- proteins %>%
  exprs() %>%
  as.data.frame() %>%
  set_names(pData(proteins)$Cell.Line.R) %>%
  cbind(proteins %>% fData()) %>%
  filter(gene_symbol %in% gois) %>%
  gather(key = "Cell.Line.R", value = "expression_value", 1:ncol(proteins)) %>%
  left_join(proteins %>% pData(), by = "Cell.Line.R") %>%
  mutate(is_cell_line_replicate = grepl("\\.BR[0-9]$", Cell.Line.R)) %>%
  mutate(cell_line = gsub("\\.BR[0-9]$", "", Cell.Line.R))

drugsens_df <- drugsens %>%
  exprs() %>%
  as.data.frame() %>%
  cbind(drugsens %>% fData()) %>%
  filter(FIMM.ID %in% drug) %>%
  gather(key = "Cell.Line.R", value = "sdss_value", 1:ncol(drugsens))

joined_df <- proteins_df %>%
  full_join(drugsens_df, by = "Cell.Line.R") %>%
  filter(!is.na(expression_value)) %>%
  filter(!is.na(sdss_value))

# plot
p <- joined_df %>%
  ggplot(aes(x = expression_value,
             y = sdss_value,
             color = Subtype_Paper)) +
  geom_smooth(method = "lm", formula = "y ~ x", color = "#555555", size = 0.6, fill = "#CCCCCC") +
  geom_point(size = 2, stroke = 0, alpha = 0.8) +
  stat_cor(label.sep = ",",
           method = "pearson",
           geom = "text",
           color = "#555555",
           size = 2,
           label.y = max(joined_df$sdss_value) * 1.2) +
  xlab(xLabel) +
  ylab(yLabel) +
  ggtitle(gois %>% unique() %>% paste(collapse = ", ")) +
  facet_wrap(~ DRUG.NAME) +
  scale_color_manual(name = "Subtype", values = colorScheme$subtype.alt) +
  theme(legend.position = "bottom")

# view plot
p

# export to svg
ggsave(
  "figure_SI05_drugsens_corr_Mubritinib_PSMC_pvals_rev.svg",
  plot = p,
  device = NULL,
  path = outputFolder,
  scale = 1,
  width = 10,
  height = 8,
  units = c("cm"),
  dpi = 600,
  limitsize = FALSE)
