library(Biobase)
library(tidyverse)
library(ggpubr)
library(ggthemes)

source(file.path("figures", "all_theme.R"))
theme_set(all_theme())

source(file.path("functions", "getColorScheme.R"))

# constants
useVersion   <- "publication"
colorScheme  <- getColorScheme()

#Choose output folder
outputFolder <- file.path("figures", "output", "figure_05")

# create major output folder, if it doesn't exist
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder)

# load proteins and drugs
proteins <- readRDS(file = file.path("output", useVersion, "consensus_leukemic_clustering", "B", "proteins.B.RDS"))
#proteins <- readRDS(file = file.path("output", useVersion, "consensus_leukemic_clustering", "all", "proteins.RDS"))
# filter example to remove B-ALL
proteins <- proteins[, !pData(proteins)$Cell.Line.R %in% c("TANOUE", "TANOUE.BR2", "MN.60", "MN.60.BR2", "SEM.BR.NOPS")]
# load drug sensitivity for drugs just BCP-ALL
drugsens <- readRDS(file = file.path("data", "dsrt", "drugsens_2021-03-22.RDS"))
#drugsens <- readRDS(file = file.path("output", useVersion, "dsrt", "drugsens_2020-12-18.RDS")) #alternatively, adjust the file to look at just BCP-ALL correlations or other data subset

#pick protein and drugs to plot
gois <- c("TCP1", "HSP90AB1")
drug <- c("FIMM133877")

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
  geom_smooth(method = "lm", formula = "y ~ x", color = "#555555") +
  geom_point(aes(shape = DRUG.NAME), size = 3, stroke = 0, alpha = 0.8) +
  stat_cor(label.sep = ",",
           method = "pearson",
           geom = "text",
           color = "#555555",
           size = rel(3),
           label.y = max(joined_df$sdss_value) * 1.2) +
  xlab(xLabel) +
  ylab(yLabel) +
  ggtitle(joined_df$DRUG.NAME %>% unique() %>% paste(collapse = ", ")) +
  facet_wrap(~ gene_symbol) +
  scale_color_manual(name = "Subtype", values = colorScheme$subtype.alt) +
  theme(legend.position = "bottom")

# view plot
p

# export to svg
library(svglite)
ggsave(
  "figure_05_drugSensExpressionScatterPlot_OSU-03012_TCP1_HSP90AB1.svg",
  plot = p,
  device = NULL,
  path = outputFolder,
  scale = 1,
  width = 12,
  height = 10,
  units = c("cm"),
  dpi = 600,
  limitsize = FALSE)

# export to png
ggsave(
  "figure_05_drugSensExpressionScatterPlot_OSU-03012_TCP1_HSP90AB1.png",
  plot = p,
  device = NULL,
  path = outputFolder,
  scale = 1,
  width = 8,
  height = 8,
  units = c("in"),
  dpi = 600,
  limitsize = FALSE)

