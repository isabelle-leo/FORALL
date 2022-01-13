library(Biobase)
library(tidyverse)
library(ggrepel)
library(ggpmisc)
library(ggthemes)
library(svglite)

source(file.path("functions", "getColorScheme.R"))
source(file.path("figures", "all_theme.R"))
theme_set(all_theme())

# constants
useVersion   <- "publication"
colorScheme  <- getColorScheme()

#Choose output folder
outputFolder <- file.path("figures", "output", "figure_04")

# create major output folder, if it doesn't exist
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder)

# parameters
labelColor <- "#555555"
segmentColor <- "#333333"

# labeling including T-ALL
drugsToLabel <- read_tsv(file = file.path("meta", "drugs_to_label_fig4_variance.txt"), col_names = F) %>% .[[1]]

# drugsToLabel <- c("FIMM000227",
#                   "FIMM136384",
#                   "FIMM136514",
#                   "FIMM115484",
#                   "FIMM000484",
#                   "FIMM023797",
#                   "FIMM000406",
#                   "FIMM003775",
#                   "FIMM000519",
#                   "FIMM003747",
#                   "FIMM136458",
#                   "FIMM003730",
#                   "FIMM136513",
#                   "FIMM136521",
#                   "FIMM000321",
#                   "FIMM003761",
#                   "FIMM136412",
#                   "FIMM003771",
#                   "FIMM133851",
#                   "FIMM133860",
#                   "FIMM003742",
#                   "FIMM023798",
#                   "FIMM136469",
#                   "FIMM003763",
#                   "FIMM133795",
#                   "FIMM023808",
#                   "FIMM000264",
#                   "FIMM136388",
#                   "FIMM003707",
#                   "FIMM136466",
#                   "FIMM136567",
#                   "FIMM136615")

#drugsToLabel <- NULL #use to label none

# load data
drugsens <- readRDS(file = file.path("output", useVersion, "dsrt", "drugsens_2021-03-22.RDS")) %>% 
  .[, !pData(.)$Cell.Line.R %in% c("SEM.NOPS")]

modify_class_names <- function (names) {
  names %>%
    gsub("^[A-Z]\\.", "", .) %>%
    gsub("\\.", " ", .)
}

# mutate pheno data
fData(drugsens) <- drugsens %>%
  fData() %>%
  mutate(class_explained_paper = modify_class_names(Class.explained))

# calculate median and variance
rowVars <- function (x) {
  rowSums((x - rowMeans(x))^2) / (dim(x)[2] - 1)
}

medians_vars <- drugsens %>%
  exprs() %>%
  data.frame(median = rowMedians(.),
             variance = rowVars(.),
             .,
             fData(drugsens))


# plot
p <- medians_vars %>%
  ggplot(aes(x = variance, y = median, color = class_explained_paper, label = DRUG.NAME)) +
  geom_hline(yintercept = 0, color = "#777777", linetype = "dotted") +
  geom_point(data = medians_vars %>% filter(!FIMM.ID %in% drugsToLabel), color = "#AAAAAA", size = 1.5, stroke = 0, alpha = 0.8) +
  geom_point(data = medians_vars %>% filter(FIMM.ID %in% drugsToLabel), shape = 21, size = 2, color = "black", stroke = 0.8, fill = "white", alpha = 0.8) +
  geom_point(data = medians_vars %>% filter(FIMM.ID %in% drugsToLabel), aes(fill = class_explained_paper), shape = 21, size = 1.2, stroke = 0, alpha = 1) +
  scale_color_manual(values = colorScheme$Class.explained %>% set_names(modify_class_names(names(.)))) +
  scale_fill_manual(values = colorScheme$Class.explained %>% set_names(modify_class_names(names(.)))) +
  geom_label_repel(data = medians_vars %>% filter(FIMM.ID %in% drugsToLabel),
                   color = labelColor,
                   size = rel(4),
                   direction = "both",
                   box.padding = 0.05,
                   point.padding = 0.4,
                   label.size = 0.15,
                   segment.size = 0.6,
                   segment.color = segmentColor,
                   min.segment.length = 0.2,
                   max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                   force = 300) +
  xlab("sDSS variance across cell lines") +
  ylab("sDSS medians across cell lines") +
  theme(legend.position = "right", legend.title = element_blank())


# view plot
p

# export to svg
ggsave(
  "figure_04_drugsens_variance_rev.svg",
  plot = p,
  device = NULL,
  path = outputFolder,
  scale = 1,
  width = 14,
  height = 6,
  units = c("cm"),
  dpi = 600,
  limitsize = FALSE)
