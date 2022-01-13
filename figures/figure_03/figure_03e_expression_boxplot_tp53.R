library(Biobase)
library(tidyverse)
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
outputFolder <- file.path("figures", "output", "figure_03")

# create maijor output folder, if it doesn't exist
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder)

# load proteins
proteins <- readRDS(file = file.path("output", useVersion, "consensus_leukemic_clustering", "B", "proteins.B.RDS")) %>%
  .[, !pData(.)$Cell.Line.R %in% c("TANOUE", "TANOUE.BR2", "MN.60", "MN.60.BR2", "SEM.BR.NOPS")] %>%
  .[, grepl("^KMT2A", pData(.)$Subtype)]

# parameters
goi    <- "TP53"
groupA <- pData(proteins)$Subtype.Alt == "KMT2A.AFF1" 
groupALabel <- "KMT2A-AFF1"
groupAColor <- "#009acd"
groupAOpacity <- 1
groupB <- pData(proteins)$Subtype.Alt == "KMT2A.MLLT1" 
groupBLabel <- "KMT2A-MLLT1"
groupBColor <- "#009acd"
groupBOpacity <- 0.5
xLabel <- ""
yLabel <- expression(paste("rel. ", log[2], " expr."))
comparisons <- list(c("groupA", "groupB"))
comparisonMethod <- "t.test"


# set the groups
p_data <- proteins %>%
  pData() %>%
  mutate(group = if_else(groupA, "groupA",
                         if_else(groupB, "groupB", "none"))) %>%
  mutate(group = factor(group))

# prepare data
data <- proteins %>%
  exprs() %>%
  as.data.frame() %>%
  cbind(proteins %>% fData()) %>%
  filter(gene_symbol == goi) %>%
  gather(key = "proteomics_id", value = "value", 1:ncol(proteins)) %>%
  left_join(p_data, by = "proteomics_id")


# plot
TP53 <- data %>%
  ggplot(aes(x = group, y = value, color = group)) +
  # geom_violin(aes(alpha = group), color = "#444444", fill = "#009acd", linetype = 0) +
  scale_color_manual(values = c("groupA" = groupAColor, "groupB" = groupBColor)) +
  scale_alpha_manual(values = c("groupA" = groupAOpacity, "groupB" = groupBOpacity)) +
  scale_x_discrete(breaks = c("groupA", "groupB"), labels = c(groupALabel, groupBLabel)) +
  geom_jitter(size = 1.5, stroke = 0, alpha = 0.8, color = "gray10") +
  stat_compare_means(comparisons = comparisons,
                     label = "p.format",
                     method = comparisonMethod,
                     label.y = max(data$value, na.rm = TRUE) * 1.2,
                     bracket.size = 0.5) +
  xlab(xLabel) +
  ylab(yLabel) +
  ggtitle(goi) +
  theme(legend.position = "none")


#view plot
TP53

# export to svg
ggsave(
  "tp53_expression_boxplot.svg",
  plot = TP53,
  device = NULL,
  path = outputFolder,
  scale = 1,
  width = 4,
  height = 4,
  units = c("cm"),
  dpi = 600,
  limitsize = FALSE)
