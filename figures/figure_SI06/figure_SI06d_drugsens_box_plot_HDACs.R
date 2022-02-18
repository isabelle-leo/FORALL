library(Biobase)
library(tidyverse)
library(ggpubr)
library(ggthemes)
library(svglite)
library(viridis)
library(scales)
library(readxl)

source(file.path("functions", "getColorScheme.R"))
source(file.path("figures", "all_theme.R"))
theme_set(all_theme())

# show_col(viridis_pal(option = "B")(20))

# constants
useVersion   <- "publication"
colorScheme  <- getColorScheme()

#Choose output folder
outputFolder <- file.path("figures", "output", "figure_SI06")

# create major output folder, if it doesn't exist
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder)

# load drugsens
# drugsens <- readRDS(file = file.path("output", useVersion, "dsrt", "drugsens_2020-12-19.RDS"))

# parameters
drug <- "FIMM003787" #FIMM003783	Panobinostat FIMM100399	Quisinostat FIMM003787	Vorinostat
groupALabel <- "MEF2D-HNRNPUL1"
groupAColor <- "#b46cf0"
groupAOpacity <- 1
groupBLabel <- "other BCP-ALL and B-ALL"
groupBColor <- "#b46cf0"
groupBOpacity <- 0.5
xLabel <- ""
yLabel <- "sDSS"
comparisons <- list(c("groupA", "groupB"))
comparisonMethod <- "t.test"

# load from Excel table to include Kasumi-7
meta <- readRDS(file = file.path("output", useVersion, "proteins.RDS")) %>%
  pData() %>%
  dplyr::select(Cell.Line.R, Subtype_Paper, Type_Paper) %>%
  rbind(data.frame(Cell.Line.R = "KASUMI.7",
                   Subtype_Paper = "MEF2D-HNRNPUL1",
                   Type_Paper = "BCP-ALL"))

data <- read_excel(path = file.path("output", useVersion, "dsrt", "sDSS_2D_matrix_31_cells_with_annotation_RJ.xlsx"),
                   sheet = "sDSS_2D_matrix_31_cells_wit (2)") %>%
  .[, -(3:11)] %>%
  gather(key = "Cell.Line.R", value = "value", -FIMM.ID, -DRUG.NAME) %>%
  left_join(meta, by = "Cell.Line.R") %>%
  filter(Type_Paper != "T-ALL") %>%
  filter(FIMM.ID == drug) %>%
  mutate(group = if_else(Subtype_Paper == "MEF2D-HNRNPUL1", "groupA", "groupB"))

# plot
p <- data %>%
  ggplot(aes(x = group, y = value, alpha = group)) +
  geom_hline(yintercept = 0, color = "#444444") +
  geom_violin(aes(alpha = group), color = "#444444", fill = "#b46cf0", linetype = 0) +
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
  ggtitle(data$DRUG.NAME %>% unique() %>% paste(collapse = ", ")) +
  theme(legend.position = "none")

#view plots
p

# export to svg
ggsave(
  "bryostatin_1_drugsens_boxplot.svg",
  plot = p,
  device = NULL,
  path = outputFolder,
  scale = 1,
  width = 4,
  height = 4,
  units = c("cm"),
  dpi = 600,
  limitsize = FALSE)
