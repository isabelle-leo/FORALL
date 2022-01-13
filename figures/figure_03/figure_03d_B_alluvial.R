library(Biobase)
library(d3Network)
library(networkD3)
library(tidyverse)
library(readxl)
library(RColorBrewer)
library(ggalluvial)
library(ggthemes)
library(svglite)

source(file.path("figures", "all_theme.R"))
theme_set(all_theme())

#get functions
sourceDir <- function(path, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    source(file.path(path, nm), ...)
  }
}
sourceDir(path = "functions")

color_scheme <- getColorScheme()
useVersion <- "publication"

# set output folder
outputFolder <- file.path("figures", "output", "figure_03")

# create major output folder, if it doesn't exist
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder)


#================
proteinsMeta <- readRDS(file = file.path("output", useVersion, "consensus_leukemic_clustering", "all", "proteins.RDS")) %>%
  pData() %>%
  filter(!isReplicate) %>%
  filter(Type_Paper != "EBV") %>%
  dplyr::select(Cell.Line.R, Subtype_Paper, Type_Paper)

proteinClusters <- read_tsv(file = file.path("meta", "cell_line_clustering.txt"))

rnaClusters_all <- readRDS(file = file.path("data", "mrna_protein", "RNA.all.RDS")) %>% 
  pData() %>%
  dplyr::select(Cell.Line.R, starts_with("hierarCluster_pearson_ward.D2"))
#rename the column
rnaClusters_all$RNA_all <- rnaClusters_all$hierarCluster_pearson_ward.D2
rnaClusters_all$hierarCluster_pearson_ward.D2 <- NULL

rnaClusters_B <- readRDS(file = file.path("data", "mrna_protein", "RNA.B.RDS")) %>%
  pData() %>%
  dplyr::select(Cell.Line.R, starts_with("hierarCluster_pearson_ward.D2"))
#rename the column
rnaClusters_B$RNA_B <- rnaClusters_B$hierarCluster_pearson_ward.D2
rnaClusters_B$hierarCluster_pearson_ward.D2 <- NULL

rnaClusters_T <- readRDS(file = file.path("data", "mrna_protein", "RNA.T.RDS")) %>%
  pData() %>%
  dplyr::select(Cell.Line.R, starts_with("hierarCluster_pearson_ward.D2"))
#rename the column
rnaClusters_T$RNA_T <- rnaClusters_T$hierarCluster_pearson_ward.D2
rnaClusters_T$hierarCluster_pearson_ward.D2 <- NULL

allmeta <- proteinsMeta %>%
  left_join(proteinClusters, by = "Cell.Line.R") %>%
  left_join(rnaClusters_all, by = "Cell.Line.R") %>%
  left_join(rnaClusters_B, by = "Cell.Line.R") %>%
  left_join(rnaClusters_T, by = "Cell.Line.R") %>%
  `row.names<-`(.$Cell.Line.R)

# Run this to plot B cells only
fileName <- "B_alluvial_plot.svg"
sankeymeta <- allmeta %>%
  dplyr::filter(Type_Paper %in% c("B-ALL", "BCP-ALL")) %>%
  dplyr::select(Subtype_Paper, proteinCluster = hierarCluster_pearson_ward.D2_B, rnaCluster = RNA_B)


#================
#ggalluvial plotting

stage_to_subtype_and_cluster <- sankeymeta %>%
  group_by(proteinCluster, Subtype_Paper, rnaCluster) %>%
  summarise(connection_count = n())

p <- ggplot(as.data.frame(stage_to_subtype_and_cluster),
            aes(y = connection_count,
                axis1 = proteinCluster, axis2 = Subtype_Paper, axis3 = rnaCluster)) +
  geom_alluvium(aes(fill = as.factor(Subtype_Paper)), alpha = .6) +
  geom_stratum(fill = NA, color = "gray60", width = 0.1) + 
  geom_label(stat = "stratum", aes(label = after_stat(stratum)), size = 1, color = "gray14", label.size = 0) +
  scale_x_discrete(limits = c("Proteomics cluster", "Subtype", "RNA cluster"), expand = c(.05, .05)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 140, 20),  labels = seq(0, 140, 20)) +
  scale_fill_manual(values = color_scheme$subtype.alt) +
  theme(panel.grid = element_blank(),
        line = element_blank()) +
  labs(x = "", y = "", title = '', color = "gray14") +
  guides(fill =  FALSE) +
  theme(axis.line=element_blank(),
        axis.text.y=element_blank())

p

# export to svg
ggsave(
  filename = fileName,
  plot = p,
  device = NULL,
  path = outputFolder,
  scale = 1,
  width = 6,
  height = 7,
  units = c("cm"),
  dpi = 600,
  limitsize = FALSE)

