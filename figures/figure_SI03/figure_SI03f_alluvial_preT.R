library(Biobase)
library(d3Network)
library(networkD3)
library(tidyverse)
library(ggalluvial)
library(svglite)
options(stringsAsFactors = FALSE, scipen = 9999)

#get functions
sourceDir <- function(path, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    source(file.path(path, nm), ...)
  }
}
sourceDir(path = "functions")

theme_set(all_theme())


useVersion   <- "publication"
color_scheme <- getColorScheme()

# set output folder
outputFolder <- file.path("figures", "output", "figure_SI03")

# create major output folder, if it doesn't exist
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder)


#================
#Get metadata for sankey plot

allmeta <- read.csv(file = file.path("meta", "all_meta_reclassified_tall.csv"), row.names = "Cell.Line.Name") %>% #read in file created during T-ALL  classification script
  dplyr::select(cell_line, Type, AssignedStages, Subtype_Paper)

# get Cell.Line.R
cellLineTranslator <- readRDS(file = file.path("output", useVersion, "proteins.RDS")) %>%
  pData() %>%
  dplyr::select(cell_line, Cell.Line.R) %>%
  dplyr::mutate(cell_line = gsub("\\.0", "", cell_line))

# get cluster information
clusters <- read_tsv(file = file.path("meta", "cell_line_clustering.txt"))

sankeymeta_2 <- allmeta %>%
  left_join(cellLineTranslator, by = "cell_line") %>%
  left_join(clusters, by = "Cell.Line.R") %>%
  filter(Type == "T-ALL") %>%
  dplyr::select(AssignedStages, Subtype_Paper, hierarCluster_pearson_ward.D2_T)


#================
# method 2: ggalluvial


#function to test if they are alluvia format
#is_alluvia_form(stage_to_subtype)

stage_to_subtype_and_cluster <- sankeymeta_2 %>%
  group_by(AssignedStages, Subtype_Paper, hierarCluster_pearson_ward.D2_T) %>%
  summarise(connection_count = n())

color_scheme <- getColorScheme()

color_scheme$AssignedStages <- c("T-lineage other" = "#A786A3", "pro-T DN" = "#C5C5CF",
                                 "LYL1+ DN T-precursor ALL" = "#9EAF9F", "pre-T DN" = "#F3C2A1", "TAL1- cortical" ="#7B7B54", 
                                 "TAL1+ cortical" = "#7B7B54",  "TAL1+ other" = "#BE8E56")

p1 <- ggplot(as.data.frame(stage_to_subtype_and_cluster),
             aes(y = connection_count,
                 axis1 = Subtype_Paper, axis2 = AssignedStages, axis3 = hierarCluster_pearson_ward.D2_T)) +
  geom_alluvium(aes(fill = as.factor(AssignedStages)), alpha = .6) +
  geom_stratum(fill = NA, color = "gray60", width = 0.1) + 
  geom_label(stat = "stratum", aes(label = after_stat(stratum)), size = 1, color = "gray14", label.size = 0) +
  scale_x_discrete(limits = c("Cytogenetic Subtype","Assigned Stage", "Pearson Ward D2 Cluster"), expand = c(.05, .05)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 140, 20),  labels = seq(0, 140, 20)) +
  scale_fill_manual(values = color_scheme$AssignedStages) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        line = element_blank()) +
  labs(x = "", y = "", title = '', color = "gray14") +
  guides(fill =  FALSE) 

p1

# export to svg
ggsave(
  filename = "alluvial.svg",
  plot = p1,
  device = NULL,
  path = outputFolder,
  scale = 1,
  width = 7.5,
  height = 7,
  units = c("cm"),
  dpi = 600,
  limitsize = FALSE)
