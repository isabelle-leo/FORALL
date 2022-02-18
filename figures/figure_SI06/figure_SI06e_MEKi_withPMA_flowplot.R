library(RColorBrewer)
library(tidyverse)
library(readxl)
library(viridis)
library(ggthemes)
library(svglite)

library(RColorBrewer)
options(stringsAsFactors = F, scipen = 9999)

source(file.path("figures", "all_theme.R"))
theme_set(all_theme())

# constants
useVersion   <- "publication"
colorScheme  <- getColorScheme()

#Choose output folder
outputFolder <- file.path("figures", "output", "figure_SI06")

# create major output folder, if it doesn't exist
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder)

cell.line.palette <- c("KASUMI 7" = "#E41A1C", "KASUMI 9" = "#377EB8","LC41" = "#4DAF4A", "P30-OKHUBO" = "#984EA3")

# ==============================================

# read in the raw data
raw_counts <- read_xlsx(path = file.path("data", "flow_cytometry", "20210305-flow-data-meki-excel.xlsx"))
raw_counts <- filter(raw_counts, !SubType %in% c("Other"))
raw_counts <- raw_counts[complete.cases(as.numeric(raw_counts$Count)),]
raw_counts$Count <- as.numeric(raw_counts$Count)

#mean and SD
mean_counts <- raw_counts %>%
  group_by(CellType, DrugType, SubType) %>%
  summarise(MeanCount = mean(Count),
            StdDev = sd(Count))
mean_counts <- mutate(mean_counts, SubGroup = DrugType) 
mean_counts$SubGroup <- factor(mean_counts$SubGroup, levels = c("DMSO",
                                                                "1uM Tramatenib",
                                                                "1uM UO126",
                                                                "1uM Selumetinib",
                                                                "1uM ERK11e",
                                                                    "100nM Bryostatin-1",
                                                                "100nM Bryostatin-1, 1uM UO126", 
                                                                "100nM Bryostatin-1, 1uM Tramatenib",
                                                                "100nM Bryostatin-1, 1uM Selumetinib",
                                                                "100nM Bryostatin-1, 1uM Erk11e",
                                                                "25nM PMA",
                                                                "25nM PMA, 1uM UO126",
                                                                    "25nM PMA, 1uM Tramatenib",
                                                                "25nM PMA, 1uM Selumetinib",
                                                                "25nM PMA, 1uM Erk11e"))


# Violin Plots --------------------------------------------------

#define SubGroup for the raw_counts data.frame
raw_counts_plus <- mutate(raw_counts, SubGroup = paste0(DrugType))
raw_counts_plus$SubGroup <- factor(raw_counts_plus$SubGroup, levels = c("DMSO",
                                                                        "1uM Tramatenib",
                                                                        "1uM UO126",
                                                                        "1uM Selumetinib",
                                                                        "1uM ERK11e",
                                                                        "100nM Bryostatin-1",
                                                                        "100nM Bryostatin-1, 1uM UO126", 
                                                                        "100nM Bryostatin-1, 1uM Tramatenib",
                                                                        "100nM Bryostatin-1, 1uM Selumetinib",
                                                                        "100nM Bryostatin-1, 1uM Erk11e",
                                                                        "25nM PMA",
                                                                        "25nM PMA, 1uM UO126",
                                                                        "25nM PMA, 1uM Tramatenib",
                                                                        "25nM PMA, 1uM Selumetinib",
                                                                        "25nM PMA, 1uM Erk11e"))

#plot only MEF2D fusion cell lines
subset_raw_counts_plus <- raw_counts_plus

#fix the labels by subbing strings in the metadata
subset_raw_counts_plus$CellType <- gsub("KAS9", "KASUMI 9", subset_raw_counts_plus$CellType)
subset_raw_counts_plus$CellType <- gsub("KAS7", "KASUMI 7", subset_raw_counts_plus$CellType)

#supplemental plot: PMA

#just pma, normalized to dmso
subset_raw_counts_plus_pma <- subset_raw_counts_plus[subset_raw_counts_plus$SubGroup != "1uM Tramatenib",]
subset_raw_counts_plus_pma <- subset_raw_counts_plus_pma[subset_raw_counts_plus_pma$SubGroup != "1uM UO126",]
subset_raw_counts_plus_pma <- subset_raw_counts_plus_pma[subset_raw_counts_plus_pma$SubGroup != "1uM Selumetinib",]
subset_raw_counts_plus_pma <- subset_raw_counts_plus_pma[subset_raw_counts_plus_pma$SubGroup != "1uM ERK11e",]
subset_raw_counts_plus_pma <- subset_raw_counts_plus_pma[subset_raw_counts_plus_pma$SubGroup != "100nM Bryostatin-1",]
subset_raw_counts_plus_pma <- subset_raw_counts_plus_pma[subset_raw_counts_plus_pma$SubGroup != "100nM Bryostatin-1, 1uM UO126",]
subset_raw_counts_plus_pma <- subset_raw_counts_plus_pma[subset_raw_counts_plus_pma$SubGroup != "100nM Bryostatin-1, 1uM Tramatenib",]
subset_raw_counts_plus_pma <- subset_raw_counts_plus_pma[subset_raw_counts_plus_pma$SubGroup != "100nM Bryostatin-1, 1uM Selumetinib",]
subset_raw_counts_plus_pma <- subset_raw_counts_plus_pma[subset_raw_counts_plus_pma$SubGroup != "100nM Bryostatin-1, 1uM Erk11e",]


experiments <- unique(subset_raw_counts_plus_pma$ExperimentDate)
#make empty data frame for filling
normalized.to.dmso <- subset_raw_counts_plus_pma
normalized.to.dmso$Count <- NULL
#everything will be one by default, make the column
normalized.to.dmso$NormalizedDMSO <- 1

#Fill normalized.to.dmso with normalized values
for (i in 1:length(experiments)){

  one_experiment_subset <- subset_raw_counts_plus_pma[subset_raw_counts_plus_pma$ExperimentDate == experiments[i],]
  one.cell.type <- unique(one_experiment_subset$CellType)
  for (k in 1:length(one.cell.type)){

    one_exp_cell_subset <- one_experiment_subset[one_experiment_subset$CellType == one.cell.type[k],]
    dmso.mean <- mean(one_exp_cell_subset[one_exp_cell_subset$DrugType == "DMSO",]$Count)
    one.experiment.conditions <- unique(one_exp_cell_subset$SubGroup)

    for (j in 1:length(one.experiment.conditions)){

      single.condition <- one_exp_cell_subset[one_exp_cell_subset$SubGroup == one.experiment.conditions[j],]
      normalized.to.dmso[normalized.to.dmso$ExperimentDate == experiments[i] &
                           normalized.to.dmso$SubGroup == one.experiment.conditions[j] &
                           normalized.to.dmso$CellType == one.cell.type[k],]$NormalizedDMSO <- single.condition$Count/dmso.mean
    }#end for

  } #end for

} #end for


p <- ggplot(data = normalized.to.dmso, aes(x = SubGroup, y = NormalizedDMSO), show.legend = F) +
  geom_violin(scale = "width", kernel = "gaussian", linetype = "blank", draw_quantiles = T, show.legend = F, fill = "#CBD5E8") +
  ylab("Viable Treated Cell Count/Viable DMSO Cell Count") +
  xlab("") +
  geom_jitter(size = 4.5, stroke = 0, alpha = 0.8, aes(x = SubGroup, y = NormalizedDMSO, colour = CellType)) +
  scale_color_manual(values = cell.line.palette) + #scale later when there is multiple cell lines to inferno(8) +
  scale_y_continuous(limits = c(0,1.9)) +
  theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1)) +
  theme(legend.position = "bottom")
p

#add p value bars

library("ggpubr")
my_comparisons <- list(c("25nM PMA", "25nM PMA, 1uM UO126"),
                       c("25nM PMA", "25nM PMA, 1uM Tramatenib"),
                       c("25nM PMA", "25nM PMA, 1uM Selumetinib"),
                       c("25nM PMA", "25nM PMA, 1uM Erk11e"),
                       c("25nM PMA", "DMSO"))

p <- p + stat_compare_means(comparisons = my_comparisons, method = "t.test", paired = FALSE)

p

# export to svg
ggsave(
  "flow_suppPMA_meki.svg",
  plot = p,
  device = NULL,
  path = outputFolder,
  scale = 1,
  width = 4,
  height = 8,
  units = c("cm"),
  dpi = 600,
  limitsize = FALSE)
