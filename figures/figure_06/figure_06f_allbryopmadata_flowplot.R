library(RColorBrewer)
library(tidyverse)
library(readxl)
library(viridis)
library(ggthemes)
library(svglite)
options(stringsAsFactors = F, scipen = 9999)

source(file.path("figures", "all_theme.R"))
theme_set(all_theme())

source(file.path("functions", "getColorScheme.R"))

# constants
useVersion   <- "publication"
colorScheme  <- getColorScheme()

#Choose output folder
outputFolder <- file.path("figures", "output", "figure_06")

# create major output folder, if it doesn't exist
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder)

cell.line.palette <- c("KAS7" = "#E41A1C", "KAS9" = "#377EB8","LC41" = "#4DAF4A", "P30" = "#984EA3", 
                       "ALLPO" = "#80B1D3", "COG355" = "#FDB462", "COG394" = "#B3DE69", "MHHCALL2" = "#FCCDE5","NALM6" = "#D9D9D9",
                       "RCHACV" = "#BC80BD","REH" = "#CCEBC5","SUPB15" = "#FFED6F")

# ==============================================

# read in the raw data
raw_counts <- read_xlsx(path = file.path("data", "flow_cytometry", "20210305-ALL-flow-data-bryopma-excel.xlsx"))

# fix labels
raw_counts <- raw_counts[complete.cases(as.numeric(raw_counts$Count)),]
raw_counts$Count <- as.numeric(raw_counts$Count)
raw_counts$DrugType <- gsub("100nM Bryostatin-1", "BRYO", raw_counts$DrugType)
raw_counts$DrugType <- gsub("25nM PMA", "PMA", raw_counts$DrugType)
raw_counts$CellType <- gsub("P30-OKHUBO", "P30", raw_counts$CellType)

# Violin Plots --------------------------------------------------

#define SubGroup for the raw_counts data.frame
raw_counts_plus <- mutate(raw_counts, SubGroup = paste0(SubType, ".", DrugType))
#factors and order
raw_counts_plus$X2 <- factor(raw_counts_plus$DrugType, levels = c("DMSO","BRYO", "PMA"))
raw_counts_plus$X3 <- factor(raw_counts_plus$SubGroup, levels = c("Mef2d.Hnrnpul1.DMSO",
                                                                  "Mef2d.Hnrnpul1.BRYO", 
                                                                  "Mef2d.Hnrnpul1.PMA",
                                                                  "Other.DMSO",
                                                                  "Other.BRYO", 
                                                                  "Other.PMA"))
raw_counts_plus$Count <- as.numeric(raw_counts_plus$Count)
raw_counts_plus <- raw_counts_plus[complete.cases(raw_counts_plus$Count),]

#make a violin plot of raw_counts_plus

p <- ggplot(data = raw_counts_plus, aes(x = X3, y = Count, fill = SubType), show.legend = F) +
  geom_violin(scale = "width", kernel = "gaussian", linetype = "blank", draw_quantiles = T, show.legend = F) +
  scale_fill_manual(values = c("#CBD5E8", "#FDCDAC")) +
  ylab("Viable cells per 15 uL") +
  xlab("") +
  geom_jitter(size = 1.5, stroke = 0, alpha = 0.8, aes(x = SubGroup, y = Count)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")


#add p value bars

library("ggpubr")
compare_means(Count ~ X3, data = raw_counts_plus, method = "t.test", paired = FALSE)
my_comparisons <- list(c("Other.DMSO", "Other.BRYO"), c("Mef2d.Hnrnpul1.DMSO", "Mef2d.Hnrnpul1.BRYO"),
                       c("Other.DMSO", "Other.PMA"), c("Mef2d.Hnrnpul1.DMSO", "Mef2d.Hnrnpul1.PMA"))

p <- p + stat_compare_means(comparisons = my_comparisons, method = "t.test", paired = FALSE)

p

#----------re-make with dmso normalized data----------

experiments <- unique(raw_counts_plus$ExperimentDate)
#make empty data frame for filling
normalized.to.dmso <- raw_counts_plus
normalized.to.dmso$Count <- NULL
#everything will be one by default, make the column
normalized.to.dmso$NormalizedDMSO <- 1

#Fill normalized.to.dmso with normalized values
for (i in 1:length(experiments)){
  
  one_experiment_subset <- raw_counts_plus[raw_counts_plus$ExperimentDate == experiments[i],]
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


p <- ggplot(data = normalized.to.dmso, aes(x = X3, y = NormalizedDMSO, fill = SubType), show.legend = F) + 
  geom_violin(scale = "width", kernel = "gaussian", linetype = "blank", draw_quantiles = T, show.legend = F) +
  scale_fill_manual(values = c("#CBD5E8", "#FDCDAC")) +
  ylab("Viable Treated Cell Count/Viable DMSO Cell Count") + 
  xlab("") + 
  geom_jitter(size = 1.5, stroke = 0, alpha = 0.8, aes(x = SubGroup, y = NormalizedDMSO)) +
  theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1)) +
  theme(legend.position = "none") 

#add error bars
my_comparisons <- list(c("Other.DMSO", "Other.BRYO"), c("Mef2d.Hnrnpul1.DMSO", "Mef2d.Hnrnpul1.BRYO"),
                       c("Other.DMSO", "Other.PMA"), c("Mef2d.Hnrnpul1.DMSO", "Mef2d.Hnrnpul1.PMA"),
                       c("Other.PMA", "Mef2d.Hnrnpul1.PMA"))

p <- p + stat_compare_means(comparisons = my_comparisons, method = "t.test", paired = FALSE)


p

# export to svg
ggsave(
  "flow_allbryopma.svg",
  plot = p,
  device = NULL,
  path = outputFolder,
  scale = 1,
  width = 9,
  height = 7,
  units = c("cm"),
  dpi = 600,
  limitsize = FALSE)
