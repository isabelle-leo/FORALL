library(tidyverse)
library(ggrepel)
library(ggthemes)

source(file.path("functions", "getColorScheme.R"))
source(file.path("figures", "all_theme.R"))
theme_set(all_theme())

# constants
useVersion   <- "production"
colorScheme  <- getColorScheme()

#Choose output folder
outputFolder <- file.path("figures", "output", "figure_04")

# create major output folder, if it doesn't exist
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder)

# load data - latest drug screen results as a 2D matrix of drugs versus protein correlations, indexing relevant pearson corr values
corrMatrix <- get(load(file.path("data", "dsrt", "protein_all_samples_sDSS_pearson_2D_cor_matrix_min_0.75_with_drug_annotation.rda")))
rm(pearson_cor_matrix_with_annotation)

# read meta
meta <- read_tsv(file = file.path("meta", "Metadata_Drugs_DSRT_2020-11-15_528drugs_TargetCorrRankPlot_RJ - Metadata_Drugs_DSRT_2020-11-15.txt")) %>%
  column_to_rownames("...1")

# reorder correlation matrix
corrMatrix <- corrMatrix[, na.exclude(match(row.names(meta), colnames(corrMatrix)))]
corrMatrix <- corrMatrix[rownames(corrMatrix) != "Mechanism.Targets" &
                           rownames(corrMatrix) != "Class.explained" &
                           rownames(corrMatrix) != "Putative.Target.Protein",]
# reorder / filter meta
meta <- meta[colnames(corrMatrix), , drop = FALSE]

# extract genes
genes <- gsub("_ENS*.+", "" ,row.names(corrMatrix))

# get rid of duplicated genes
corrMatrix <- corrMatrix[!duplicated(genes), ]
# reset genes list
genes <- row.names(corrMatrix)

# set the drugs
drugs <- row.names(meta)

# set correlation and data type
corrType <- "pearson"
dataType <- "protein"

# corrType <- "spearman"
# dataType <- "protein"

# drugs of interest
dois <- c("FIMM000227", "FIMM000484", "FIMM000172", "FIMM003775",
          "FIMM003794", "FIMM133798", "FIMM003763", "FIMM136500",
          "FIMM136412", "FIMM136384", "FIMM136531", "FIMM100374",
          "FIMM003790", "FIMM136388", "FIMM000173", "FIMM001355",
          "FIMM136466", "FIMM136498", "FIMM115484", "FIMM136514", 
          "FIMM136513", "FIMM000490")

#dois <- c(NULL)

# target correlation
targetCorrList <- lapply(X = 1:nrow(meta),
                         FUN = function (x) {
                           tmpTargets <- unlist(strsplit(meta$Putative.Target.Protein[x], ";", fixed = TRUE))
                           tmpCorrVec <- as.numeric(as.character(corrMatrix[, x]))
                           tmpTargetCorr <- tmpCorrVec[match(tmpTargets, genes)]
                           tmpTargetCorrRank <- unlist(lapply(tmpTargetCorr, function (y) {
                             return(round(sum(tmpCorrVec <= y, na.rm = TRUE) / length(tmpCorrVec), 2))
                           }))
                           return(cbind(drugs[x], tmpTargets, tmpTargetCorr, tmpTargetCorrRank))
                         })

targetCorrMatrix <- do.call(rbind, targetCorrList)

targetCorrMatrix <- targetCorrMatrix[!is.na(as.numeric(targetCorrMatrix[, 3])), ]
targetCorrMatrix <- targetCorrMatrix[!is.na(as.numeric(targetCorrMatrix[, 4])), ]

targetCorrDf <- data.frame("Drug" = as.character(targetCorrMatrix[, 1]),
                           "Target" = as.character(targetCorrMatrix[, 2]),
                           "Corr" = as.numeric(as.character(targetCorrMatrix[, 3])),
                           "Rank" = as.numeric(as.character(targetCorrMatrix[, 4])),
                           "Drug_class" = as.character(meta$Drug_class[match(targetCorrMatrix[, 1], row.names(meta))]),
                           "Drug_Target" = as.character(paste(gsub("*.+\\.\\.\\.", "", targetCorrMatrix[, 1]), targetCorrMatrix[, 2], sep = " > ")), stringsAsFactors = FALSE) %>%
  mutate(FIMM.ID = Drug %>% gsub("\\.\\.\\..*$", "", .))

# filter for drugs of interest
plotDfDois <- targetCorrDf[!is.na(match(targetCorrDf$FIMM.ID, dois)), ]

# plot
p <- targetCorrDf %>%
  filter(!Drug %in% plotDfDois$Drug) %>%
  ggplot(aes(x = Rank, y = Corr)) +
  geom_point(aes(fill = Drug_class), size = 2, stroke = 0, alpha = 0.5, shape = 21) +
  geom_point(data = targetCorrDf %>% filter(Drug %in% plotDfDois$Drug),
             aes(fill = Drug_class), size = 3, color = "white", stroke = 1, alpha = 1, shape = 21) +
  geom_text_repel(data = plotDfDois,
                  aes(color = Drug_class, label = Drug_Target),
                  show.legend = FALSE,
                  na.rm = TRUE,
                  segment.size = 0.2,
                  parse = TRUE,
                  segment.color = "#AAAAAA",
                  size = 3,
                  point.padding = 0.2,
                  force = 50,
                  max.overlaps = 30) +
  scale_fill_manual(values = colorScheme$Class.explained %>% set_names(gsub("^[A-Z]\\.", "", names(.)))) +
  xlab("Rank") +
  ylab("Spearman correlation") + #Pearson
  labs(color = "Drug class")

#max(spearman_cor_matrix, na.rm=TRUE)

# view plot
p

ggsave(
  "figure_05_drug_target_corr_rank_plot_rev.svg",
  plot = p,
  device = NULL,
  path = outputFolder,
  scale = 1,
  width = 18,
  height = 6,
  units = c("cm"),
  dpi = 600,
  limitsize = FALSE)

library(openxlsx)
#to use for labeling in illustrator
write.xlsx(plotDfDois, col.names=TRUE, row.names=TRUE, asTable = TRUE, file = "plotDfDois.xlsx")
write.table(plotDfDois, file="plotDfDois.txt",row.names=T,quote=F,sep="\t")