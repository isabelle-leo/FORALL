####################################################################################
####################################################################################
# drug-sensitivity correlations
####################################################################################
####################################################################################

# dependencies
library(Biobase)
library(tidyverse)

theme_set(theme_tufte())

useVersion <- "publication"
outputFolder <- file.path("output", useVersion, "dsrt")

####################### LOAD AND PREPARE DATA ####################### 

proteins <- readRDS(file = file.path("output", useVersion, "proteins.RDS"))

clusters <- read_tsv(file = file.path("meta", "cell_line_clustering.txt"))

assay_data <- read_delim(file = file.path("data", "dsrt", "sDSS_2D_matrix_27_Cell_Lines_RJ_2020-12-18.txt"),
                         delim = "\t",
                         quote = "",
                         col_names = TRUE) %>%
  column_to_rownames("FIMM.ID") %>%
  as.matrix()

drugs_meta <- read_delim(file = file.path("meta", "Metadata_Drugs_DSRT_2020-11-15_528drugs_TargetCorrRankPlot_RJ - Metadata_Drugs_DSRT_2020-11-15.txt"),
                         delim = "\t",
                         quote = "",
                         col_names = TRUE)

feature_data <- drugs_meta %>%
  mutate(rownames = FIMM.ID) %>%
  column_to_rownames("rownames") %>%
  .[row.names(assay_data), ]

pheno_data <- proteins %>%
  pData() %>%
  distinct(Cell.Line.R, .keep_all = TRUE) %>%
  dplyr::select(-starts_with("hierarCluster")) %>%
  left_join(clusters, by = "Cell.Line.R") %>%
  mutate(rownames = Cell.Line.R) %>%
  column_to_rownames("rownames") %>%
  .[colnames(assay_data), ]

drugsens <- ExpressionSet(assayData = assay_data,
                          phenoData = AnnotatedDataFrame(pheno_data),
                          featureData = AnnotatedDataFrame(feature_data))

# filter out TANOUE and MN60
drugsensFilt <- drugsens %>%
  .[, !pData(.)$Cell.Line.R %in% c("TANOUE", "MN.60")]

saveRDS (object = drugsensFilt,
         file = file.path(outputFolder, "drugsens.RDS"))
