####################################################################################
####################################################################################
# initialize constants and eSets
####################################################################################
####################################################################################

# dependencies
library(data.table)
library(readxl)
library(tidyverse)
library(matrixStats)
library(Biobase)

sourceDir <- function(path, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    source(file.path(path, nm), ...)
  }
}
sourceDir(path = "functions")


# constants
name                  <- "[Initialisation] "
version               <- "publication"
outputFolder          <- file.path("output", version)

proteinsFileName      <- file.path("data", "12683_STD_ALL_Cell_Line_Panel_Set1-10_BothIPG_Nextflow_ENS99_2020-05-16", "symbols_table.txt") #NOTE: Obtained from NextFlow pipeline, tidy version is available in Supplementary Data 3
peptidesFileName      <- file.path("data", "12683_STD_ALL_Cell_Line_Panel_Set1-10_BothIPG_Nextflow_ENS99_2020-05-16", "peptides_table.txt") #NOTE: Obtained from NextFlow pipeline
metaFileName          <- "meta/ALL cell lines Scilife_CCK_CMB_2020-12-10.xlsx" #NOTE: in the paper, can be obtained from Supplementary Data 1
metaSheetName         <- "Commercial ALL cell lines"
subcellLocFile        <- file.path("meta", "SCInfo.txt")
geneAnnotationFolder  <- file.path("meta", "gene_annotation")

tmtString             <- "_tmt10plex_"
proteinQuantBaseRegex <- paste0("Set[0-9]+", tmtString, "1[2,3][0,1,6,7,8,9][N,C]?")
proteinQuantRegex     <- paste0(proteinQuantBaseRegex, "$")

# create output folder
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder, recursive = TRUE)

####################################################################################
message(name, "Make protein and peptide ExpressionSets")
####################################################################################

# load the gene-symbol centric table
proteins.raw <- makeExpressionSet(expressionsFileName = proteinsFileName,
                                  type = "protein",
                                  mode = "gene_symbol_centric",
                                  quantRegex = proteinQuantRegex,
                                  geneSymbolColumn = "Gene Name",
                                  metaFileName = metaFileName,
                                  metaSheetName = metaSheetName,
                                  tmtString = tmtString)

# load the peptides
peptides.raw <- makeExpressionSet(expressionsFileName = peptidesFileName,
                                  type = "peptide",
                                  mode = "gene_symbol_centric",
                                  quantRegex = proteinQuantRegex,
                                  metaFileName = metaFileName,
                                  metaSheetName = metaSheetName,
                                  tmtString = tmtString)


####################################################################################
message(name, "Tidy up the ExpressionSets")
####################################################################################

proteins <- tidyExpressionSet(expressions.raw = proteins.raw,
                              idCol = "cell_line",
                              poolName = "__POOL__",
                              doLog2 = FALSE)

peptides <- tidyExpressionSet(expressions.raw = peptides.raw,
                              idCol = "cell_line",
                              poolName = "__POOL__",
                              doLog2 = FALSE)


####################################################################################
message(name, "Annotate the ExpressionSets")
####################################################################################

proteins <- annotateExpressionSet(expressions = proteins,
                                  subcellLocFile = subcellLocFile,
                                  geneAnnotationFolder = geneAnnotationFolder)

peptides <- annotateExpressionSet(expressions = peptides,
                                  subcellLocFile = subcellLocFile,
                                  geneAnnotationFolder = geneAnnotationFolder)


####################################################################################
message(name, "Customisation")
####################################################################################

# split off the SEM cell lines for future analysis
sems <- proteins %>%
  .[, pData(.)$cell_line == "SEM"]

# handle the replicates and take SEM median
proteins <- customise(expressions = proteins)

peptides <- customise(expressions = peptides)


####################################################################################
message(name, "Store the proteins and peptides objects as RDS")
####################################################################################

saveRDS(object = proteins,
        file = file.path(outputFolder, "proteins.RDS"))

saveRDS(object = sems,
        file = file.path(outputFolder, "sems.RDS"))

saveRDS(object = peptides,
        file = file.path(outputFolder, "peptides.RDS"))


####################################################################################
message(name, "Store the proteins as TXT")
####################################################################################

saveESetAsTxt(eSet = proteins,
              file = file.path(outputFolder, "proteins_all.txt"),
              colNameColumn = "Cell.Line.R",
              rowNameColumn = "gene_symbol",
              additionalRows = c("Type", "Gender", "Tissue", "Age", "Subtype"),
              additionalColumns = NULL,
              extendedColnames = FALSE)

saveESetAsTxt(eSet = removeNAsFromESet(proteins, na_ratio = 0),
              file = file.path(outputFolder, "proteins_full_overlap.txt"),
              colNameColumn = "Cell.Line.R",
              rowNameColumn = "gene_symbol",
              additionalRows = c("Type", "Gender", "Tissue", "Age", "Subtype"),
              additionalColumns = NULL,
              extendedColnames = FALSE)

saveESetAsTxt(eSet = proteins[, grepl("B", pData(proteins)$Type)] %>% removeNAsFromESet(na_ratio = 0),
              file = file.path(outputFolder, "proteins_B_full_overlap.txt"),
              colNameColumn = "Cell.Line.R",
              rowNameColumn = "gene_symbol",
              additionalRows = c("Type", "Gender", "Tissue", "Age", "Subtype"),
              additionalColumns = NULL,
              extendedColnames = FALSE)

saveESetAsTxt(eSet = proteins[, grepl("T", pData(proteins)$Type)] %>% removeNAsFromESet(na_ratio = 0),
              file = file.path(outputFolder, "proteins_T_full_overlap.txt"),
              colNameColumn = "Cell.Line.R",
              rowNameColumn = "gene_symbol",
              additionalRows = c("Type", "Gender", "Tissue", "Age", "Subtype"),
              additionalColumns = NULL,
              extendedColnames = FALSE)


####################################################################################
message(name, "Clear environment")
####################################################################################

rm(list = ls())
