####################################################################################
####################################################################################
# differential expression analysis
####################################################################################
####################################################################################

# dependencies
library(Biobase)
library(tidyverse)
library(ggrepel)
library(org.Hs.eg.db)
library(DEqMS)


sourceDir <- function(path, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    source(file.path(path, nm), ...)
  }
}
sourceDir(path = "functions")


# constants
name         <- "[Differential expression analysis] "
useVersion   <- "publication"
outputFolder <- file.path("output", useVersion, "de")
colorScheme  <- getColorScheme()


# create maijor output folder, if it doesn't exist
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder)


####################################################################################
message(name, "Load data")
####################################################################################

proteins <- readRDS(file = file.path("output", useVersion, "consensus_leukemic_clustering", "all", "proteins.RDS"))
proteins.B <- readRDS(file = file.path("output", useVersion, "consensus_leukemic_clustering", "B", "proteins.B.RDS"))
proteins.T <- readRDS(file = file.path("output", useVersion, "consensus_leukemic_clustering", "T", "proteins.T.RDS"))


####################################################################################
message(name, "KMT2A-AFF1 vs KMT2A-MLLT1")
####################################################################################

o <- file.path(outputFolder, "KMT2A_AFF1_vs_KMT2A_MLLT1")

eSet <- proteins.B %>%
  .[, pData(.)$Type_Paper == "BCP-ALL"]

deqms <- runDEqMS(eSet = eSet,
                  conditions = as.factor(pData(eSet)$Subtype.Alt),
                  contrasts = "KMT2A.AFF1-KMT2A.MLLT1",
                  outputFolder = o)

volcanoPlot(data = deqms$outputTables$`KMT2A.AFF1-KMT2A.MLLT1`,
            minXLabel = 1.1,
            minYLabel = -log10(.05),
            outputFolder = NULL)


####################################################################################
message(name, "MEF2D vs. other preB")
####################################################################################

o <- file.path(outputFolder, "MEF2D_HNRNPUL1_vs_BCP-ALL")

eSet <- proteins.B %>%
  .[, pData(.)$Type_Paper == "BCP-ALL"]

deqms <- runDEqMS(eSet = eSet,
                  conditions = as.factor(if_else(pData(eSet)$Subtype.Alt == "MEF2D.HNRNPUL1", "MEF2D.HNRNPUL1", "other")),
                  contrasts = "MEF2D.HNRNPUL1-other",
                  outputFolder = o)

volcanoPlot(data = deqms$outputTables$`MEF2D.HNRNPUL1-other`,
            minXLabel = 1.1,
            minYLabel = -log10(.05),
            outputFolder = NULL)
