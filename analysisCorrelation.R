####################################################################################
####################################################################################
# correlation analysis
####################################################################################
####################################################################################

library(Biobase)
library(tidyverse)
library(Hmisc)
library(ggthemes)

theme_set(theme_tufte())

sourceDir <- function(path, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    source(file.path(path, nm), ...)
  }
}
sourceDir(path = "functions")


# constants
name          <- "[Correlation analysis] "
useVersion    <- "2020-09-02"
outputFolder  <- file.path("output", useVersion, "correlation")
colorScheme   <- getColorScheme()

if (!is.null(outputFolder) && !dir.exists(outputFolder)) {
  dir.create(outputFolder, recursive = TRUE)
  dir.create(file.path(outputFolder, "all"))
  dir.create(file.path(outputFolder, "B"))
  dir.create(file.path(outputFolder, "T"))
}


####################################################################################
message(name, "Load data")
####################################################################################

proteins <- readRDS(file = file.path("output", useVersion, "consensus_leukemic_clustering", "all", "proteins.RDS"))
proteins_B <- readRDS(file = file.path("output", useVersion, "consensus_leukemic_clustering", "B", "proteins.B.RDS"))
proteins_T <- readRDS(file = file.path("output", useVersion, "consensus_leukemic_clustering", "T", "proteins.T.RDS"))


####################################################################################
message(name, "Calculate correlations")
####################################################################################

# calculate correlations
proteins_corr <- correlate(mat = proteins %>% exprs(), type = "pearson")
write_tsv(x = proteins_corr$df,
          file = file.path(outputFolder, "all", "correlations_df.txt"))
write_tsv(x = as.data.frame(proteins_corr$mat),
          file = file.path(outputFolder, "all", "correlations_mat.txt"))
rm(proteins_corr)

proteins_corr_B <- correlate(mat = proteins_B %>% exprs(), type = "pearson")
write_tsv(x = proteins_corr_B$df,
          file = file.path(outputFolder, "B", "correlations_df.txt"))
write_tsv(x = as.data.frame(proteins_corr_B$mat),
          file = file.path(outputFolder, "B", "correlations_mat.txt"))
rm(proteins_corr_B)

proteins_corr_T <- correlate(mat = proteins_T %>% exprs(), type = "pearson")
write_tsv(x = proteins_corr_T$df,
          file = file.path(outputFolder, "T", "correlations_df.txt"))
write_tsv(x = as.data.frame(proteins_corr_T$mat),
          file = file.path(outputFolder, "T", "correlations_mat.txt"))
rm(proteins_corr_T)
