library(Biobase)
library(tidyverse)
library(Seurat)
library(patchwork)
library(Rsamtools)
library(stringr)
library(viridis)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
options(stringsAsFactors = FALSE, scipen = 9999)

source(file.path("functions", "getColorScheme.R"))
source(file.path("functions", "removeNAsFromESet.R"))
source(file.path("figures", "all_theme.R"))
theme_set(all_theme())

useVersion  <- "publication"
colorscheme <- getColorScheme()

outputFolder <- file.path("output", useVersion, "figures", "figure_04")
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder)

#===================================================================================================================
# PREPROCESSING
#===================================================================================================================

# load proteins
proteins <- readRDS(file = file.path("output", useVersion, "consensus_leukemic_clustering", "all", "proteins.RDS"))

# remove all NA values
proteins_noNA <- proteins %>% removeNAsFromESet()

# extract meta data
meta <- proteins_noNA %>%
  pData() %>%
  `row.names<-`(.$Cell.Line.R)

# extract a matrix with quantitative values
proteins_mat <- proteins_noNA %>%
  exprs() %>%
  `colnames<-`(meta$Cell.Line.R)

# create the Seurat object
imm.seur <- CreateSeuratObject(counts = proteins_mat,
                               project = "ProtCLCSeurat",
                               meta.data = meta)

#===================================================================================================================
# PRE-B ONLY
#===================================================================================================================
imm.seur.mef2d <- subset(imm.seur, subset = Type != "T-ALL") #exclude T-ALL
#imm.seur.mef2d <- subset(imm.seur, subset = Subtype == "MEF2D.HNRNPUL1") #to just plot mef2d

imm.seur.mef2d@meta.data$Cell.Line.Factor <- gsub(".BR2", "", imm.seur.mef2d@meta.data$Cell.Line.R)
imm.seur.mef2d@meta.data$Logical.Mef2d <- imm.seur.mef2d@meta.data$Subtype  == "MEF2D.HNRNPUL1"
imm.seur.mef2d@meta.data$Logical.Mef2d <- gsub("TRUE","MEF2D.HNRNPUL1",imm.seur.mef2d@meta.data$Logical.Mef2d)
imm.seur.mef2d@meta.data$Logical.Mef2d <- gsub("FALSE","Other",imm.seur.mef2d@meta.data$Logical.Mef2d)

#plot violin plots of subtype specific genes
#MEF2D versus other
plot5a <- VlnPlot(imm.seur.mef2d, features = c("CD19", "DNTT", "VPREB1", "RAG2", "IL7R", "SPN"), group.by = "Logical.Mef2d", combine = F, cols = c(colorscheme[["subtype.alt"]][["MEF2D-HNRNPUL1"]], "gray90"))

#by cell line
#plot5a <- VlnPlot(imm.seur.mef2d, features = c("CD19", "DNTT", "VPREB1", "RAG2", "IL7R", "SPN"), group.by = "Cell.Line.Factor", combine = F, cols = brewer.pal(3, "Pastel1"))
yLabel <- expression(paste("relative ", log[2], " abundance"))

plot5a <- lapply(plot5a, function (x){
  
  rethemed_plot <- x + theme_minimal() + #reset to ALL theme
    theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust=.5)) + theme(legend.position = "none") + ylab(yLabel) + xlab("")
  
  return(rethemed_plot)
  
})

#pdf export
pdf(file="supp_fig_6c.pdf", width=10, height=8)
(plot5a[[1]] | plot5a[[2]] | plot5a[[3]])/(plot5a[[4]] | plot5a[[5]] | plot5a[[6]])
dev.off()

# export to svg
for (i in 1:length(plot5a)) {
  ggsave(filename = paste0("seurat_classification_", i, ".svg"),
         plot = plot5a[[i]],
         device = NULL,
         path = outputFolder,
         scale = 1,
         width = 6,
         height = 5.5,
         units = c("cm"),
         dpi = 600,
         limitsize = FALSE)
}
