library(Biobase)
library(tidyverse)
library(Seurat)
library(patchwork)
library(Rsamtools)
library(stringr)
library(viridis)
library(ggplot2)
library(ggthemes)
options(stringsAsFactors = FALSE, scipen = 9999)

source(file.path("functions", "getColorScheme.R"))
source(file.path("functions", "removeNAsFromESet.R"))
source(file.path("figures", "all_theme.R"))

useVersion  <- "publication"
colorscheme <- getColorScheme()

outputFolder <- file.path("figures", "output", "figure_SI03")
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder)

#===================================================================================================================
# PREPROCESSING
#===================================================================================================================

type <- "T-ALL"

# load proteins
proteins <- readRDS(file = file.path("output", useVersion, "consensus_leukemic_clustering", "all", "proteins.RDS"))

# select samples and remove all NA values
proteins_noNA <- proteins %>% 
 # .[, pData(.)$Type == type] %>%
  removeNAsFromESet()

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
# T-ALL ONLY
#===================================================================================================================

imm.seur.tall <- subset(imm.seur, subset = Type == "T-ALL")

tall.lines.toassign <- rownames(imm.seur.tall@meta.data)
assignedstage <- c()
staging.criteria <- c("CD4", "DNTT", "RAG2", "RAG1", "CD5", "SPN", "CD33", "CD1A", "CD2", "HLA-DRB1", "IL2RG","CD3G", "CD3E", "CD7", "LYL1", "TAL1", "CD44") #missing! HOX11, CD8 (Check manually!!!), pre-TCR (PTCRA), CD25, CD10, CD56
staging.model <- list()

#plot staging criteria to identify breakpoints
VlnPlot(imm.seur.tall, features = staging.criteria) #only TALL
VlnPlot(imm.seur, features = staging.criteria, group.by = "Type") #look at whole dataset, with B lineages

seurat.data.mat <- imm.seur.tall@assays$RNA@data
for (i in 1:length(tall.lines.toassign)) {
  staging.model[[tall.lines.toassign[i]]] <- seurat.data.mat[rownames(seurat.data.mat) %in% staging.criteria, colnames(seurat.data.mat) %in% tall.lines.toassign[i]]
  listindx <- staging.model[[tall.lines.toassign[i]]]
  if (max(c(listindx["CD5"]))< 0.3 & max(c(listindx["CD1A"]))< 0.2 & max(c(listindx["CD4"]))< -0.2 & max(c(listindx["CD2"]))< 0  & max(c(listindx["CD7"]))> 0  &  max(c(listindx["CD44"]))> 0) {
    assignedstage[i] <- "pro-T DN"
    cat("done\n")
  } else if (max(c(listindx["LYL1"])) > 0 & max(c(listindx["CD4"]))< -0.2 & max(c(listindx["CD2"]))< 0  &  max(c(listindx["CD7"]))> 0  &  max(c(listindx["CD44"]))> 0) {
    assignedstage[i] <- "LYL1+ DN T-precursor ALL"
    cat("done\n")
  } else if (max(c(listindx["TAL1"])) > 0.5 & max(c(listindx["CD4"]))> -0.2 & max(c(listindx["CD2"]))> 0 & max(c(listindx["CD1A"])) > 0.2 & max(c(listindx["CD44"]))< 0) {
    assignedstage[i] <- "TAL1+ cortical"
    cat("done\n")
  } else if (max(c(listindx["TAL1"]))> 0.5) {
    assignedstage[i] <- "TAL1+ other"
    cat("done\n")
  } else if (max(c(listindx["CD4"]))< -0.2 & (max(c(listindx["CD2"]))> 0 | max(c(listindx["CD5"]))> 0.3)& max(c(listindx["CD7"]))> 0 & max(c(listindx["CD44"]))< 0) {
    assignedstage[i] <- "pre-T DN"
    cat("done\n")
  } else if (max(c(listindx["CD1A"])) > 0.2 & max(c(listindx["CD4"]))> -0.2 & max(c(listindx["CD44"]))< 0) {
    assignedstage[i] <- "TAL1- cortical"
    cat("done\n")
  } else {assignedstage[i] <- "T-lineage other" #ok to bin and confirm T-lineage commitment, since pre-check of markers showed all cell lines T-ALL contain clear high CD3
  cat("unassigned\n")
  } #end if statements
}#end for

#merge to metadata
imm.seur.tall@meta.data[["AssignedStages"]] <- assignedstage

#write out file with assignments for future analyses
write.csv(x = data.frame(cell_line = imm.seur.tall@meta.data$cell_line, 
                     Cell.Line.Name = imm.seur.tall@meta.data$Cell_Line_Name_Paper,
                     Subtype = imm.seur.tall@meta.data$Subtype,
                     Subtype.Alt = imm.seur.tall@meta.data$Subtype.Alt,
                     Subtype_Paper = imm.seur.tall@meta.data$Subtype_Paper,
                     hierarCluster_pearson_ward.D2_T = imm.seur.tall@meta.data$hierarCluster_pearson_ward.D2_T,
                     RNAseq = imm.seur.tall@meta.data$RNAseq,
                     AssignedStages = imm.seur.tall@meta.data$AssignedStages,
                     Type = imm.seur.tall@meta.data$Type),
          file = file.path("meta", "all_meta_reclassified_tall.csv"))

#make them factors
imm.seur.tall@meta.data[["AssignedStagesF"]] <- factor(imm.seur.tall@meta.data[["AssignedStages"]], levels = c("T-lineage other", "pro-T DN", "LYL1+ DN T-precursor ALL", "pre-T DN",
                                                                                                           "TAL1- cortical", "TAL1+ cortical", "TAL1+ other"))

#now, plot violin plots of subtype specific genes
tall.plot <- VlnPlot(imm.seur.tall,
                  features = c("CD1A", "TAL1", "LYL1", "CD4", "CD2", "CD44"),
                  group.by = "AssignedStagesF",
                  combine = F,
                  cols = c("T-lineage other" = "#F3C2A1", "pro-T DN" = "7e3042",
                           "LYL1+ DN T-precursor ALL" = "#A786A3", "pre-T DN" = "#C5C5CF", "TAL1- cortical" ="#f9f49f", 
                           "TAL1+ cortical" = "#f3ea44",  "TAL1+ other" = "#BE8E56"))

yLabel <- expression(paste("relative ", log[2], " protein level"))

tall.plot <- lapply(tall.plot, function (x){
  
  rethemed_plot <- x +
    all_theme() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    theme(axis.text=element_text(size=4))+
    theme(legend.position = "none") +
    ylab(yLabel) +
    xlab("")
  
  
  
  return(rethemed_plot)
  
})

#change color of points/lines from black to grey
for(i in 1:length(tall.plot)){
  
  tall.plot[[i]][["layers"]][[2]][["geom"]][["default_aes"]][["colour"]] <- "grey30"
  tall.plot[[i]][["theme"]][["text"]][["colour"]] <- "grey30"
  tall.plot[[i]][["theme"]][["line"]][["colour"]] <- "grey30"
  
} #end for

(tall.plot[[1]] | tall.plot[[2]] | tall.plot[[3]])/(tall.plot[[4]] | tall.plot[[5]] | tall.plot[[6]])

#export to svg
library(svglite)
outputFolder <- file.path("figures", "output")
for (i in 1:length(tall.plot)) {
  ggsave(filename = paste0("seurat_classification_", i, ".svg"),
         plot = tall.plot[[i]],
         device = NULL,
         path = outputFolder,
         scale = 1,
         width = 6,
         height = 6,
         units = c("cm"),
         dpi = 600,
         limitsize = FALSE)
}

