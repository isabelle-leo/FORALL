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

type <- "preB"

# load proteins
proteins <- readRDS(file = file.path("output", useVersion, "consensus_leukemic_clustering", "all", "proteins.RDS"))

# select samples and remove all NA values
proteins_noNA <- proteins %>% 
  .[, pData(.)$Type == type] %>%
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
# PRE-B ONLY
#===================================================================================================================

imm.seurpb <- imm.seur
# imm.seurpb <- subset(imm.seur, subset = Type == "preB")
# imm.seurpb <- subset(imm.seur, subset = Type != "T-ALL") #To run with B-ALL included

preb.lines.toassign <- rownames(imm.seurpb@meta.data)
assignedstage <- c()
staging.criteria <- c("CD19", "DNTT", "VPREB1", "RAG2", "RAG1", "IGHM")
staging.model <- list()

seurat.data.mat <- imm.seurpb@assays$RNA@data
for (i in 1:length(preb.lines.toassign)) {
  staging.model[[preb.lines.toassign[i]]] <- seurat.data.mat[rownames(seurat.data.mat) %in% staging.criteria, colnames(seurat.data.mat) %in% preb.lines.toassign[i]]
  listindx <- staging.model[[preb.lines.toassign[i]]]
  if (max(c(listindx["IGHM"]))> 1) {
    assignedstage[i] <- "immature B"
    cat("done\n")
  } else if (listindx["DNTT"] < 0 & max(c(listindx["RAG1"], listindx["RAG2"])) > -.05 & max(c(listindx["IGHM"]))< 1) {
    assignedstage[i] <- "late pre-B"
    cat("done\n")
  } else if (listindx["DNTT"] > 0 & listindx["CD19"] > -.5 & max(c(listindx["RAG1"], listindx["RAG2"])) > -.05) {
    assignedstage[i] <- "pro-B"
    cat("done\n")
  } else if (listindx["DNTT"] > 0 & listindx["CD19"] < -.5 & max(c(listindx["RAG1"], listindx["RAG2"])) > -.05 & listindx["VPREB1"] < 0) {
    assignedstage[i] <- "pre-pro-B"
    cat("done\n")
  } else if (max(c(listindx["RAG1"], listindx["RAG2"])) < 0 & listindx["VPREB1"] >0) {
    assignedstage[i] <- "early pre-B"
    cat("done\n")
  } else {assignedstage[i] <- "pre-B other"
  cat("unassigned\n")
  } #end if statements
}#end for

#merge to metadata
imm.seurpb@meta.data[["AssignedStages"]] <- assignedstage

#write out file with assignments for future analyses
write_csv(x = data.frame(cell_line = imm.seurpb@meta.data$cell_line, 
                     Cell.Line.Name = imm.seurpb@meta.data$Cell_Line_Name_Paper,
                     Subtype = imm.seurpb@meta.data$Subtype,
                     Subtype.Alt = imm.seurpb@meta.data$Subtype.Alt,
                     Subtype_Paper = imm.seurpb@meta.data$Subtype_Paper,
                     hierarCluster_pearson_ward.D2_B = imm.seurpb@meta.data$hierarCluster_pearson_ward.D2_B,
                     RNAseq = imm.seurpb@meta.data$RNAseq,
                     AssignedStages = imm.seurpb@meta.data$AssignedStages,
                     Type = imm.seurpb@meta.data$Type),
          file = file.path("meta", "all_meta_reclassified.csv"))

#make them factors
imm.seurpb@meta.data[["AssignedStagesF"]] <- factor(imm.seurpb@meta.data[["AssignedStages"]], levels = c("pre-B other", "pre-pro-B",
                                                                                                         "pro-B", "early pre-B", "late pre-B", "immature B"))

#now, plot violin plots of subtype specific genes
plot_si03 <- VlnPlot(imm.seurpb,
                  features = c("CD19", "DNTT", "VPREB1", "RAG2", "IL7R", "SPN"),
                  group.by = "AssignedStagesF",
                  combine = F,
                  cols = c("#fb8072", "#80b1d3", "#fdb462", "#fccde5","#b3b300", "gray20"))

yLabel <- expression(paste("relative ", log[2], " protein level"))

plot_si03 <- lapply(plot_si03, function (x){
  
  rethemed_plot <- x +
    all_theme() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    theme(legend.position = "none") +
    ylab(yLabel) +
    xlab("")
  
  
  
  return(rethemed_plot)
  
})

#change color of points/lines from black to grey
for(i in 1:length(plot_si03)){
  
  plot_si03[[i]][["layers"]][[2]][["geom"]][["default_aes"]][["colour"]] <- "grey30"
  plot_si03[[i]][["theme"]][["text"]][["colour"]] <- "grey30"
  plot_si03[[i]][["theme"]][["line"]][["colour"]] <- "grey30"
  
} #end for

(plot_si03[[1]] | plot_si03[[2]] | plot_si03[[3]])/(plot_si03[[4]] | plot_si03[[5]] | plot_si03[[6]])

# export to svg
for (i in 1:length(plot_si03)) {
  ggsave(filename = paste0("seurat_classification_", i, ".svg"),
         plot = plot_si03[[i]],
         device = NULL,
         path = outputFolder,
         scale = 1,
         width = 4.5,
         height = 5,
         units = c("cm"),
         dpi = 600,
         limitsize = FALSE)
}
