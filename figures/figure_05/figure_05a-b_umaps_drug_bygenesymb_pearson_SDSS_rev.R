library(dplyr)
library(Seurat)
library(patchwork)
library(stringr)
library(viridis)
library(ggplot2)
library(readxl)
library(rjson)
library(uwot)
library(Biobase)
library(plyr)
library(ggpubr)
library(ggthemes)
library(extrafont)
options(stringsAsFactors = FALSE, scipen = 9999)

source(file.path("functions", "getColorScheme.R"))
source(file.path("figures", "all_theme.R"))
theme_set(all_theme())

# constants
useVersion   <- "publication"
colorScheme  <- getColorScheme()

#=========
#INITIALIZE
#setwd("/Users/lab/Documents/githubspot/childhood_ALL_cell_lines/")
font_import()
loadfonts()
#==========
#get functions
sourceDir <- function(path, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    source(file.path(path, nm), ...)
  }
}
sourceDir(path = "functions")

colors_all <- getColorScheme()

#Choose output folder
outputFolder <- file.path("figures", "output", "figure_05")

# create major output folder, if it doesn't exist
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder)

#-------------------------------------
#get files

#Newest - B and T both
drugs <- readRDS(file = file.path("output", useVersion, "dsrt", "sDSS.RDS"))
#sDSS above 8
cells_by_drug <- as.data.frame(exprs(drugs))
cells_by_drug$max <- apply(cells_by_drug[, 1:ncol(cells_by_drug)], 1, max) #new col named max with highest sdss per row
cells_by_drug_filt <- cells_by_drug[cells_by_drug$max > 8,]
cells_by_drug_filt$max <- NULL
cells_by_drug_filt <- as.data.frame(t(cells_by_drug_filt))

#read in metadata
drugmeta <- fData(drugs)

#read in pearson correlation per gene symbol: pearson_cor_matrix_with_annotation
#File downloaded from: Luay's analysis, Childhood ALL Microsoft teams: KI.SE\GRP_Childhood ALL - /General/Luay/Data_Analysis_2021/exp_sDSS_cor/
load(file.path("data", "dsrt", "protein_all_samples_sDSS_pearson_2D_cor_matrix_min_0.75_with_drug_annotation.rda"))

#remove annotation
prots_by_drug <- pearson_cor_matrix_with_annotation[rownames(pearson_cor_matrix_with_annotation) != "Mechanism.Targets" &
                                                      rownames(pearson_cor_matrix_with_annotation) != "Class.explained" &
                                                      rownames(pearson_cor_matrix_with_annotation) != "Putative.Target.Protein",]

#no NA values
cells_by_drug_filt <- cells_by_drug_filt[complete.cases(cells_by_drug_filt), ]
prots_by_drug <- prots_by_drug[complete.cases(prots_by_drug), ]

#filter drug metadata by available protein data
drugmeta <- drugmeta[rownames(drugmeta) %in% colnames(cells_by_drug_filt),]
cells_by_drug_filt <- cells_by_drug_filt[,colnames(cells_by_drug_filt) %in% rownames(drugmeta)]
prots_by_drug <- prots_by_drug[,colnames(prots_by_drug) %in% rownames(drugmeta)]
#total - 336 drugs

#put them in the same order
drugmeta <- drugmeta[order(rownames(drugmeta)),]
cells_by_drug_filt <- cells_by_drug_filt[,order(colnames(cells_by_drug_filt))]
prots_by_drug <- prots_by_drug[,order(colnames(prots_by_drug))]

#make seurat.object
seurat.object <- CreateSeuratObject(counts = prots_by_drug,
                                    project = "Drug.Seurat",
                                    meta.data = drugmeta)
#=================================
#Working from Seurat object
#=================================
#can use mvp for finding variable features - mean, dispersion, and bin for each, ID top bin for clustering
#log variance, log mean, loess: use "vst" for selection.method
all.prots <- rownames(seurat.object)
seurat.object <- Seurat::FindVariableFeatures(object = seurat.object, selection.method = "mvp", num.bin = 40, nfeatures = all.prots) 
seurat.object <- ScaleData(seurat.object, features = all.prots) #need this for PCA: center and scale, add to file, keeps raw data


prot.pca <- prcomp(data.matrix(t(seurat.object@assays[["RNA"]]@counts)), scale. = TRUE)
score.df <- as.data.frame(prot.pca$x) *-1 #transform to stay consistent - different programs (prcomp, irlba) make the pca in different rotation directions
score.df$type <- drugmeta[rownames(score.df),"Class.explained"]

#Skree plot
Contributions <- round(summary(prot.pca)$importance[2,] * 100, 2)
barplot(Contributions[1:20], las=3, main="Screeplot Drugs by Pearson, PCA", xlab="component", ylab="Contribution(%)", col="lightblue")
PC1contr <- Contributions[1]
PC2contr <- Contributions[2]

#compare PCA within seurat function
seurat.object <- RunPCA(seurat.object, features = all.prots, npcs = 55) #extra pcs for jack straw
#plot PCA by Seurat to check that their internal math was the same as the other script
DimPlot(seurat.object, reduction = "pca", group.by = "Class.explained", pt.size = 4) + theme_minimal()

#or elbow plot - Seurat function
#pdf(file="elbow_plot_pcs.pdf", width=7, height=6)
ElbowPlot(seurat.object) + theme_minimal()
#dev.off()

#slow!!!! be careful, any more replicates it will take too much time
#jack straw plot - Seurat function
seurat.object <- JackStraw(seurat.object, num.replicate = 100, dims = 55)
seurat.object <- ScoreJackStraw(seurat.object, dims = 1:55)
JackStrawPlot(seurat.object, dims = 1:55)
#use 27

#plot PCA by prcomp {stats}
#pdf(file = file.path(outputFolder, "pca_plot_general.pdf"), width=12, height=8)
ggplot(data = score.df, aes(x = PC1, y = PC2, label = rownames(score.df))) +
  geom_point(aes(color=type), size = 4) +
  xlab(sprintf("PC1 (%s%%)", PC1contr)) + ylab(sprintf("PC2 (%s%%)", PC2contr)) + theme_minimal()
#dev.off()

#umap with uwot package-------------------------------------------------------------
#will plot each drug by cell line

#run on all working drugs
sparse.matrix.filtered <- data.matrix(seurat.object@assays$RNA@data)
set.seed("123")
umap.output <- uwot::umap(data.matrix(t(sparse.matrix.filtered)),
                          n_neighbors = 25, local_connectivity = 1, n_components = 2,
                          repulsion_strength = 1, negative_sample_rate = 5, #can change repulsion strength to weigh low values less than high values
                          spread = 3.5, min_dist = .05, metric = "euclidean",
                          pca = 27, init = "normlaplacian", #SET NUMBER OF PCS IN PCA ARGUMENT
                          init_sdev = .0001, ret_nn = T,
                          nn_method = "annoy", n_trees = 9000, #trees for annoy index, more is better
                          pca_center = F, scale = "z", search_k = 9000*3*80) #perplexity times 3 times n_trees
#18 for sdss 8 cutoff
#nn_umap.output <- umap.output$nn #not actually calculating it, just for n_neighbors in umap
umap.output <- umap.output$embedding
colnames(umap.output) <- paste0("UMAP_", 1:ncol(umap.output))
rownames(umap.output) <- colnames(seurat.object@assays[["RNA"]]@data)
umap.output <- as.matrix(umap.output)

#nice umap - read in the file
saveRDS(umap.output, file.path(outputFolder, "umap_output_drug_by_genesymb.RDS"))
umap.output <- readRDS(file.path(outputFolder, "umap_output_drug_by_genesymb.RDS"))

seurat.object@reductions[["umap"]] <- CreateDimReducObject(embeddings = umap.output, assay = "RNA")

#plot umap with drugs  by drug class
pdf(file=file.path(outputFolder, "overall_point_plot_umap_drugs.pdf"), width=6, height=3)
DimPlot(seurat.object, reduction = "umap", group.by = "Class.explained", pt.size = 1.5, cols = colors_all$Class.explained)
dev.off()

#=====
#Plot drugs by targets
#=====
#=====
#MAKE PLOT WITH ANNOTATED GROUPS BY GREP 
#=====

seurat.object@meta.data$Specific.annotation <- "Other"
for (i in 1:length(seurat.object@meta.data$Specific.annotation)){
  
  if(grepl("ABL",seurat.object@meta.data$Mechanism.Targets[i])){
    seurat.object@meta.data$Specific.annotation[i] <- "ABL inhibitors"
  } #end if
  if(grepl("Bcl",seurat.object@meta.data$Mechanism.Targets[i])){
    seurat.object@meta.data$Specific.annotation[i] <- "Bcl-2 family inhibitors"
  } #end if
  if(grepl("HDAC",seurat.object@meta.data$Mechanism.Targets[i])){
    seurat.object@meta.data$Specific.annotation[i] <- "HDAC inhibitors"
  } #end if
  if(grepl("Mitotic",seurat.object@meta.data$Mechanism.Targets[i])){
    seurat.object@meta.data$Specific.annotation[i] <- "Mitotic Inhibitors"
  } #end if
  if(grepl("PI3K",seurat.object@meta.data$Mechanism.Targets[i])){
    seurat.object@meta.data$Specific.annotation[i] <- "PIK3 family inhibitors"
  } #end if
  if(grepl("TOP", seurat.object@meta.data$Putative.Target.Protein[i])){
    seurat.object@meta.data$Specific.annotation[i] <- "Topoisomerase II inhibitors"
  } #end if
  if(grepl("VEGF",seurat.object@meta.data$Mechanism.Targets[i])){
    seurat.object@meta.data$Specific.annotation[i] <- "VEGFR inhibitors"
  } #end if
  
} #end for

#DimPlot(seurat.object, reduction = "umap", group.by = "Specific.annotation", pt.size = 1.5)

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#plot by grepped target annotation
plot <- ggplot(data = df.plot, aes(x = UMAP_1, y = UMAP_2, colour=Specific.annotation, fill = Specific.annotation), cols = ) +
  geom_point(aes(size = 4, alpha = .5)) +
  scale_color_manual(values = c("Other" = "grey80", "ABL inhibitors" = "#05afa6", "HDAC inhibitors" = "#C01d67", 
                       "Mitotic Inhibitors" = "blue", "PIK3 family inhibitors" = "orange",
                       "Bcl-2 family inhibitors" = "pink", "VEGFR inhibitors" = "#C17af1",
                       "Topoisomerase II inhibitors" = "#873e23")) +
  theme_pubr()

  

seur.plot <- DimPlot(seurat.object, reduction = "umap", group.by = "Specific.annotation", pt.size = 3.5,
                     cols = c("Other" = "grey80", "ABL inhibitors" = "#05afa6", "HDAC inhibitors" = "#C01d67", 
                              "Mitotic Inhibitors" = "blue", "PIK3 family inhibitors" = "orange",
                              "Bcl-2 family inhibitors" = "pink", "VEGFR inhibitors" = "#C17af1",
                              "Topoisomerase II inhibitors" = "#873e23"))

#add specific points to plot, from drugs mentioned in text
drug.from.text <- c("Mubritinib","Tacedinaline","BAY_87-2243", "NMS-873", "Romidepsin", "Panobinostat", "Vorinostat", "Belinostat")

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#Make plot with annotated points
pdf(file=file.path(outputFolder, "point_plot_umap_drugs_labels.pdf"), width=6, height=6)
LabelPoints(seur.plot, points = WhichCells(seurat.object, expression = DRUG.NAME %in% drug.from.text), repel = T, xnudge = .3, ynudge = .05, labels = NULL) +
  geom_point(data = df.plot[df.plot$DRUG.NAME %in% drug.from.text,], aes(x = UMAP_1, y = UMAP_2), size = 2.5, stroke = 2) + theme_pubr() +theme(text = element_text(family = "Arial"))
dev.off()

