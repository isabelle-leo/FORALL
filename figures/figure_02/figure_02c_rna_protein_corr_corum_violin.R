library(Biobase)
library(vioplot)

# dependencies
sourceDir <- function(path, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    source(file.path(path, nm), ...)
  }
}
sourceDir(path = "functions")

# constants
useVersion   <- "publication"
colorScheme  <- getColorScheme()

# Set outputfolder
outputFolder <- file.path("figures", "output", "figure_02")

# constants
name         <- "[Produce figure 2C] "
useVersion   <- "publication"
colorScheme  <- getColorScheme()

# Set outputfolder
outputFolder <- file.path("figures", "output", "figure_02")

# create maijor output folder, if it doesn't exist
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder)

protein_RDS_obj <- file.path("output", useVersion, "consensus_leukemic_clustering", "all", "proteins.RDS")
RNA_RDS_obj <- file.path("data", "mrna_protein", "RNA.RDS")
CORUM_data <- file.path("data", "mrna_protein", "CORUM_coreComplexes_ComplexID_subunits_Genes_name")

correlation <- "spearman"
left_box_color <- "yellow3"
right_box_color <- "steelblue4"

protein_obj <-readRDS(protein_RDS_obj) %>%
  .[, pData(.)$Cell.Line.R != "SEM.BR.NOPS"]
colnames(protein_obj) <- pData(protein_obj)$Cell.Line.R
RNA_obj <-readRDS(RNA_RDS_obj)


NA_n <- apply(exprs(protein_obj), 1, function(x){sum(is.na(x))})
protein_obj <- protein_obj[NA_n == 0, ]
common_genes <- intersect(rownames(RNA_obj), rownames(protein_obj))
common_samples <- intersect(colnames(RNA_obj), colnames(protein_obj))

RNA_obj <- RNA_obj[common_genes, common_samples]
protein_obj <- protein_obj[common_genes, common_samples]

CORUM_data_human <- read.table(CORUM_data, sep = "\t", header = T, stringsAsFactors = F, quote = "")

dim(CORUM_data_human)
CORUM_data_human <- unique(CORUM_data_human)

CORUM_data_human <- CORUM_data_human[!is.na(match(CORUM_data_human$subunits_Genes_name, common_genes)), ]
dim(CORUM_data_human)
number_of_genes_per_complex <- table(CORUM_data_human$ComplexID)
number_of_genes_per_complex <- number_of_genes_per_complex[number_of_genes_per_complex > 1]

CORUM_data_human <- CORUM_data_human[!is.na(match(CORUM_data_human$ComplexID, names(number_of_genes_per_complex))), ]
dim(CORUM_data_human)
rm("number_of_genes_per_complex", "common_genes")

complex_IDs <- sort(unique(CORUM_data_human$ComplexID))

corum_pairs_cor <- function(matrix_1 = "", matrix_2 = ""){
  cor_vec <-  lapply(1:dim(matrix_1)[[1]], function(x){
    mate_1 <- as.numeric(matrix_1[x, ])
    mate_2 <- as.numeric(matrix_2[x, ])
    if((sd(mate_1, na.rm = T)>0) & (sd(mate_2, na.rm = T)>0)){
      temp_cor <- round(cor(x = mate_1, y = mate_2, method = correlation, use = "pairwise.complete.obs"), 3)
    }else{temp_cor <- NA}; return(temp_cor)})
  
  return(do.call("rbind", cor_vec))}

genes_pairs <- lapply(complex_IDs, function(x){
  
  temp_genes_list <- unique(sort(CORUM_data_human$subunits_Genes_name[which(CORUM_data_human$ComplexID == x)]))
  if(length(temp_genes_list) > 2){
    temp_genes_pairs <- expand.grid(temp_genes_list, temp_genes_list)
    temp_genes_pairs <- unique(t(apply(temp_genes_pairs, 1, sort)))
    temp_genes_pairs <- temp_genes_pairs[temp_genes_pairs[, 1] != temp_genes_pairs[, 2], ]
    return(cbind(x, temp_genes_pairs))
  }else{
    return(matrix(c(x, unlist(temp_genes_list)), ncol=3, nrow=1))
  }})
genes_pairs <- do.call("rbind", genes_pairs)
genes_pairs[, 1] <- 0
genes_pairs <- unique(genes_pairs)

for(i in 2:3){genes_pairs <- cbind(genes_pairs,match(genes_pairs[,i],rownames(RNA_obj)))}


samples_type <- pData(protein_obj)$Type
RNA_matrix <- exprs(RNA_obj)
protein_matrix <- exprs(protein_obj)

i <- 1 
while(i !=0){
  rand_1 <- sample(1:dim(protein_matrix)[[1]], dim(genes_pairs)[[1]], replace = T)
  rand_2 <-  sample(1:dim(protein_matrix)[[1]], dim(genes_pairs)[[1]], replace = T)
  i <- sum(rand_1 == rand_2)
}

genes_pairs <- cbind(genes_pairs, rand_1, rand_2)

cor_group <- c("ALL", "B", "T")
regular_expression <- c("B-ALL|preB|T-ALL","B-ALL|preB","T-ALL")

pdf(paste(outputFolder,"/CORUM_complexes_RNA_and_protein_", correlation, "_corr_", gsub("-|:| ","_", Sys.time()), ".pdf", sep =""))
par(mar=c(6,5,4,4))

for( i in 1:3){
  
  RNA_cor  <- corum_pairs_cor(matrix_1 = RNA_matrix[as.numeric(genes_pairs[, 4]), grep(regular_expression[i], samples_type) ], matrix_2 = RNA_matrix[as.numeric(genes_pairs[, 5]), grep(regular_expression[i], samples_type) ])[, 1]
  
  RNA_cor_random  <-  corum_pairs_cor(matrix_1 = RNA_matrix[as.numeric(genes_pairs[, 6]), grep(regular_expression[i], samples_type) ], matrix_2 = RNA_matrix[as.numeric(genes_pairs[, 7]), grep(regular_expression[i], samples_type) ])[, 1]
  
  protein_cor  <- corum_pairs_cor(matrix_1 = protein_matrix[as.numeric(genes_pairs[, 4]), grep(regular_expression[i], samples_type) ], matrix_2 = protein_matrix[as.numeric(genes_pairs[, 5]), grep(regular_expression[i], samples_type) ])[, 1]
  
  protein_cor_random  <- corum_pairs_cor(matrix_1 = protein_matrix[as.numeric(genes_pairs[, 6]), grep(regular_expression[i], samples_type) ], matrix_2 = protein_matrix[as.numeric(genes_pairs[, 7]), grep(regular_expression[i], samples_type) ])[, 1]
  
  plot_df <- data.frame("within_complex_pairs" =  c(protein_cor, RNA_cor),
                        "random_pairs" = c(protein_cor_random, RNA_cor_random),
                        "exp" = as.factor(c(rep("protein", length(RNA_cor)),
                                            rep("RNA", length(protein_cor)))))
  
  # plot_df$group <- relevel(plot_df$group, ref = "random pairs")
  
  vioplot(within_complex_pairs ~ exp, data = plot_df, col = left_box_color, plotCentre = "line", side = "left", ylim = c(-1, 1.2), main = paste("Type:", cor_group[i]),  cex.axis = 1.5, cex.main = 1.5, xlab = NULL, ylab = NULL)
  
  vioplot(random_pairs ~ exp, data = plot_df, col = right_box_color, plotCentre = "line", side = "right", add = T)
  
  legend("bottomleft", fill = c(left_box_color, right_box_color), legend = c("within complex pairs", "random pairs"), title = "CORUM protein complexes", xpd = T, bty = "n",  horiz = TRUE )
  
  # to avoid zero p-value, do t-test for 1000 points only
  random_1000 <- sample(1:20000, 1000, replace = F)
  
  
  mtext(text =  paste("p = ", c(signif(t.test(plot_df$within_complex_pairs[plot_df$exp == "protein"][random_1000] , plot_df$random_pairs[plot_df$exp == "protein"][random_1000] )$p.value,2), signif(t.test(plot_df$within_complex_pairs[plot_df$exp == "RNA"][random_1000] , plot_df$random_pairs[plot_df$exp == "RNA"][random_1000] )$p.value,2))),side = 3, line = -2, at = 1:2, cex = 1)
  mtext(text = paste("n=",as.numeric(table(plot_df$exp)[1]),"paris"),side = 1,line = 3,at = 1.5,cex = 1.5)
  
  
  mtext(text = paste(correlation, "corr."), side = 2, line = 3, cex = 1.5)
}
dev.off()