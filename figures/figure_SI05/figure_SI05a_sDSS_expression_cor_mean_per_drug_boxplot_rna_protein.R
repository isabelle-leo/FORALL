## Libraries
library(tidyverse)
library(ggrepel)
library(ggpubr)
library(ggthemes)
library(svglite)

# dependencies
sourceDir <- function(path, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    source(file.path(path, nm), ...)
  }
}
sourceDir(path = "functions")

source(file.path("functions", "getColorScheme.R"))
source(file.path("figures", "all_theme.R"))
theme_set(all_theme())

# constants
useVersion   <- "publication"
colorScheme  <- getColorScheme()

# Set outputfolder
outputFolder <- file.path("figures", "output", "figure_SI05")

# create maijor output folder, if it doesn't exist
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder)

###################################################################################################
########################----- sDSS expression cor mean per drug boxplot -----######################
###################################################################################################

RNA_box_color <- "#cecc08"
protein_box_color <- "#36658b"

load(file.path("data", "dsrt", "drugs_matched_RNA_protein_cor_2021_03_23.rda"))

mean_list <- lapply(drugs_matched_RNA_protein_cor,function(x){
  x1 <- apply(x[,2:5],2,function(y){
    y <- as.numeric(y)
    y1 <- round(mean(abs(y),na.rm=T),4)
    y2 <- round(mean(y[y>0],na.rm=T),4)
    y3 <- round(mean(y[y<0],na.rm=T),4)
    return(c(y1,y2,y3))
  })
  return(c(x1,as.character(x[1,6])))
})

mean_matrix <- do.call("rbind",mean_list)

abs_mean_df <- data.frame("Expression" = as.character(rep(c("RNA", "protein"), each = dim(mean_matrix)[[1]])),
                          "pearson_all_abs" = c(as.numeric(mean_matrix[,1]), as.numeric(mean_matrix[,7])),
                          "pearson_pos" = c(as.numeric(mean_matrix[,2]), as.numeric(mean_matrix[,8])),
                          "pearson_neg" = c(as.numeric(mean_matrix[,3]), as.numeric(mean_matrix[,9])),
                          "spearman_all_abs" = c(as.numeric(mean_matrix[,4]), as.numeric(mean_matrix[,10])),
                          "spearman_pos" = c(as.numeric(mean_matrix[,5]), as.numeric(mean_matrix[,11])),
                          "spearman_neg" = c(as.numeric(mean_matrix[,6]), as.numeric(mean_matrix[,12])),stringsAsFactors = T)

# postscript(paste(outputFolder,"/","sDSS_exp_cor_mean_per_drug_boxplot_", gsub(":| |-", "_", Sys.time()), ".ps", sep = ""), width = 8.3, height = 11.7)
i <- "all"
j <- "pearson"

# for( i in c("all", "pos", "neg")){
#   for( j in c("spearman", "pearson")){

temp_df <- abs_mean_df[ , c(1, grep(paste(j, i, sep = "_"), colnames(abs_mean_df)))]

colnames(temp_df)[2] <- "corr"
temp_df <- temp_df[ !is.na(temp_df$corr), ]
ttest_pvalue <- signif(t.test(corr ~ Expression, data = temp_df, paired = T)$p.value, 2)

if(i=="all"){ylabel=paste("Average abs.",j,"corr.")}else{ylabel=paste("Average ",i,". ",j," corr.",sep = "")}

p <- temp_df %>%
  ggplot(aes(x = Expression, y = corr, color = Expression)) +
  geom_jitter(size = 1, stroke = 0, alpha = 0.8) +
  geom_violin(alpha = 0.1, color = "#222222") +
  stat_compare_means(method = "t.test", paired = TRUE, size = 5) +
  scale_color_manual(values = c("RNA" = RNA_box_color, "protein" = protein_box_color)) +
  xlab("") +
  ylab("Pearson correlation") +
  theme(legend.position = "none")

#view the plot
p
#   }
# }
# dev.off()

# export to svg
ggsave(
  "figure_SI05_sDSS_expression_cor_mean_per_drug_boxplot_rna_protein_rev.svg",
  plot = p,
  device = NULL,
  path = outputFolder,
  scale = 1,
  width = 5,
  height = 5,
  units = c("cm"),
  dpi = 600,
  limitsize = FALSE)

