library(Biobase)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(RColorBrewer)
library(reshape2)
library(ggthemes)
library(svglite)

source(file.path("figures", "all_theme.R"))
theme_set(all_theme())

useVersion <- "publication"

complex_name <- "MCM complex"
complex_name <- "Arp2/3 protein complex"
RNA_RDS <- file.path("data", "mrna_protein", "RNA_batch_corr.RDS") 
protein_RDS <- file.path("output", useVersion, "consensus_leukemic_clustering", "all", "proteins.RDS")

color_start <- "#090979"
color_end <- "#00d4ff"

corum_data <- file.path("data", "mrna_protein", "CORUM_coreComplexes_ComplexID_subunits_Genes_name")
outputFolder <- file.path("figures", "output", "figure_SI02")

plot_function<-function(plot_df,ylabel="",plot_title=""){
  colfunc <- colorRampPalette(c(color_start,color_end))
  ggplot_object<- ggplot(data = plot_df, aes(x = cell, y =  Abundance, colour =  gene)) +
    geom_line(size = 0.6, alpha = .8) +
    scale_colour_manual( values =  colfunc(length(unique(plot_df$gene)))) +
    labs(x = "Cell lines", y = ylabel) +
    ggtitle(plot_title) +
    theme(legend.position = "bottom")
  return(ggplot_object)
}


df_prep<-function(order_vec){
  temp_RNA_matrix_ordered <- RNA_centralized_matrix[,order_vec]
  temp_protein_matrix_ordered <- protein_centralized_matrix[,order_vec]
  
  RNA_df<-melt(temp_RNA_matrix_ordered,varnames = c("gene","cell"),value.name = "Abundance")
  protein_df<-melt(temp_protein_matrix_ordered,varnames = c("gene","cell"),value.name = "Abundance")
  
  not_na<- sort(as.numeric(intersect(rownames(RNA_df[!is.na(RNA_df$Abundance), ]), rownames(protein_df[!is.na(protein_df$Abundance), ]))))
  
  RNA_df_to_plot<-RNA_df[not_na, ]
  protein_df_to_plot<-protein_df[not_na, ]
  
  RNA_df_to_plot$cell<-as.numeric(factor(RNA_df_to_plot$cell,levels = unique(as.character(RNA_df_to_plot$cell))))
  protein_df_to_plot$cell<-as.numeric(factor(protein_df_to_plot$cell,levels = unique(as.character(protein_df_to_plot$cell))))
  return(cbind(protein_df_to_plot, RNA_df_to_plot))
}
# reading data
corum_data<- read.table(corum_data,sep = "\t",stringsAsFactors = F,header = T,quote = "")
RNA_obj<-readRDS(RNA_RDS)
protein_obj<-readRDS(protein_RDS) %>%
  .[, pData(.)$Cell.Line.R != "SEM.BR.NOPS"]
colnames(protein_obj) <- pData(protein_obj)$Cell.Line.R

# prepare expression matrices
# prepare expression matrices
RNA_centralized_matrix <- RNA_obj %>% exprs() %>% + 0.1 %>% log2() %>% t() %>% scale() %>% t()
protein_centralized_matrix <- protein_obj %>% exprs() %>% t() %>% scale() %>% t()


print(identical(colnames(RNA_centralized_matrix) ,colnames(protein_centralized_matrix)))

common_genes <- sort(intersect(rownames(RNA_centralized_matrix),rownames(protein_centralized_matrix)))

RNA_centralized_matrix <- RNA_centralized_matrix[common_genes, ]
protein_centralized_matrix <- protein_centralized_matrix[common_genes, ]

not_na<-apply(protein_centralized_matrix,1,function(x){sum(!is.na(x))})
proteins_passed <- (not_na> (dim(protein_centralized_matrix)[[2]]*0.75))


RNA_centralized_matrix <- RNA_centralized_matrix[proteins_passed, ]
protein_centralized_matrix <- protein_centralized_matrix[proteins_passed, ]

print(identical(rownames(RNA_centralized_matrix), rownames(protein_centralized_matrix)))
print(identical(dim(RNA_centralized_matrix), dim(protein_centralized_matrix)))

temp_genes<-intersect(rownames(RNA_centralized_matrix),corum_data$subunits_Genes_name[grep(paste("^", complex_name, "$",sep = ""), corum_data$ComplexName)])
print(temp_genes)

output_pdf<-paste(outputFolder,"/CORUM_complex_",gsub("-|:|;| |/|\\|","_",complex_name),".pdf",sep = "")
# pdf(output_pdf,width = 15 ,height = 8.5)
RNA_centralized_matrix <- RNA_centralized_matrix[temp_genes, ]
protein_centralized_matrix <- protein_centralized_matrix[temp_genes, ]


# order_vec <- order(unlist(c(apply(RNA_centralized_matrix,2,sum,na.rm=T))))
order_vec <- 1:ncol(RNA_centralized_matrix)
df_to_plot<-df_prep(order_vec =order_vec)
# protein_plot_1<-plot_function(plot_df =df_to_plot[,1:3],ylabel = "Scaled protein abund.",plot_title = "") 
# RNA_plot_1<-plot_function(plot_df =df_to_plot[,4:6],ylabel = "Scaled RNA abund.",plot_title = "")  

# to get mean corr coef and SD
protein_plot_1 <- plot_function(plot_df = df_to_plot[, 1:3],
                                ylabel = "protein",
                                plot_title = paste("mean cor = ",
                                                   round(mean(c(cor(t(protein_centralized_matrix)))), 2),
                                                   ", sd = ", round(sd(c(cor(t(protein_centralized_matrix)))), 2), sep = ""))

RNA_plot_1 <- plot_function(plot_df = df_to_plot[, 4:6],
                            ylabel = "RNA",
                            plot_title = paste("mean cor = ",
                                               round(mean(c(cor(t(RNA_centralized_matrix)))), 2),
                                               ", sd = ", round(sd(c(cor(t(RNA_centralized_matrix)))), 2), sep = ""))

# order_vec <- order(unlist(c(apply(protein_centralized_matrix,2,sum,na.rm=T))))
# df_to_plot<-df_prep(order_vec =order_vec) 
# protein_plot_2<-plot_function(plot_df =df_to_plot[,1:3],ylabel = "protein",plot_title = "ordered by protein")  
# RNA_plot_2<-plot_function(plot_df =df_to_plot[,4:6],ylabel = "RNA",plot_title = "ordered by protein")  

# print(ggarrange(protein_plot_1,RNA_plot_1,protein_plot_2,RNA_plot_2,ncol = 2,nrow = 2))

# dev.off()
# print(paste("The output plots are saved in ","/",output_pdf,sep = ""))

ggsave(filename = "corum_complex_rna_line_plot_mcm.svg",
       plot = RNA_plot_1,
       device = NULL,
       path = outputFolder,
       scale = 1,
       width = 6,
       height = 7,
       units = c("cm"),
       dpi = 600,
       limitsize = FALSE)

ggsave(filename = "corum_complex_protein_line_plot_mcm.svg",
       plot = protein_plot_1,
       device = NULL,
       path = outputFolder,
       scale = 1,
       width = 6,
       height = 7,
       units = c("cm"),
       dpi = 600,
       limitsize = FALSE)


ggsave(filename = "corum_complex_rna_line_plot_arp23.svg",
       plot = RNA_plot_1,
       device = NULL,
       path = outputFolder,
       scale = 1,
       width = 6,
       height = 7,
       units = c("cm"),
       dpi = 600,
       limitsize = FALSE)

ggsave(filename = "corum_complex_protein_line_plot_arp23.svg",
       plot = protein_plot_1,
       device = NULL,
       path = outputFolder,
       scale = 1,
       width = 6,
       height = 7,
       units = c("cm"),
       dpi = 600,
       limitsize = FALSE)

###################################################################################################