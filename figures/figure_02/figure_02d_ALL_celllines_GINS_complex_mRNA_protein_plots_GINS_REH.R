### ALL GINS complex


## Libraries
library(Biobase)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(openxlsx)


# Functions
plot.predicted <- function(dt1, dt2) {
  common_samples <- names(dt1)
  shared_samples <- common_samples[!(is.na(dt1) | is.na(dt2))]
  non_shared_samples <- setdiff(common_samples, shared_samples)
  
  val1 <- as.numeric(dt1)
  val2 <-  as.numeric(dt2)
  
  
  mod <- lm(val1 ~ val2)
  
  modpred <- predict(mod)
  modpred <- c(modpred, rep(NA, length(non_shared_samples)))
  names(modpred) <- c(shared_samples, non_shared_samples)
  
  
  modresiduals <- residuals(mod)
  modresiduals <- c(modresiduals, rep(NA, length(non_shared_samples)))
  names(modresiduals) <- c(shared_samples, non_shared_samples)
  return(list('residuals' = modresiduals[common_samples],
              'predictions' = modpred))
  
}

## Data
proteomics <- readRDS(file.path("output", useVersion, "consensus_leukemic_clustering", "all", "proteins.RDS"))
transcriptomics <- readRDS('data/mrna_protein/RNA.RDS')

proteomics_pdat <- pData(proteomics)
protein_table <- exprs(proteomics)

transcriptomics_pdat <- pData(transcriptomics)
rna_table <- exprs(transcriptomics)
colnames(rna_table) <- transcriptomics_pdat$Proteomics.sample

## CCLE Mutation
mutation <- read.delim('data/mrna_protein/ALL_cell_line_depmap_mutation_non_silent_all.txt', sep= '\t')
mutation <- mutation[grepl('GINS', mutation$Hugo_Symbol) & mutation$DepMap_ID %in%  proteomics_pdat$DepMap_ID & mutation$Variant_Classification == 'Missense_Mutation', ]

proteomics_pdat$mutation_GINS <- sapply(proteomics_pdat$DepMap_ID, function(i) {
  idx <- which(mutation$DepMap_ID %in% i)
  paste(paste(mutation$Hugo_Symbol[idx], mutation$Protein_Change[idx], sep = '_'), collapse = ', ')})

## Filter
## Remove all NAs
# toremove <- which(apply(rna_table, 2, function(i) all(is.na(i))))
# rna_table <- rna_table[, -toremove]


# Log2 transformation
rna_table <- log2(rna_table + 1)

common_samples <- Reduce(intersect, list(colnames(protein_table), unlist(strsplit(colnames(rna_table), ';')),
                                         proteomics_pdat$proteomics_id[proteomics_pdat$Type_Paper != 'EBV']))

proteomics_pdat <- proteomics_pdat[match(common_samples,proteomics_pdat$proteomics_id),  ]

rna_table <- rna_table[, sapply(common_samples, function(i) grep(i, colnames(rna_table))) ]
protein_table <- protein_table[, common_samples]

## RNA - protein
gins_complex <- c('GINS1', 'GINS2', 'GINS3', 'GINS4')

res_rna_protein  = lapply(gins_complex, function(i) {
  
  predicted_val <- plot.predicted( as.numeric(protein_table[i, ]), as.numeric(rna_table[i, ]))$predictions
  
  datPlot <- data.frame('mRNA' = as.numeric(rna_table[i, ]),
                        'Protein' =   as.numeric(protein_table[i, ]),
                        'Predicted' = predicted_val, 
                        'cellline' = proteomics_pdat$Cell.Line.R, 
                        'mutation' = ifelse(proteomics_pdat$mutation_GINS == '','', paste(proteomics_pdat$Cell.Line.R,
                                                                                          proteomics_pdat$mutation_GINS, sep = ':'))
  )
  
  cor_p <- cor.test(datPlot$mRNA, datPlot$Protein, method = 'spearman')
  cor_p$p.value <- ifelse(cor_p$p.value == 0, '< 2.2e-16', paste0('= ', format(cor_p$p.value, digits = 2, scientific = TRUE)))
  
  p <- ggplot(data = datPlot, mapping = aes(x = mRNA, y = Protein)) +
    geom_smooth(method = 'lm' , se = FALSE, col = '#D0CD4A', size = 0.2) + 
    geom_point(col = 'black', size = 0.05) + 
    geom_point(data = subset(datPlot, mutation != ''), mapping = aes(fill = mutation), size = 0.4, col = '#2A638B') + 
    geom_segment(data = subset(datPlot, mutation != ''), aes(xend = mRNA, yend = Predicted), lty = 2, col = '#2A638B', size = 0.2) +
    geom_label_repel(data = subset(datPlot, mutation != ''), label = c(4,2,1,3), force = 10, col = '#2A638B', size = 1, label.size = 0.1, label.padding = 0.1, box.padding = 0.1, nudge_y = 0.05) +
    scale_fill_manual(values = rep('black', 4), labels = sort(paste(c(1,2,3,4), sort(subset(datPlot, mutation != '')$mutation)))) +
    annotate(geom = 'text', x = -Inf, y= Inf,  label = paste('Spearman\'s rho =', round(cor_p$estimate,digits = 2)),  vjust = 1.2, hjust = -0.01, size = 3/2.54) +
    annotate(geom = 'text', x = -Inf, y= Inf,  label = paste0('p-value ', cor_p$p.value),  vjust = 3.2, hjust = -0.02,  size = 3/2.54) + 
    annotate(geom = 'text', x = -Inf, y= Inf,  label = paste0('n = ', nrow(datPlot)),  vjust = 5.2, hjust = -0.04,  size = 3/2.54) + 
    theme_bw(base_size=4) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          legend.position = 'bottom',
          legend.text = element_text(size=2.4)) +
    labs(x= 'RNA - log2 TPM', y = 'Protein - log2 ratio', title  = i) + 
    theme(plot.title = element_text(hjust = 0, face = 'bold'))  + 
    guides(fill = guide_legend(title = '', ncol = 2, override.aes = list(col = NA), 
                               keywidth = 0.05,
                               keyheight = 0.05,
                               default.unit="inch"))
  
  colnames(datPlot) <- paste(i, colnames(datPlot), sep = '_')
  
  list(p, datPlot)
  
  # pdf(paste0('Figures/RNA_protein_correlation_GINScomplex_', i, '.pdf'), width = 5, height = 6)
  # plot(p)
  # dev.off()
  
})

# saveRDS(res_rna_protein, 'Figures/RNA_protein_correlation_GINScomplex')

pdf('Figures/RNA_protein_correlation_GINScomplex_merged.pdf', width = 8/2.54, height = 8/2.54)
p <- res_rna_protein[[1]][[1]] + res_rna_protein[[2]][[1]] + res_rna_protein[[3]][[1]] + res_rna_protein[[4]][[1]] + plot_annotation(tag_levels = 'a')  & theme(plot.tag = element_text(face ='bold'))
plot(p)
dev.off()

# Source data 
source_dat <- do.call(cbind, lapply(1:length(gins_complex), function(i) res_rna_protein[[i]][[2]]))


# Remove superfluous columns
source_dat <- source_dat[, -grep('_cellline', colnames(source_dat))[-1]]
source_dat <- source_dat[, -grep('_mutation', colnames(source_dat))[-1]]
source_dat <- source_dat[, -grep('_Predicted', colnames(source_dat))]


source_dat <- source_dat[, c('GINS1_cellline', 'GINS1_mutation', setdiff(colnames(source_dat),  c('GINS1_cellline', 'GINS1_mutation')))]

colnames(source_dat)[1:2] <- c('Cell_line', 'GINS_complex_mutation')

wb <- createWorkbook()
modifyBaseFont(wb, fontSize = 12, fontName = "Calibri")

addWorksheet(wb, sheetName = "Figure 2d", gridLines = TRUE)
writeData(wb, sheet = "Figure 2d", 'Figure 2d', startCol = "A", startRow = 1, colNames = FALSE, rowNames = FALSE,  keepNA = TRUE, na.string = 'NA')
writeData(wb, sheet = "Figure 2d", source_dat, startCol = "A", startRow = 2, colNames = TRUE, rowNames = FALSE, keepNA = TRUE, na.string = 'NA')
setColWidths(wb, sheet = "Figure 2d", cols = 1:(ncol(source_dat)), widths = 15)

# Save the source data as an excel document
#saveWorkbook(wb,  'Data/Source_data_figure2D.xlsx', overwrite = TRUE) ### 

###