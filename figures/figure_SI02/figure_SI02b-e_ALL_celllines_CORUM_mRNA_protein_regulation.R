### ALL CORUM post transcriptional regulation

## Libraries
library(reshape2)
library(DescTools)
library(dplyr)
library(ggplot2)
library(ggsignif)
library(GSA)
library(openxlsx)
library(data.table)



## Functions
# Find mean, sum and label position for each group 
boxplot.top <- function(value) {
  thres <- quantile(value,0.75, na.rm = TRUE) + 1.5*IQR(value, na.rm = TRUE)
  idx1 <- which(value - thres <= 0)
  val <- value[idx1]
  idx2 <- which.min(abs(val - thres))
  new_val <- val[idx2]
  return(new_val)
}


## Correlations
cis_mrna_protein <- read.delim("data/mrna_protein/ALL_corresponding_RNA_protein_correlation.txt", sep = '\t')


## CORUM complexes
protein_complexes <- read.delim('data/mrna_protein/allComplexes.txt')
protein_complexes <- protein_complexes[protein_complexes$Organism == 'Human', ]

protein_complexes_list <- strsplit(protein_complexes$subunits.Gene.name., ';')
names(protein_complexes_list) <- protein_complexes$ComplexName

protein_complexes_list <- reshape2::melt(protein_complexes_list)
protein_complexes_list <- protein_complexes_list[protein_complexes_list$value %in% cis_mrna_protein$gene, ]
toremove <- unique(protein_complexes_list$L1[grep(paste(c('Spliceosome', 'spliceosome', 'Nop56p-associated pre-rRNA complex', 'Ribosome', 'ribosom', 'proteasome'), collapse = '|'), protein_complexes_list$L1)])

toremove <- protein_complexes_list$value[protein_complexes_list$L1 %in% toremove]
protein_complexes_list <- protein_complexes_list[!protein_complexes_list$value %in% toremove, ]
protein_complexes <- unique(protein_complexes_list$value)

cis_mrna_protein$Zcor_Spearman <-  FisherZ(cis_mrna_protein$spearman_corr)
cis_mrna_protein <- cis_mrna_protein[cis_mrna_protein$gene %in% protein_complexes, ]
############################################################################################################################
## Subcell data
subcell <- read.delim("data/mrna_protein/CellAtlasAnnotationsVsSubCellBarCode.txt")

combined_data <- data.frame('genes' <- cis_mrna_protein$gene, 
                            'cor_Spearman' = cis_mrna_protein$spearman_corr,
                            'Zcor_Spearman' = cis_mrna_protein$Zcor_Spearman,
                            'localityHPA' = subcell$HPAbin[match(cis_mrna_protein$gene, subcell$Gene)],
                            'localitySC' = subcell$SCneigh[match(cis_mrna_protein$gene, subcell$Gene)],
                            'multilocalityHPA' = subcell$Type[match(cis_mrna_protein$gene, subcell$Gene)],
                            'multilocalitySC' = subcell$SCtype[match(cis_mrna_protein$gene, subcell$Gene)])


per_col_sum <- as.data.frame(combined_data[!is.na(combined_data$localitySC), ] %>% 
                                       dplyr::group_by(localitySC) %>%
                                       dplyr::summarise(n = n(), m = median(cor_Spearman, na.rm = TRUE),
                                                        q3 = boxplot.top(cor_Spearman)))


combined_data$localitySC <- factor(combined_data$localitySC, levels = per_col_sum$localitySC[order(per_col_sum$m)])

pvalAnova <- anova(lm(Zcor_Spearman ~ localitySC, combined_data))$`Pr(>F)`[1]
pvalAnova <- ifelse(pvalAnova < 2.2 * 10^-16, 2.2*10^-16, pvalAnova)

tab_order <- 1:length(levels(combined_data$localitySC))
names(tab_order) <- levels(combined_data$localitySC)
ttest_groups <- as.character(levels(combined_data$localitySC))
ttest_groups <- as.data.frame(combn(ttest_groups, 2))
ttest_groups <- ttest_groups[, order(apply(ttest_groups,2, function(i) diff(tab_order[i])))]


ttest_res <- do.call(rbind, apply(ttest_groups,2, function(l) {
  pval <- t.test(combined_data$Zcor_Spearman[combined_data$localitySC == as.character(l[1])], 
                 combined_data$Zcor_Spearman[combined_data$localitySC == as.character(l[2])])
  pval <- ifelse(pval$p.value < 2.2 * 10^-16, 2.2*10^-16, pval$p.value)
  data.frame('Comparison1' = l[1],'Comparison2' = l[2], 'Pvalue' = pval)
  
}))

max_comp <- max(table(ttest_res$Comparison1))
y_val <- 1:max_comp * 2
ttest_res$y_val <- unlist(sapply(1:max_comp, function(k) y_val[k] + seq(1:(max_comp + 1-k))/2))

overall_cor <- round(median(cis_mrna_protein$spearman_corr, na.rm = TRUE),2)

# Sort by correlation
p <- ggplot(combined_data[!is.na(combined_data$localitySC), ], aes(x = localitySC, y = cor_Spearman, fill = localitySC)) + 
  geom_boxplot(notch = TRUE) + 
  geom_hline(yintercept = overall_cor, lty =2) + 
  annotate(geom = 'text', x = -Inf, y = Inf, label = paste0('Anova p-value = ', format(pvalAnova, digits = 2, scientific = TRUE)), hjust = -0.2, vjust = 1.2) + 
  geom_text(data = per_col_sum, aes(x = localitySC, y = q3, label = paste0("n=",n)), size = 3, 
            position=position_dodge(width=0.9), vjust=-1, hjust = 0.5) +
  geom_text(data = per_col_sum, aes(x = localitySC, y = q3, label = paste0("mean=",round(m,2))), size = 3, 
            position=position_dodge(width=0.9),  vjust=-3, hjust = 0.5) + 
  geom_signif(data = ttest_res, mapping = aes(xmin=Comparison1, xmax=Comparison2, y_position = y_val, annotations= paste0('P = ', format(Pvalue, digits = 2, scientific = TRUE))), tip_length=0.01, manual=TRUE,  inherit.aes = FALSE) +
  scale_fill_manual(name = 'Neighbourhood', values = c("4_Mitochondria" = "#E0B891",
                                                       "3_Cytosol" = "#A6DDEB",
                                                       "2_Nucleus" = "#D1D1D1",
                                                       "1_Secretory" = "#FCC22C")) +
  scale_x_discrete(name = "", labels =  sapply(strsplit(levels(combined_data$localitySC), '_'), function(i) i[2])) +
  scale_y_continuous(name = "mRNA-Protein correlation (Spearman)", limits = c(-0.5, 8), breaks = seq(-1, 1, 0.5)) +
  theme_bw() + 
  ggtitle(paste0("Subcellular localization")) + 
  guides(fill = FALSE) +
  theme(plot.title = element_text(hjust=0.5),
        panel.grid.major.x  = element_blank(),
        axis.text.x = element_text(hjust = 1, angle = 45),
        legend.position = 'none')

# saveRDS(p, 'Figures/RNA_protein_correlation_SC')
pdf("Figures/mRNA_Protein_correlation_SC_thres.pdf", width = 5, height = 7, paper = "special")
plot(p)
dev.off()
#################################################################################################
### mirTarbase
mirTarbase <- read.xlsx('data/mrna_protein/hsa_MTI.xlsx')
mirTarbase <- mirTarbase[!grepl(paste(c('Non-Functional MTI', 'Weak'), collapse = '|'), mirTarbase$Support.Type), ]
mirTarbase <- mirTarbase[!duplicated(paste0(mirTarbase$miRTarBase.ID, '_', mirTarbase$Target.Gene)), ]
cis_mrna_protein$mirTarbase <- ifelse(cis_mrna_protein$gene %in% mirTarbase$Target.Gene, 'Yes', 'No')

per_col_sum <- as.data.frame(cis_mrna_protein %>%
                               dplyr::group_by(mirTarbase) %>%
                               dplyr::summarise(n = n(), m = median(spearman_corr, na.rm = TRUE),
                                                q3 = boxplot.top(spearman_corr)))


w <- t.test(cis_mrna_protein$Zcor_Spearman[cis_mrna_protein$mirTarbase == "Yes"],
            cis_mrna_protein$Zcor_Spearman[cis_mrna_protein$mirTarbase  == "No"])

pval <- ifelse(w$p.value < 2.2 * 10^-16, 2.2*10^-16, w$p.value)


cis_mrna_protein$mirTarbase <- factor(cis_mrna_protein$mirTarbase,
                                      levels = c('No', 'Yes'), ordered = TRUE)

p <- ggplot(data = cis_mrna_protein, aes(x = mirTarbase, y = spearman_corr, fill = mirTarbase)) +
  geom_boxplot(notch = FALSE,show.legend = FALSE,
               position = position_dodge(width=0.9)) +
  geom_hline(yintercept = overall_cor, lty = 2) + 
  scale_fill_manual(values = c('grey','black')) +
  scale_x_discrete(name = "", breaks = c('No', 'Yes')) +
  scale_y_continuous(name = "mRNA-Protein correlation (Spearman)", limits = c(-0.5, 1.5), breaks = seq(-1, 1, 0.5)) +
  geom_text(data = per_col_sum, aes(x = mirTarbase, y = 1, label = paste0("n=",n)), size = 3,
            position=position_dodge(width=0.9), vjust= 0, hjust = 0.5) +
  geom_text(data = per_col_sum, aes(x = mirTarbase, y = 1, label = paste0("median=",round(m,2))), size = 3,
            position=position_dodge(width=0.9),  vjust=-2, hjust = 0.5) +
  geom_signif(y_position= 1.5, xmin=c(1), xmax=c(2),
              annotation=  paste0('P = ', format(pval,digits = 2, scientific = TRUE)), tip_length = 0.01, textsize = 3) +
  theme_bw() +
  ggtitle(paste0("miRNA targets")) + 
  theme(plot.title = element_text(hjust = 0.5), panel.grid.major.x  = element_blank())

# saveRDS(p,"Figures/mRNA_protein_miRNA_mirTarbase")

pdf("Figures/mRNA_protein_miRNA_mirTarbase.pdf", width = 3, height = 6, paper = "special")
plot(p)
dev.off()

# mirTarbaseSites <- as.data.frame(mirTarbase %>% group_by(Target.Gene) %>%
#                                    summarise(n = n()))
# 
# cis_mrna_protein$mirTarbase_freq <- mirTarbaseSites[match(cis_mrna_protein$gene, mirTarbaseSites$Target.Gene), 2]
# 
# p <- ggplot(cis_mrna_protein, aes(x = mirTarbase_freq, y  = spearman_corr)) +
#   geom_hex(bins = 70) +
#   scale_fill_continuous(name = 'Count', type = "viridis") +
#   geom_smooth(col = '#FDE725FF', method = 'gam') +
#   theme_bw() +
#   labs(x = 'Number of miRNAs', y = 'RNA-Protein correlation (Spearman)', 
#        title = paste0('mirTarbase')) + 
#   theme(plot.title = element_text(hjust = 0.5))
# 
# 
# pdf("Figures/mRNA_protein_miRNA_mirTarbase_freq.pdf", width = 6, height = 5, paper = "special")
# plot(p)
# dev.off()
##########################################################################################
## Degradation profile
degradation_profile <- read.xlsx('data/mRNA_protein/NED_1-s2.0-S009286741631248X-mmc4.xlsx')
degradation_profile <- degradation_profile[degradation_profile$Degradation.profile %in% c('NED', 'ED'), ]

cis_mrna_protein$degradation <- degradation_profile$Degradation.profile[match(cis_mrna_protein$gene, degradation_profile$Gene.names)]

per_col_sum <- as.data.frame(cis_mrna_protein[!is.na(cis_mrna_protein$degradation),] %>%
                               dplyr::group_by(degradation) %>%
                               dplyr::summarise(n = n(), m = median(spearman_corr, na.rm = TRUE),
                                                q3 = boxplot.top(spearman_corr)))


w <- t.test(cis_mrna_protein$Zcor_Spearman[cis_mrna_protein$degradation == "NED"],
            cis_mrna_protein$Zcor_Spearman[cis_mrna_protein$degradation  == "ED"])

pval <- ifelse(w$p.value < 2.2 * 10^-16, 2.2*10^-16, w$p.value)


cis_mrna_protein$degradation <- factor(cis_mrna_protein$degradation,
                                       levels = c('ED', 'NED'), ordered = TRUE)

p <- ggplot(data = cis_mrna_protein[!is.na(cis_mrna_protein$degradation), ], aes(x = degradation, y = spearman_corr, fill = degradation)) +
  # stat_boxplot(geom ='errorbar', linetype = 2, position = position_dodge(width=0.9)) +
  geom_boxplot(notch = FALSE,show.legend = FALSE,
               position = position_dodge(width=0.9)) +
  # geom_jitter(pch=21, colour = "white", alpha = 0.3, position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.2),
  #             size = 1) +
  geom_hline(yintercept = overall_cor, lty = 2) + 
  scale_fill_manual(values = c('white','red')) +
  scale_x_discrete(name = "", breaks = c('ED', 'NED')) +
  scale_y_continuous(name = "mRNA-Protein correlation (Spearman)", limits = c(-0.5, 1.5), breaks = seq(-1, 1, 0.5)) +
  geom_text(data = per_col_sum, aes(x = degradation, y = 1, label = paste0("n=",n)), size = 3,
            position=position_dodge(width=0.9), vjust= 0, hjust = 0.5) +
  geom_text(data = per_col_sum, aes(x = degradation, y = 1, label = paste0("median=",round(m,2))), size = 3,
            position=position_dodge(width=0.9),  vjust=-2, hjust = 0.5) +
  geom_signif(y_position= 1.5, xmin=c(1), xmax=c(2),
              annotation=  paste0('P = ', format(pval,digits = 2, scientific = TRUE)), tip_length = 0.01, textsize = 3) +
  theme_bw() +
  ggtitle("Degradation profile") +
  theme(plot.title = element_text(hjust=0.5),
        panel.grid.major.x  = element_blank())

# saveRDS(p, "Figures/mRNA_protein_degradation_profile")

pdf("Figures/mRNA_protein_degradation_profile.pdf" , width = 3.5, height = 6, paper = "special")
plot(p)
dev.off()
##########################################################################################
## Combined miRNA - degradation profile
cis_mrna_protein_miRNA_deg <- cis_mrna_protein[!is.na(cis_mrna_protein$degradation), ]
cis_mrna_protein_miRNA_deg$miRNA_deg <- as.factor(paste(cis_mrna_protein_miRNA_deg$mirTarbase, cis_mrna_protein_miRNA_deg$degradation, sep = '_'))

ttest_groups <- levels(cis_mrna_protein_miRNA_deg$miRNA_deg)
ttest_groups <- as.data.frame(combn(ttest_groups, 2))
ttest_groups <- ttest_groups[, c(1, 6)]

ttest_res <- do.call(rbind, apply(ttest_groups,2, function(l) {
  pval <- t.test(cis_mrna_protein_miRNA_deg$Zcor_Spearman[cis_mrna_protein_miRNA_deg$miRNA_deg == as.character(l[1])], 
                 cis_mrna_protein_miRNA_deg$Zcor_Spearman[cis_mrna_protein_miRNA_deg$miRNA_deg == as.character(l[2])])
  pval <- ifelse(pval$p.value < 2.2 * 10^-16, 2.2*10^-16, pval$p.value)
  data.frame('Comparison1' = l[1],'Comparison2' = l[2], 'Pvalue' = pval)
  
}))

ttest_res$y_val <- 1.3
ttest_res$mirTarbase <- c('No','Yes')


per_col_sum <- as.data.frame(cis_mrna_protein_miRNA_deg %>%
                               dplyr::group_by(miRNA_deg) %>%
                               dplyr::summarise(n = n(), m = median(spearman_corr, na.rm = TRUE),
                                                q3 = boxplot.top(spearman_corr)))

per_col_sum$mirTarbase <- c(rep('No', 2), rep('Yes', 2))


overall_cor <- round(median(cis_mrna_protein$spearman_corr, na.rm = TRUE),2)

facet_labs <- c("miRNA target - No", "miRNA target - Yes")
names(facet_labs) <- c('No', 'Yes')

label_text <- data.frame(text =c("miRNA target - No", "miRNA target - Yes"), 
                         mirTarbase = c('No', 'Yes'))

p <- ggplot(cis_mrna_protein_miRNA_deg, aes(x = miRNA_deg, y = spearman_corr, fill = miRNA_deg)) + 
  geom_boxplot() + 
  geom_hline(yintercept = overall_cor, lty =2) + 
  geom_text(data = per_col_sum, aes(x = miRNA_deg, y = q3, label = paste0("n=",n)), size = 3, 
            position=position_dodge(width=0.9), vjust=-1, hjust = 0.5) +
  geom_text(data = per_col_sum, aes(x = miRNA_deg, y = q3, label = paste0("mean=",round(m,2))), size = 3, position=position_dodge(width=0.9),  vjust=-3, hjust = 0.5) + 
  geom_label(label_text, x = 1.5, y = 1.5, mapping = aes(label = text), inherit.aes = FALSE, fill = c('grey', 'black'), color= 'white') + 
  geom_signif(data = ttest_res,
              mapping = aes(xmin=Comparison1, xmax=Comparison2, y_position = y_val, 
                            annotations = paste0('P = ', format(Pvalue, digits = 2, scientific = TRUE))),
              tip_length=0.01, manual=TRUE, inherit.aes = FALSE, textsize = 3)  +
  scale_fill_manual(values = c('Yes_ED' = 'white',  'No_ED' = 'white', 'Yes_NED' = 'red',  'No_NED' = 'red')) + 
  facet_grid(.~mirTarbase, scales = 'free', labeller = labeller(mirTarbase = facet_labs)) + 
  scale_x_discrete(name = "", labels = c('ED', 'NED', 'ED', 'NED')) +
  scale_y_continuous(name = "mRNA-Protein correlation (Spearman)", limits = c(-0.5, 1.5), breaks = seq(-1, 1, 0.5)) +
  theme_bw() + 
  ggtitle(paste0("miRNA targets / Degradation profile")) + 
  theme(plot.title = element_text(hjust=0.5),
        panel.grid.major.x  = element_blank(),panel.spacing = unit(0,'cm'),
        legend.position = 'none', 
        strip.background = element_blank(), strip.text = element_blank())


# saveRDS(p, 'Figures/mRNA_protein_degradation_profile_miRNA')

pdf("Figures/mRNA_protein_degradation_profile_miRNA.pdf", width = 6, height = 6, paper = "special")
plot(p)
dev.off()

###