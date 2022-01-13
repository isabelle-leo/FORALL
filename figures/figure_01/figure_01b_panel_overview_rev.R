library(ComplexHeatmap)

# Set outputfolder
outputFolder <- file.path("figures", "output", "figure_01")

# create major output folder, if it doesn't exist
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder)

df <- as.data.frame(read.table(file = file.path("meta", "overview_HM_v2_T-ALL.txt"), sep ="\t",header=F))
df_anno <- as.data.frame(read.table(file = file.path("meta", "overview_annotation_T-ALL.txt"), sep ="\t",header=F))

colnames(df)<- df[1,]
df<-df[-1,]
rownames(df)<-df[,1]
df<-df[,-1]

colnames(df_anno)<- df_anno[1,]
df_anno<-df_anno[-1,]
rownames(df_anno)<-df_anno[,1]
df_anno<-df_anno[,-1]

#colors
categorical_col = c("TRUE" = "#6e96c9", "FALSE" = "#ffffff")
categorical_col_ms = c("TRUE" = "#742c55", "FALSE" = "#e3e6e8")
categorical_col_rna = c("TRUE" = "#28a1bc", "FALSE" = "#e3e6e8")
categorical_col_dsrt = c("TRUE" = "#66c1bf", "FALSE" = "#e3e6e8")

type_col = c("T-ALL" = "#f6b2af", "B-ALL" = "#b2cde2", "BCP-ALL" = "#cce3c3", "EBV" = "#949994")
tissue_col = c("BM" = "#676767", "PB" = "#b0b0b0", "PE" = "#eaeaea")
subtype_col = c("ABL1-ZMIZ1"= "#66b66e", "B-Other"="#f094a6", "BCL11B-TLX3"="#7c817d","BCR-ABL1"="#ea5c55","EBV"="#ea5c55","ETV6-RUNX1"="#2370b7","IGH-MYC"="#0ca9d4","KMT2A-AFF1"="#4bb796","KMT2A-MLLT1"="#26ac79","MEF2D-HNRNPUL1"="#fbd773","NUP214-ABL1"="#f9e683","PAX5-ETV6"="#f4eb73","T-Other"="#b7b5db","TCF3-PBX1"="#c5a0ca")
#annotation
col_anno<- HeatmapAnnotation(df = df_anno,
                             col = list("RNAseq" = categorical_col_rna, "LC-MS/MS-Replicate" = categorical_col_ms, "DSRT" = categorical_col_dsrt,"RNAseq-replicate" = categorical_col_rna , "Type" = type_col, "Tissue" = tissue_col,"Subtype" = subtype_col),
                             gp = gpar(col = "White"),
                             annotation_name_side = "left",
                             annotation_legend_param = list(
                                     "Subtype" = list(
                                             title = "Subtype",
                                             nrow = 2
                                     ),
                                     "Tissue" = list(
                                             nrow = 2
                                     ),
                                     "Type" = list(
                                             nrow = 2
                                     )),
                             show_legend = c(T,T,F,F,F,F,T),
                             which = 'col')
#print file
pdf(file = file.path(outputFolder, "panel_overview_Rev.pdf"), width = 12, height = 8)
draw(Heatmap(as.matrix(df),col=categorical_col,  cluster_rows = F, cluster_columns = F,show_row_dend = F,show_column_dend = F,row_names_side = "left",column_names_side = "top",
             rect_gp = gpar(col = "#e3e6e8", lwd = 0.2),
             row_title = "Gene Fusion", column_title_gp = gpar(fontsize = 11),
             column_names_gp = gpar(fontsize = 9),
             show_heatmap_legend = F,
             top_annotation = col_anno),
     annotation_legend_side = "bottom")
dev.off()
