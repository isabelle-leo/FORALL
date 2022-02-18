#' Cluster samples to core-sample clusters (hierarchical)
#'
#' @param eSet ExpressionSet
#' @param k number of clusters
#' @param method the correlation method (pearson, spearman, kendall)
#' @param linkage the lonkage method (default: average)
#' @param heatmapAnnotation the pData column for heat map annotation
#' @param outputFolder output folder
#'
#' @return updated ExpressionSet with cluster number in pData
clusterSamplesHierarchical <- function (eSet,
                                        k = 2,
                                        method = "pearson",
                                        linkage = "average",
                                        heatmapAnnotation = "cell_line",
                                        outputFolder = NULL) {
  
  o <- NULL
  if (!is.null(outputFolder)) o <- file.path(outputFolder, "hierarchical", method, linkage)
  if (!is.null(o) && !dir.exists(o)) dir.create(o, recursive = TRUE)
  
  # get matrix of full overlap protein quants and without replicates
  mat <- eSet[, !pData(eSet)$isReplicate] %>%
    removeNAsFromESet() %>%
    exprs()
  
  # hierarchical clustering of samples with k groups
  hClustersSamples <- as.dist(1 - cor(mat, method = method)) %>%
    hclust(method = linkage)
  
  # the cluster assignments
  clusters <- cutree(hClustersSamples, k)
  
  # update the eSet
  eSet.updated <- eSet
  identifier <- paste0("hierarCluster_", method, "_", linkage)
  dfPheno <- pData(eSet.updated) %>%
    left_join(data.frame(position = names(clusters),
                                   clusters = clusters,
                                   stringsAsFactors = FALSE),
              by = "position")
  
  # loop through NA clusters and find replicate
  naClusterSamples <- dfPheno %>%
    filter(is.na(clusters)) %>%
    .$Cell.Line.R
  
  for (sample in naClusterSamples) {
    # get the non-replicate sample
    pairSample <- gsub(".BR2", "", sample)
    
    # get it's cluster number
    clusterNumber <- dfPheno %>%
      filter(Cell.Line.R == pairSample) %>%
      .$cluster
    
    # assign this number to the replicate pair
    dfPheno[dfPheno$Cell.Line.R == sample, "clusters"] <- clusterNumber
    
    message("Cluster number expanded from ", pairSample, " to ", sample, ".")
  }
  
  # finally correct the column name for convenience
  colnames(dfPheno)[colnames(dfPheno) == "clusters"] <- identifier
  
  # and back to the expression set
  pData(eSet.updated) <- dfPheno
  
  # make a PCA for the samples
  matPca <- eSet.updated %>%
    removeNAsFromESet() %>%
    exprs()
  
  pca <- prcomp(x = t(matPca),
                scale. = TRUE)
  
  # convert the PCs into a dataframe
  pca$x <- as.data.frame(pca$x)
  
  # check if the PCA df kept the same sample order
  if (!identical(row.names(pca$x), pData(eSet.updated)$position)) {
    stop("Rownames inconsistency between PCA dataframe and pData.")
  }
  
  # add the cluster info to the PCA df
  pca$x$cluster <- factor(pData(eSet.updated)[[identifier]])
  
  # plot the PCA
  if (!is.null(o)) pdf(file = file.path(o, "sample_pca.pdf"), width = 12)
  print(ggplot(pca$x, aes(x = PC1, y = PC2, color = cluster)) +
          geom_point() +
          xlab(label = paste0("PC1 (", round(pca$sdev / sum(pca$sdev) * 100, 2)[[1]], "%)")) +
          ylab(label = paste0("PC2 (", round(pca$sdev / sum(pca$sdev) * 100, 2)[[2]], "%)")) +
          ggtitle(label = paste0("Sample PCA (k=", k, ")")) +
          guides(color = guide_legend(title = "CLC number")))
  if (!is.null(o)) dev.off()
  
  # plot a heat map
  ha.top <- HeatmapAnnotation(clusters = factor(pData(eSet.updated)[[identifier]]),
                              type = factor(pData(eSet.updated)$Type))
  
  ha.bottom <- HeatmapAnnotation(
    text = anno_text(x = pData(eSet.updated)[[heatmapAnnotation]],
                     rot = 90,
                     just = "right",
                     gp = gpar(fontsize = 7),
                     location = unit(2.5, "cm"),
                     height = unit(2.5, "cm")))
  
  hm <- Heatmap(
    name = "log2 expression ratio",
    matrix = exprs(removeNAsFromESet(eSet.updated)),
    col = colorRamp2(c(-1.5, 0, 1.5), c("#3869f2", "#ffffff", "#f23841")),
    show_row_names = nrow(mat) < 200,
    row_names_gp = gpar(fontsize = 5),
    show_column_names = FALSE,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_heatmap_legend = TRUE,
    clustering_distance_columns = method,
    clustering_distance_rows = "pearson",
    clustering_method_columns = linkage,
    clustering_method_rows = "average",
    column_split = factor(pData(eSet.updated)[[identifier]]),
    cluster_column_slices = TRUE,
    top_annotation = ha.top,
    bottom_annotation = ha.bottom
  )
  
  if (!is.null(o)) pdf(file = file.path(o, "heatmap.pdf"), width = 12)
  print(hm)
  if (!is.null(o)) dev.off()
  
  return(eSet.updated)
}
