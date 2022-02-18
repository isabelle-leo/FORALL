#' Cluster samples to core-sample clusters (hierarchical)
#'
#' @param seur Seurat object sparse matrix
#' @param k number of clusters
#' @param method the correlation method (pearson, spearman, kendall)
#' @param linkage the linkage method (default: average)
#' @param heatmapAnnotation the pData column for heat map annotation
#' @param outputFolder output folder
#'
#' @return updated ExpressionSet with cluster number in pData
clusterSamplesHierarchicalSeurat <- function (seur,
                                        k = 2,
                                        method = "pearson",
                                        linkage = "average",
                                        heatmapAnnotation = "cell_line",
                                        outputFolder = NULL) {
  
  o <- NULL
  if (!is.null(outputFolder)) o <- file.path(outputFolder, "hierarchical", method, linkage)
  if (!is.null(o) && !dir.exists(o)) dir.create(o, recursive = TRUE)
 
   mat <- data.matrix(seur)
  
  # hierarchical clustering of samples with k groups
  hClustersSamples <- as.dist(1 - cor(mat, method = method)) %>%
    hclust(method = linkage)
  
  # the cluster assignments
  clusters <- cutree(hClustersSamples, k)
  
  return(clusters)
}
