clusterFeatures <- function (eSet,
                             method = "pearson",
                             linkage = "average",
                             k = 2) {
  
  # get matrix of full overlap protein quants
  eSet_filt <- eSet[, !pData(eSet)$isReplicate] %>%
    removeNAsFromESet()
  
  mat <- eSet_filt %>%
    exprs() %>%
    t()
  
  # hierarchical clustering of features with k groups
  clusters_features <- as.dist(1 - cor(mat, method = method)) %>%
    hclust(method = linkage)
  
  # write cluster number to fData
  df <- data.frame(gene_symbol = fData(eSet_filt)$gene_symbol,
                   order = clusters_features$order,
                   stringsAsFactors = FALSE)
  
  eSet_updated <- addSparseFData(eSet = eSet,
                                 df = df,
                                 by = "gene_symbol")
  
  return(eSet_updated)
}
