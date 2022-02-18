#' Consensus clustering of samples
#'
#' @param eSet ExpressionSet
#' @param maxK maximum number of clusters
#' @param outputFolder output folder
#'
#' @return only plots in output folder
#' @importFrom ConsensusClusterPlus ConsensusClusterPlus
#' @export
#'
#' @examples
consensusClusterSamples <- function (eSet,
                                     maxK = 15,
                                     outputFolder) {
  
  if (!dir.exists(outputFolder)) dir.create(outputFolder, recursive = TRUE)
  
  methods <- c("pearson", "spearman")
  linkages <- c("average", "ward.D2", "complete")
  
  res <- list()
  for (method in methods) {
    method_res <- list()
    for (linkage in linkages) {
      d <- file.path(outputFolder, method, linkage)
      if (!dir.exists(d)) dir.create(d, recursive = TRUE)
      
      method_res[[linkage]] <- ConsensusClusterPlus(d = data.matrix(exprs(eSet)),
                                                    maxK = maxK,
                                                    reps = 2000,
                                                    pFeature = 0.8,
                                                    pItem = 0.8,
                                                    clusterAlg = "hc",
                                                    innerLinkage = linkage,
                                                    finalLinkage = linkage,
                                                    distance = method,
                                                    seed = 123,
                                                    plot = "pdf",
                                                    title = d)
    }
    res[[method]] <- method_res
  }
  return(res)
}
