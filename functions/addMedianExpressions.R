#' Calculate median expression per category
#'
#' @param eSet ExpressionSet (with CoSC annotation)
#' @param use which pdata column should be used for cluster retrieval
#' @param removeReplicates should replicates be removed (by isReplicate)
#' @param recognitionString string for quantitative column names in fData
#'
#' @return updated ExpressionSet
#' @importFrom Biobase pData exprs fData
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom matrixStats rowMedians
#' @export
#'
#' @examples
addMedianExpressions <- function (eSet,
                                  use,
                                  removeReplicates = TRUE,
                                  recognitionString = "$$cosc_cluster.medianExpression.") {
  
  # get all possible clusters
  clusters <- eSet %>%
    pData() %>%
    dplyr::select(!!use) %>%
    unlist() %>%
    unique() %>%
    as.character()
  
  # calculate all median expressions
  medianExpressions <- sapply(X = clusters,
                              FUN = function (cluster) {
                                # get expression values of cluster
                                if (removeReplicates) {
                                  sub <- eSet[, pData(eSet)[[use]] == cluster & !pData(eSet)$isReplicate] %>%
                                    exprs()
                                } else {
                                  sub <- eSet[, pData(eSet)[[use]] == cluster] %>%
                                    exprs()
                                }
                                
                                # calculate row medians
                                rowMedians(sub, na.rm = TRUE)
                              },
                              USE.NAMES = TRUE)
  colnames(medianExpressions) <- paste0(recognitionString, colnames(medianExpressions))
  
  # append them to fData
  eSet.updated <- eSet
  fData(eSet.updated) <- cbind(fData(eSet.updated), medianExpressions)
  
  return(eSet.updated)
}
