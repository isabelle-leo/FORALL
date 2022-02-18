#' Tidy up and transform an ExpressionSet
#'
#' @param expressions.raw ExpressionSet
#' @param idCol column name for pool identification
#' @param poolName identifier of pool samples
#' @param doLog2 log2 the values
#'
#' @return tidied ExpressionSet
#' @importFrom magrittr %>%
#' @importFrom Biobase pData exprs esApply fData fData<-
#' @importFrom dplyr filter
#' @export
#'
#' @examples
tidyExpressionSet <- function (expressions.raw,
                               idCol = "proteomics_id",
                               poolName = "Pool",
                               doLog2 = FALSE) {

  # remove pool from eSet
  samplesWoPool <- pData(expressions.raw)[which(pData(expressions.raw)[[idCol]] != poolName), ]
  expressions <- expressions.raw[, samplesWoPool$position]

  # log2 transformation
  if (doLog2) {
    exprs(expressions) <- log2(exprs(expressions))
  }

  # add number of NAs
  numNAsPerRow <- esApply(X = expressions,
                          MARGIN = 1,
                          FUN = function (row) {
                            sum(is.na(row))
                          })

  fData(expressions)$numNAs <- numNAsPerRow

  # check if removeNAs function is working
  if (!identical(fData(expressions)[which(fData(expressions)$numNAs == 0), ]$gene_symbol,
                 fData(removeNAsFromESet(expressions))$gene_symbol)) {
    stop("RemoveNA function not working properly!")
  }

  return(expressions)
}
