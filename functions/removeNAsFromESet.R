#' Remove NA containing protein rows from ExpressionSet
#'
#' \code{removeNAsFromESet} returns an \code{ExpressionSet} with proteins
#'   containing a given ratio of NAs.
#'
#' @param eSet ExpressionSet with quantification values
#' @param na_ratio fraction of NAs allowed per feature
#'
#' @return an ExpressionSet that fulfills the NA ratio criterium
#' @importFrom Biobase esApply exprs
#' @export
#'
#' @examples
removeNAsFromESet <- function (eSet, na_ratio = 0) {

  # count NAs per row
  number_NAs_row <- esApply(X = eSet,
                            MARGIN = 1,
                            FUN = function (row) {
                              sum(is.na(row))
                            })

  # get threshold number
  threshold <- na_ratio * ncol(exprs(eSet))

  # get rows that fulfill threshold
  keep <- number_NAs_row <= threshold

  # subset eSet
  eSet_filtered <- eSet[keep, ]

  return(eSet_filtered)
}
