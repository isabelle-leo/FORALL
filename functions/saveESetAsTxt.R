#' Save an ExpressionSet as plain but annotated txt
#'
#' @param eSet ExpressionSet
#' @param file destination
#' @param colNameColumn identifier for final column names
#' @param rowNameColumn identifier for final row names
#' @param additionalRows names of additional row annotations
#' @param additionalColumns names of additional column annotations
#' @param extendedColnames use extended column names?
#'
#' @return only file output
#' @importFrom magrittr "%>%"
#' @importFrom Biobase exprs pData fData
#' @export
#'
#' @examples
saveESetAsTxt <- function (eSet,
                           file,
                           colNameColumn = "proteomics_id",
                           rowNameColumn = "gene_symbol",
                           additionalRows = NULL,
                           additionalColumns = NULL,
                           extendedColnames = TRUE) {

  directory <- gsub(pattern = paste0(basename(file), "$"),
                    replacement = "",
                    x = file)
  if (!dir.exists(directory)) dir.create(directory, recursive = TRUE)

  num.add.rows <- 0

  df.save <- eSet %>%
    exprs() %>%
    as.data.frame()

  if (extendedColnames) {
    colnames(df.save) <- paste(pData(eSet)[[colNameColumn]], colnames(df.save), sep = "_")
  } else {
    colnames(df.save) <- pData(eSet)[[colNameColumn]]
  }

  if (!is.null(additionalRows)) {
    # check if columns exist
    if (!all(additionalRows %in% colnames(pData(eSet)))) stop("Column name not existing!")

    add.rows <- as.data.frame(pData(eSet)[, additionalRows])
    num.add.rows <- ncol(add.rows)
    add.rows.t <- t(add.rows)
    row.names(add.rows.t) <- additionalRows
    colnames(add.rows.t) <- colnames(df.save)
    df.save <- rbind(add.rows.t, as.matrix(df.save))
  }

  if (!is.null(additionalColumns)) {
    # check if columns exist
    if (!all(additionalColumns %in% colnames(fData(eSet)))) stop("Column name not existing!")

    add.cols <- as.data.frame(fData(eSet)[, additionalColumns])
    colnames(add.cols) <- additionalColumns
    dummy.matrix <- matrix(nrow = num.add.rows, ncol = ncol(add.cols), data = NA)
    colnames(dummy.matrix) <- additionalColumns
    if (num.add.rows > 0) row.names(dummy.matrix) <- row.names(df.save[1:num.add.rows, ])
    add.cols.full <- rbind(dummy.matrix, add.cols)
    df.save <- cbind(add.cols.full, df.save)
  }

  row.names(df.save)[(num.add.rows+1):length(row.names(df.save))] <- fData(eSet)[[rowNameColumn]]

  write.table(x = df.save,
              file = file,
              quote = FALSE,
              sep = "\t",
              row.names = TRUE,
              col.names = NA)
}
