#' Load large data files
#'
#' @param fileName the path and file name
#' @param sep the separator
#' @param quote the quote character(s)
#'
#' @return a dataframe
#' @importFrom data.table fread
#' @export
#'
#' @examples
getData <- function (fileName, sep = "\t", quote = "") {
  return(fread(input = fileName,
               sep = sep,
               quote = quote,
               header = TRUE,
               stringsAsFactors = FALSE,
               data.table = FALSE))
}
