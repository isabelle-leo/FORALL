#' Annotate an ExpressionSet with meta data
#'
#' @param expressions expression set
#' @param subcellLocFile file path for subcellular localization
#' @param geneAnnotationFolder folder with standardised annotations
#'
#' @return annotated ExpressionSet
annotateExpressionSet <- function (expressions,
                                   subcellLocFile,
                                   geneAnnotationFolder) {
  
  # -- add subcellular localization
  subcellLoc <- getData(subcellLocFile) %>%
    rename(gene_symbol = SYMBOL) %>%
    rename(localization = Neighbourhood) %>%
    select(gene_symbol, localization)
  
  # fuse subcell loc to features
  fData(expressions) <- fData(expressions) %>%
    left_join(subcellLoc, by = "gene_symbol") %>%
    mutate(localization = ifelse(is.na(localization), "Unknown", localization))
  
  
  # -- add other standardised annotations
  files <- list.files(path = geneAnnotationFolder, full.names = TRUE)
  getGenes <- function (file) {
    read_tsv(file = file,
             comment = "#",
             col_names = FALSE) %>%
      unlist()
  }
  
  # prepare gene lists
  genes <- lapply(X = files,
                  FUN = getGenes)
  
  names(genes) <- lapply(X = genes,
                         FUN = function (x) x[[1]])
  
  genes <- lapply(X = genes,
                  FUN = function (x) x[2:length(x)])
  
  # get the current gene list
  featureGenes <- expressions %>%
    fData() %>%
    .$gene_symbol
  
  # check for each annotation, if the gene is in
  for (category in names(genes)) {
    fData(expressions)[[paste0("is", category)]] <- featureGenes %in% genes[[category]]
  }
  
  return(expressions)
}
