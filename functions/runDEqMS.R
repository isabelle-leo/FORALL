#' Run a DEqMS analysis
#'
#' @param eSet ExpressionSet with log expression values
#' @param conditions vector with different sample conditions
#' @param contrasts the comparisons to be done
#' @param validValuesPostFilter if set, minimum number of valid values per comparison group
#' @param outputFolder an output location (not required)
#'
#' @return a result object
#' @importFrom magrittr %>%
#' @importFrom Biobase exprs fData
#' @importFrom matrixStats rowMins
#' @importFrom limma makeContrasts lmFit contrasts.fit eBayes
#' @importFrom tibble column_to_rownames
#' @importFrom DEqMS spectraCounteBayes plotFitCurve outputResult
#' @export
#'
#' @examples
runDEqMS <- function (eSet,
                      conditions,
                      contrasts,
                      validValuesPostFilter = 3,
                      outputFolder = NULL) {
  
  if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder, recursive = TRUE)
  
  # make a dataframe with gene symbols as rownames
  data <- eSet %>%
    exprs() %>%
    as.data.frame()
  
  row.names(data) <- fData(eSet)$gene_symbol
  
  # make the design
  design <- model.matrix(~ 0 + conditions)
  colnames(design) <- gsub("conditions", "", colnames(design))
  
  # merge design and contrasts
  contrasts.design <- makeContrasts(contrasts = contrasts,
                                    levels = design)
  
  # get all features with at least one valid value (group-wise)
  validFeatures <- lapply(X = unique(conditions),
                          FUN = function (condition) {
                            # get data subset of this sample group
                            data.subset <- data %>%
                              dplyr::select(colnames(.)[condition == conditions])
                            
                            # count NAs
                            data.subset$valid <- apply(X = data.subset,
                                                       MARGIN = 1,
                                                       FUN = function (x) sum(!is.na(x)))
                            
                            # filter for features with at least 1 valid value
                            data.subset.filtered <- data.subset %>%
                              filter(valid > 0)
                            
                            # return the feature names
                            return(row.names(data.subset.filtered))
                          })
  
  # add a minimum PSM count to the filter
  # only PSM counts from sets that are included in the analysis should be used
  # extract sets from position
  sets <- eSet %>%
    pData %>%
    .$position %>%
    gsub("_tmt10plex_.*$", "", .) %>%
    unique()
  
  # define PSM count columns
  psmCountCols <- paste(sets, "_quanted_psm_count", sep = "")
  
  # get min quanted PSM count
  minPsms <- eSet %>%
    fData() %>%
    dplyr::select(all_of(psmCountCols)) %>%
    dplyr::mutate(across(everything(), ~ ifelse(.x == 0, NA, .))) %>%
    dplyr::mutate(min_psms = rowMins(as.matrix(.), na.rm = TRUE)) %>%
    `row.names<-`(fData(eSet)$gene_symbol) %>%
    filter(is.finite(min_psms))
  
  validFeatures$minPsms <- row.names(minPsms)
  
  # get features fulfilling criterion in every sample group
  validFeatures.overlap <- Reduce(intersect, validFeatures)
  
  # only use data fulfilling filter criteria
  data.filtered <- data[validFeatures.overlap, ]
  
  psms.filtered <- minPsms[validFeatures.overlap, ]$min_psms
  
  # linear fit
  fit1 <- lmFit(object = as.matrix(data.filtered),
                design = design)
  
  fit2 <- eBayes(contrasts.fit(fit = fit1,
                               contrasts = contrasts.design))
  
  # add PSM counts
  fit2$count <- psms.filtered
  
  fit3 <- spectraCounteBayes(fit = fit2)
  
  # print QC stuff
  if (!is.null(outputFolder)) pdf(file = file.path(outputFolder, "DEqMS_variance_boxplot.pdf"))
  VarianceBoxplot(fit = fit3)
  if (!is.null(outputFolder)) dev.off()
  
  # output tables
  tables <- c()
  for (col in 1:length(contrasts)) {
    selectedContrast <- contrasts[[col]]
    
    # get the number of valid values per contrast group
    contrastDesign <- selectedContrast %>%
      str_split("-") %>%
      unlist() %>%
      design[, .]
    
    validValuesPerContrast <- sapply(X = contrastDesign %>% colnames(),
                                     FUN = function (contrast) {
                                       selector <- contrastDesign[, contrast] != 0
                                       eSet[, selector] %>%
                                         esApply(MARGIN = 1, function (row) sum(!is.na(row))) %>%
                                         data.frame(names(.), .) %>%
                                         set_names(c("gene", paste0("valid_values_", contrast)))
                                     },
                                     USE.NAMES = TRUE,
                                     simplify = FALSE) %>%
      Reduce(f = function (x, y) left_join(x, y, by = "gene"), x = .)
    
    outputTable <- outputResult(fit = fit3, coef_col = col) %>%
      left_join(validValuesPerContrast, by = "gene")
    
    if (is.numeric(validValuesPostFilter)) {
      outputTable <- outputTable %>%
        filter(across(starts_with("valid_values_"), ~ .x >= validValuesPostFilter))
    }
    
    tables[[selectedContrast]] <- outputTable
    if (!is.null(outputFolder)) {
      # add num valid values per contrast
      write.table(x = outputTable,
                  file = file.path(outputFolder, paste0(contrasts[[col]], "_DEqMS_table.txt")),
                  quote = FALSE,
                  sep = "\t",
                  row.names = TRUE,
                  col.names = NA)
    }
  }
  
  return(list(fit = fit3,
              outputTables = tables))
}
