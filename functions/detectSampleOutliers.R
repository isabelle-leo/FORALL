#' Detect sample outliers
#'
#' @param eSet ExpressionSet
#' @param pearsonCutOff cutoff for minimal Pearson correlation between samples
#' @param outputFolder output folder
#'
#' @return ExpressionSet with annotation which samples are outliers
detectSampleOutliers <- function (eSet,
                                  pearsonCutOff = .5,
                                  outputFolder = NULL) {

  if (!dir.exists(outputFolder)) dir.create(outputFolder, recursive = TRUE)

  # no replicates, no NAs
  eSet.filt <- eSet[, !pData(eSet)$isReplicate] %>% removeNAsFromESet()

  # calculate the pearson correlation matrix among samples
  cor.mat <- cor(exprs(eSet.filt), method = "pearson")
  diag(cor.mat) <- 0

  # get maximum inter-sample correlation
  cor.max <- apply(cor.mat, 1, max)

  # sort them ascendingly
  cor.max.sort <- sort(cor.max)

  # get outliers and non-outliers by pearson cutoff
  outliers <- cor.max.sort[cor.max.sort < pearsonCutOff]
  non.outliers <- cor.max.sort[cor.max.sort >= pearsonCutOff]

  if (!is.null(outputFolder)) pdf(file = file.path(outputFolder, paste0("sample_outliers_at_p_", pearsonCutOff, ".pdf")))
  plot(cor.max.sort,
       col = "green",
       pch = 16,
       cex = 1.5,
       main = "Maximum inter-sample correlation")
  points(outliers,
         col = "red",
         pch = 16,
         cex = 1.5)
  abline(h = pearsonCutOff)
  if (!is.null(outputFolder)) dev.off()

  # add info to pData
  pData(eSet) <- pData(eSet) %>%
    mutate(isOutlier = ifelse(position %in% names(outliers), TRUE, FALSE))

  return(eSet)
}
