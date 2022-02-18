#' Detect high variance expressions from ExpressionSet.
#' (based on https://github.com/aleferna/BCLandscape)
#'
#' \code{detectHighVarianceExpressions} returns an ExpressionSet with expression
#'   values considered as high-variance expressions
#'
#' @param eSet ExpressionSet with quantification values
#' @param identifier feature identifier
#' @param maxit maximum iterations for the Gaussian micture model
#' @param outputFolder output folder
#' @param nCores number of compute cores to use
#'
#' @return an ExpressionSet which contains annotations for high-variance features
detectHighVarianceExpressions <- function (eSet,
                                           identifier = "gene_symbol",
                                           nCores = 2,
                                           maxit = 1000,
                                           outputFolder = file.path("output")) {
  
  if (!dir.exists(outputFolder)) dir.create(path = outputFolder, recursive = TRUE)
  
  set.seed(262831)
  
  sdQ <- function (x) {
    sd(sort(x)[2:(length(x) - 1)])
  }
  
  registerDoMC(nCores)
  
  # calculate the log2 quantile standard deviation of totally overlapping expressions
  eSet.noNAs <- eSet[, !pData(eSet)$isReplicate] %>%
    removeNAsFromESet()
  
  geneSymbols <- eSet.noNAs %>%
    fData() %>%
    .[[identifier]]
  
  expressions <-  eSet.noNAs %>% exprs()
  row.names(expressions) <- geneSymbols
  expressions.sdQ <- log2(apply(X = expressions,
                                MARGIN = 1,
                                FUN = sdQ))
  
  pdf(file = file.path(outputFolder, "density.pdf"))
  plot(density(expressions.sdQ),
       xlab = "log2 Standard deviation",
       ylab = "Density")
  dev.off()
  
  # Gaussian mixture model
  gauMix <- normalmixEM(x = expressions.sdQ,
                        maxit = maxit)
  
  # compare final means and determine foreground and background
  if (gauMix$mu[[1]] > gauMix$mu[[2]]){
    idxfg <- 1
    idxbg <- 2
  } else {
    idxfg <- 2
    idxbg <- 1
  }
  
  # make a df with important data
  gmix <- data.frame(fgMu = gauMix$mu[idxfg],
                     fgSigma = gauMix$sigma[idxfg],
                     fgLambda = gauMix$lambda[idxfg],
                     bgMu = gauMix$mu[idxbg],
                     bgSigma = gauMix$sigma[idxbg],
                     bgLambda = gauMix$lambda[idxbg])
  
  # make normal distributions
  fg <- rnorm(1E5, mean = gmix$fgMu, sd = gmix$fgSigma)
  bg <- rnorm(1E5, mean = gmix$bgMu, sd = gmix$bgSigma)
  
  # start PDF
  pdf(file = file.path(outputFolder, "variance_gaussian_mix.pdf"))
  
  d <- density(expressions.sdQ)
  plot(d$x,
       d$y / max(d$y),
       type = "l",
       lwd = 5,
       xlab = "log2 Standard deviation",
       ylab = "Density",
       main = "Mix model & simulated data")
  
  d.fg <- density(fg)
  lines(d.fg$x,
        gmix$fgLambda * d.fg$y / max(d.fg$y),
        col = "blue",
        lwd = 2)
  
  d.bg <- density(bg)
  lines(d.bg$x,
        gmix$bgLambda * d.bg$y / max(d.bg$y),
        col = "red",
        lwd = 2)
  
  # simulate the data with lambda distribution
  fg.lambda <- rnorm(1E5 * gmix$fgLambda, mean = gmix$fgMu, sd = gmix$fgSigma)
  bg.lambda <- rnorm(1E5 * gmix$bgLambda, mean = gmix$bgMu, sd = gmix$bgSigma)
  
  # watch how well the simulation fits the original data
  d.comb <- density(c(fg.lambda, bg.lambda))
  lines(d.comb$x,
        d.comb$y/max(d.comb$y),
        col = "green",
        lwd = 2)
  
  # estimate ROC analysis
  # qt: quantile
  roc <- foreach(qt = 1:99 , .combine = rbind) %dopar% {
    # quantile loop
    xroc <- foreach(i = 1:10 , .combine = rbind) %do% {
      # 10 different samplings for normal distributions
      bg.iter <- rnorm(1E5 * gmix$bgLambda, mean =  gmix$bgMu, sd = gmix$bgSigma)
      fg.iter <- rnorm(1E5 * gmix$fgLambda, mean =  gmix$fgMu, sd = gmix$fgSigma)
      
      quantile.qt <- quantile(c(fg.iter, bg.iter), qt / 100.0)
      
      TP <- sum(fg.iter > quantile.qt)
      FP <- sum(bg.iter > quantile.qt)
      TN <- sum(bg.iter < quantile.qt)
      FN <- sum(fg.iter < quantile.qt)
      
      Sensitivity <- TP / (TP + FN)
      Specificty <- TN / (TN + FP)
      TPr <- TP / (TP + FP)
      score <- sum(fg.iter > quantile.qt) - sum(bg.iter > quantile.qt)
      
      data.frame(qt, quantile.qt, Sensitivity, Specificty, TP, FP, TN , FN, TPr, score)
    }
    colMeans(xroc)
  }
  
  # identify best quantile for separation
  roc <- as.data.frame(roc)
  idx <- roc$score == max(roc$score)
  quantile.qt.max <- roc$quantile.qt[idx]
  
  # round to the closest 0.5
  # quantile.qt.max <- ceiling(quantile.qt.max * 2) / 2
  
  # round to the lower .5
  quantile.qt.max <- floor(quantile.qt.max * 2) / 2
  
  # in which row is the cutoff?
  idx2 <- which.min(abs(roc$quantile.qt - quantile.qt.max))
  
  # plot the stuff
  plot(roc$quantile.qt,
       roc$score,
       type = "l",
       lwd = 2,
       xlab = "Variance Threshold",
       ylab = "Score (#TP - #FP)" ,
       main = paste("Threshold = ", quantile.qt.max))
  points(roc$quantile.qt[idx],
         roc$score[idx],
         col = "orange",
         cex = 2,
         xlim = c(0, 1),
         pch = 16)
  points(roc$quantile.qt[idx2],
         roc$score[idx2],
         col = "red",
         cex = 2,
         xlim = c(0, 1),
         pch = 16)
  
  # plot ROC curve
  plot(1 - roc$Specificty,
       roc$Sensitivity,
       type = "l",
       lwd = 3,
       ylab = "True Positive Rate",
       xlab = "False Positive Rate",
       main = 'ROC analysis')
  points(1 - roc$Specificty[idx],
         roc$Sensitivity[idx],
         col = "orange",
         pch = 16,
         cex = 2)
  points(1 - roc$Specificty[idx2],
         roc$Sensitivity[idx2],
         col = "red",
         pch = 16,
         cex = 2)
  
  expressionNames <- names(expressions.sdQ[expressions.sdQ > quantile.qt.max])
  text(0.5,
       0.5,
       labels = paste0("Selected features: ", length(expressionNames)))
  text(0.5,
       0.45,
       labels = paste0("Estimated true positives: ~",  as.integer(length(expressionNames) * roc$TPr[idx2])))
  text(0.5,
       0.4,
       labels = paste0("Estimated false positives: ~", as.integer(length(expressionNames) * (1 - roc$TPr[idx2]))))
  
  dev.off()
  
  # add feature data to eSet
  allGeneSymbols <- fData(eSet)[[identifier]]
  fData(eSet)$hasHighVariance <- allGeneSymbols %in% expressionNames

  return(eSet)
}
