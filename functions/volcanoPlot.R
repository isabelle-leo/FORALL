volcanoPlot <- function (data,
                         minXLabel = 2,
                         minYLabel = -log10(.05),
                         startColor = "rosybrown1",
                         stopColor = "firebrick",
                         outputFolder = NULL,
                         name = "volcano") {
  
  if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder, recursive = TRUE)
  
  # adapted from https://slowkow.com/notes/ggplot2-color-by-density/
  getDensity <- function(x, y, ...) {
    dens <- MASS::kde2d(x, y, ...)
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix, iy)
    return(dens$z[ii])
  }
  
  dat <- data.frame(x = data$logFC,
                    y = -log10(data$sca.adj.pval),
                    label = data$gene,
                    imputed = FALSE)
  
  # handle missing p-values
  dat$imputed <- is.infinite(dat$y)
  dat$y[dat$imputed] <- max(dat$y[!dat$imputed], na.rm = TRUE) * 1.1
  
  dat$density <- getDensity(x = dat$x,
                            y = dat$y,
                            h = c(2, 2),
                            n = 200)
  
  if (!is.null(outputFolder)) pdf(file = file.path(outputFolder, paste0(name, ".pdf")), width = 10)
  print(
    dat %>%
      ggplot(aes(x = x, y = y, color = density, shape = imputed)) + 
      geom_point(size = 2.5) +
      scale_color_gradient(low = "rosybrown1", high = stopColor) +
      geom_text_repel(data = subset(dat, abs(x) > minXLabel & y > minYLabel), aes(x = x, y = y, label = label), color = "gray24") +
      xlab(label = "log2(fold change)") +
      ylab(label = "-log10(p-value)") +
      guides(color = FALSE, shape = FALSE) +
      theme_classic()
  )
  if (!is.null(outputFolder)) dev.off()
}