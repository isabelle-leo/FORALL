analyzeSilhouettes <- function (eSet,
                                consClust,
                                outputFolder) {
  
  if (!dir.exists(outputFolder)) dir.create(outputFolder, recursive = TRUE)
  
  mat <- eSet %>% exprs()
  
  pdf(file = file.path(outputFolder, "silhouettes.pdf"))
  
  for (method in names(consClust)) {
    methodClust <- consClust[[method]]
    
    distMatrix <- 1 - cor(x = mat, method = method)
    
    for (linkage in names(methodClust)) {
      clust <- methodClust[[linkage]]
      
      for (i in 2:length(clust)) {
        sil <- silhouette(x =clust[[i]]$consensusClass,
                          dmatrix = distMatrix)
        plot(sil, main = paste0(method, " | ", linkage, " | ", i, " clusters"))
      }
    }
  }
  
  dev.off()
}