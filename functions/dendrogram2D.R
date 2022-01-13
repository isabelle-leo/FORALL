dendrogram2D <- function (proteins,
                          nodeSize = 100,
                          directShow = TRUE,
                          colorScheme,
                          outputFolder = NULL) {
  
  if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder, recursive = TRUE)
  
  # inspirational source:
  # https://github.com/aleferna/BCLandscape/blob/master/Generate_Clustering_and_Network.R
  
  # remove NAs
  eSet <- removeNAsFromESet(proteins)
  
  df <- exprs(eSet)
  colnames(df) <- pData(eSet)$cell_line
  
  # calculate sample dendrogram
  mat.distance <- as.dist(1 - cor(df, method = "pearson"))
  hierarClustering <- hclust(mat.distance, method = "ward.D2")
  data.dendrogram <- as.dendrogram(hierarClustering)
  
  if (!is.null(outputFolder)) pdf(file = file.path(outputFolder, "dendrogram.pdf"), width = 12)
  plot(data.dendrogram)
  if (!is.null(outputFolder)) dev.off()
  
  # textualise dendrogram
  ul <- function(dendro, parent) {
    if (is.list(dendro)) {
      i <- 0
      foreach (y = dendro, .combine=rbind) %do% {
        i <- i + 1 
        child = paste0(parent, ".", i)
        rbind(ul(y, child), c(parent, child, attr(y, "height")))
      }
    } else {
      name <- attr(dendro, "label")
      c(parent, name, 0)
    }
  }
  text.dendrogram <- ul(data.dendrogram, "root")
  
  # define edges
  edges <- data.frame(text.dendrogram, stringsAsFactors = FALSE)
  colnames(edges) <- c("source", "target", "weight")
  
  # find samples in the tree
  lSamples <- !grepl(edges$target, pattern = "^root")
  samples <- edges[lSamples, ]
  edges <- edges[!lSamples, ]
  
  # redefine targets
  for (i in 1:nrow(samples)){
    idx <- edges$target == samples$source[i] 
    edges$target[idx] <- samples$target[i] 
  }
  
  # define nodes
  nodes <- data.frame(name = sort(unique(c(edges$source, edges$target))), stringsAsFactors = FALSE)
  rownames(nodes) <- nodes$name
  
  # get cluster names and assign them
  # clusters <- pData(eSet)$Subtype_Paper
  clusters <- pData(eSet)$hierarCluster_pearson_ward.D2
  names(clusters) <- pData(eSet)$cell_line
  nodes$cluster <- clusters[nodes$name]
  
  nodes$cluster[is.na(nodes$cluster)] <- "Weak"
  nodes$cluster[!nodes$name %in% pData(eSet)$cell_line] <- "Group"
  nodes[nodes$name %in% pData(eSet)$cell_line, "sz"] <- nodeSize
  nodes[!nodes$name %in% pData(eSet)$cell_line, "sz"] <- 2
  
  # modify
  nodes$id <- 0:(nrow(nodes)-1)
  rownames(nodes) <- nodes$name
  edges$source <- nodes[edges$source, "id"]
  edges$target <- nodes[edges$target, "id"]
  edges$weight <- 10 * as.double(edges$weight) + 0.1
  
  # define colors
  colors <- colorScheme$samplesClusters$hierar[1:(length(unique(nodes$cluster)) - 1)]
  xvals <- paste0(paste0('"', c(labels(colors), "Group", "Weak"), '"'), collapse = ',')
  xcols <- paste0(paste0('"', c(colors, "#808080", "#808000"), '"'), collapse = ',')
  xscale <- paste0('d3.scaleOrdinal().domain([', xvals, ']).range([', xcols, "])")
  
  # build the force network
  network <- forceNetwork(Links        = edges,
                          Nodes        = nodes,
                          NodeID       = "name",
                          Group        = "cluster", 
                          legend       = TRUE, 
                          zoom         = TRUE,  
                          Value        = "weight",
                          linkDistance = JS('function(d){return d.value}') ,
                          Nodesize     = "sz",
                          fontSize     = 20,
                          opacity      = 1 ,
                          colourScale  = xscale)
  
  # show it
  if (directShow) print(network)
  
  # save it
  if (!is.null(outputFolder)) {
    saveNetwork(network = network,
                file = file.path(getwd(), outputFolder, "dendrogram2D.html"),
                selfcontained = TRUE)
  }
}