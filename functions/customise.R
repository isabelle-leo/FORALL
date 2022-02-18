customise <- function (expressions) {
  
  # take median of technical SEM replicates
  sems <- pData(expressions) %>%
    filter(Cell.Line.R == "SEM")
  
  sems.data <- expressions[, sems$position] %>%
    exprs()
  
  sems.median <- rowMedians(sems.data, na.rm = TRUE)
  
  # now remove all but one of the sem columns
  remove.these <- sems$position[2:length(sems$position)]
  expressions.updated <- expressions[, setdiff(colnames(expressions), remove.these)]
  
  # write medians to remaining SEM column
  exprs(expressions.updated)[, sems$position[[1]]] <- sems.median
  
  # update pData
  pData(expressions.updated)[sems$position[[1]], "cell_line"] <- "SEM median"
  
  
  # mark the biological duplicates
  # already done in meta file
  
  # remove the noPS SEM cell line
  expressions.updated <- expressions.updated %>%
    .[, pData(.)$Cell.Line.R != "SEM.BR.NOPS"]
  
  return(expressions.updated)
}