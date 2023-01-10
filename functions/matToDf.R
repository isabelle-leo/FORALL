matToDf <- function (mat,
                     colname = "proteinB",
                     rowname = "proteinA",
                     type = "r",
                     includeDiagonal = FALSE) {
  tmp <- mat
  
  # keep upper triangle only
  tmp[lower.tri(tmp, diag = !includeDiagonal)] <- NA
  
  res <- tmp %>%
    as.data.frame() %>%
    rownames_to_column(rowname) %>%
    gather(key = !!colname, value = !!type, -!!rowname) %>%
    filter(!is.na(.[[type]]))
  
  return(res)
}