correlate <- function (mat,
                       type = "pearson") {

  rcorr_obj <- mat %>%
    t() %>%
    rcorr(type = type)
  
  # convert to long tables
  r <- matToDf(mat = rcorr_obj$r, type = "r")
  # p <- matToDf(mat = rcorr_obj$P, type = "p")
  # n <- matToDf(mat = rcorr_obj$n, type = "n")
  
  df <- r
  # df$p <- p$p
  # df$n <- n$n
  df <- df %>% arrange(desc(r))
  
  return(list(df = df, mat = rcorr_obj$r))
}
