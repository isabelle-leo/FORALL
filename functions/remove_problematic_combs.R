remove_problematic_combs <- function(mat, min_shared_fields) {
  problematic_combs <- identify_problematic_combs(
    mat = mat, min_shared_fields = min_shared_fields)
  if (!is.null(problematic_combs$row)) {
    mat <- mat[-problematic_combs$row, ]
  }
  if (!is.null(problematic_combs$column)) {
    mat <- mat[, -problematic_combs$column]
  }
  return(mat)
}