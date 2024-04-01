vechMatrix <- function (mat, diag = TRUE){
  if (!is.matrix(mat)){
    mat <- as.matrix(mat)
  }
  if (diff(dim(mat)) != 0){
    stop("Error: Non-square matrix")
  }
  return(mat[lower.tri(mat, diag = diag)])
}
