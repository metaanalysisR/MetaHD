# Code adapted from library(mixmeta)

checkPD <- function (x, set.negeigen = sqrt(.Machine$double.eps), force = TRUE, error = FALSE, label = "x") {
  x <- getList(x)
  x <- lapply(x, function(mat) {
    eig <- eigen(mat)
    if (any(ind <- eig$values < 0) && error)
      stop(paste("Problems with positive-definiteness in '",
                 label, "'. ", sep = ""))
    if (any(ind) && force) {
      eig$values[ind] <- set.negeigen
      mat <- eig$vectors %*% diag(eig$values, ncol(mat)) %*% t(eig$vectors)
    }
    return(mat)
  })
  dropList(x)
}
