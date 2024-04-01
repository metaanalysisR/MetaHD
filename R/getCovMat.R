getCovMat <- function (sd, cor = NULL) {
  if (is.null(cor))
    cor <- 0
  if (is.data.frame(sd))
    sd <- drop(as.matrix(sd))
  if (is.data.frame(cor))
    cor <- drop(as.matrix(cor))
  if (is.vector(sd))
    sd <- t(sd)
  n <- ncol(sd)
  k <- nrow(sd)
  if (n == 1L)
    return(sd^2)
  if (is.vector(cor)) {
    cor <- if (length(cor) %in% c(1L, k))
      matrix(cor, k, n * (n - 1)/2)
    else if (length(cor) == n * (n - 1)/2)
      matrix(cor, k, n * (n - 1)/2, byrow = TRUE)
    else stop("Dimensions of 'sd' and 'cor' not consistent")
  }
  else if (is.matrix(cor)) {
    if (all(dim(cor) == n) && k == 1L)
      cor <- t(cor[lower.tri(cor)])
    else if (any(dim(cor) != c(k, n * (n - 1)/2)))
      stop("Dimensions of 'sd' and 'cor' not consistent")
  }
  if (any(cor^2 > 1))
    stop("correlations must be between -1 and 1")
  nk <- colnames(sd)
  vcov <- t(sapply(seq(k), function(i) {
    R <- diag(n)
    R[lower.tri(R)] <- cor[i, ]
    R[upper.tri(R)] <- t(R)[upper.tri(R)]
    D <- diag(sd[i, ])
    vechMatrix(D %*% R %*% D)
  }))
  if (k == 1L) {
    vcov <- xpndMatrix(vcov)
    dimnames(vcov) <- list(nk, nk)
  }
  else colnames(vcov) <- vechMatrix(outer(nk, nk, paste, sep = "."))
  vcov
}
