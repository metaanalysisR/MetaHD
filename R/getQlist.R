# Code adapted from library(mixmeta)

getQlist <- function (Zlist, nalist, rep, N, q){
  if (is.null(Zlist))
    Zlist <- lapply(nalist, function(na) list(list(diag(length(na))[!na, , drop = FALSE])))
  Qlist <- lapply(seq(Zlist), function(i) {
    Zexp <- do.call(cbind, lapply(seq(length(q)), function(j) blockDiagMat(Zlist[[i]][[j]])))
    do.call(c, lapply(seq_along(q), function(j) {
      rows <- vechMatrix(row(diag(q[j] * N)))
      cols <- vechMatrix(col(diag(q[j] * N)))
      start <- c(0, cumsum(q * N * rep[i, ]))[j]
      lapply(seq(rows), function(t) {
        sumList(lapply(seq(rep[i, j]), function(r) {
          ind1 <- start + (r - 1) * (q[j] * N) + rows[t]
          ind2 <- start + (r - 1) * (q[j] * N) + cols[t]
          if (ind1 == ind2)
            tcrossprod(Zexp[, ind1])
          else tcrossprod(Zexp[, ind1], Zexp[, ind2]) +
            tcrossprod(Zexp[, ind2], Zexp[, ind1])
        }))
      })
    }))
  })
  Qlist
}
