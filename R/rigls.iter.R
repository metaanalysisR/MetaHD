rigls.iter <- function (Psi, Qlist, Xlist, Zlist, ylist, Slist, nalist, rep, N, q){
  Sigmalist <- getSigmaList(Zlist, nalist, Psi, Slist)
  gls <- glsfit(Xlist, ylist, Sigmalist, onlycoef = FALSE)
  invSigmalist <- lapply(gls$invUlist, tcrossprod)
  invtXinvSigmaX <- solve(crossprod(gls$invtUX))
  flist <- mapply(function(y, S, X) tcrossprod(y - X %*% gls$coef) -
                    S + X %*% invtXinvSigmaX %*% t(X), ylist, Slist, Xlist, SIMPLIFY = FALSE)
  Alist <- lapply(seq(Qlist), function(i) lapply(seq(Qlist[[1L]]), function(k) Qlist[[i]][[k]] %*% invSigmalist[[i]]))
  Blist <- lapply(seq(flist), function(i) flist[[i]] %*% invSigmalist[[i]])
  ind1 <- unlist(lapply(seq(Qlist[[1L]]), ":", length(Qlist[[1L]])))
  ind2 <- rep(seq(Qlist[[1L]]), rev(seq(Qlist[[1L]])))
  XtVX <- xpndMatrix(sapply(seq(ind1), function(k) sum(sapply(seq(Qlist),
                                                           function(i) sum(diag(Alist[[i]][[ind1[k]]] %*% Alist[[i]][[ind2[k]]]))))))
  XtVy <- sapply(seq(Qlist[[1L]]), function(k) sum(sapply(seq(Qlist),
                                                          function(i) sum(diag(Alist[[i]][[k]] %*% Blist[[i]])))))
  par <- as.numeric(chol2inv(chol(XtVX)) %*% XtVy)
  Psi <- lapply(lapply(seq_along(q), function(j) par[seq(c(0,
                                                           cumsum(q * N * (q * N + 1)/2))[j] + 1, cumsum(q * N *
                                                                                                           (q * N + 1)/2)[j])]), xpndMatrix)
  Psi <- checkPD(Psi, set.negeigen = sqrt(.Machine$double.eps),force = TRUE, error = FALSE)
  dropList(Psi)
}
