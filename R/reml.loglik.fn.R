# Code adapted from library(mixmeta)

reml.loglik.fn <- function (par, Xlist, Zlist, ylist, Slist, nalist, rep, N, q, nall, const){
  d <- N*q
  Psi <- diag(exp(par), d)
  Sigmalist <- getSigmaList(Zlist, nalist, Psi, Slist)
  gls <- glsfit(Xlist, ylist, Sigmalist, onlycoef = FALSE)
  res <- -0.5 * (crossprod(gls$invtUy - gls$invtUX %*% gls$coef))
  det1 <- -sum(sapply(gls$Ulist, function(U) sum(log(diag(U)))))
  tXWXtot <- sumList(lapply(gls$invtUXlist, crossprod))
  det2 <- -sum(log(diag(chol(tXWXtot))))
  as.numeric(const + det1 + det2 + res)
}
