# Some parts of the code are adapted from library(mixmeta)

estimateBSVar <- function (Psi, Xlist, Zlist, ylist, Slist, nalist, rep, N, q, nall,rigls.maxiter) {
  const <- -0.5 * (nall - ncol(Xlist[[1L]])) * log(2 * pi) + sum(log(diag(chol(sumList(lapply(Xlist, crossprod))))))
  Qlist <- getQlist(Zlist, nalist, rep, N, q)
  niter <- 0
  converged <- FALSE
  reltol <- sqrt(.Machine$double.eps)
  while (!converged && niter < rigls.maxiter) {
    old <- unlist(Psi)
    Psi <- rigls.iter(Psi, Qlist, Xlist, Zlist, ylist, Slist, nalist, rep,N,q)
    niter <- niter + 1
    converged <- all(abs(unlist(Psi) - old) < reltol * abs(unlist(Psi) + reltol))
  }
  Psi <- reml.newton(Psi, Xlist, Zlist, ylist, Slist, nalist, rep, N, q, nall, const)
  return(diag(Psi))
}
