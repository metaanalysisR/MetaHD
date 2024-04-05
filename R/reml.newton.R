# Code adapted from library(mixmeta)

reml.newton <- function (Psi, Xlist, Zlist, ylist, Slist, nalist, rep, N, q, nall, const) {
  par <- log(diag(Psi))
  fn <- reml.loglik.fn
  gr <- NULL
  opt <- optim(par = par, fn = fn, gr = gr, Xlist = Xlist,
               Zlist = Zlist, ylist = ylist, Slist = Slist, nalist = nalist,
               rep = rep, N = N, q = q, nall = nall, const = const, method = "BFGS",
               control = list(fnscale=-1, maxit=100, reltol=sqrt(.Machine$double.eps)), hessian = FALSE)
  d <- N*q
  Psi <- diag(exp(opt$par), d)
  return(Psi)
}
