#' A Multivariate Meta-Analysis Model for Metabolomics Data
#'
#' The MetaHD function performs a multivariate meta-analysis for combining summary estimates obtained from multiple metabolomic studies by using restricted maximum likelihood estimation.
#' Assuming a meta-analysis is based on N outcomes/metabolites and K studies:
#'
#' @param Y : treatment effect sizes of the outcomes. This should be in the form of a K x N matrix
#' @param Slist : K-dimensional list of N x N matrices representing within-study variances and covariances of the treatment effects
#' @param Psi : N x N matrix representing between-study variances and covariances of the treatment effects. (optional, if not specified this will be estimated internally by "MetaHD" using "estimateBSvar" and "estimateCorMat" functions in "MetaHD" package
#' @param shrinkCor : a logical value indicating whether a shrinkage estimator should be used to estimate between-study correlation matrix. Default is TRUE
#' @param method : estimation method: "fixed" for fixed-effects models,"reml" for random-effects models fitted through restricted maximum likelihood
#' @param bscov : a character vector defining the structure of the random-effects covariance matrix. Among available covariance structures, the user can select "unstructured" to obtain between-study covariance matrix with diagonal elements (variances) estimated using restricted maximul likelihood and off-diagonal elements (co-variances) reflecting the correlations estimated via shrinkage and "diag" (diagonal) for between-study variances as diagonal elements and zero co-variances
#' @param rigls.maxiter : maximum number of iterations of the restricted iterative generalized least square algorithm. Default is set to 1
#' @param impute.na : a logical value indicating whether missing values need to be imputed or not. Default is FALSE
#' @param impute.var : multiplier for replacing the missing variances in Slist.(a large value, default is 10^4)

#' @return A list of objects containing estimate : a N-dimensional vector of the combined estimates, std.err : a N-dimensional vector of the associated standard errors, pVal : a N-dimensional vector of the p-values, I2.stat : I2 statistic
#' @export MetaHD

#' @name MetaHD
#' @useDynLib MetaHD, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

Rcpp::sourceCpp("src/cpp_XtVX.cpp")

MetaHD <- function(Y,Slist,Psi = NULL,shrinkCor = TRUE,method = c("reml","fixed"),bscov = c("unstructured","diag"),rigls.maxiter = 1,impute.na = FALSE,impute.var = 10^4){
  y <- Y
  if (!is.matrix(y)){
    y <- as.matrix(y)
  }
  nay <- is.na(y)
  nall <- sum(!nay)
  N <- ncol(y)
  K <- nrow(y)
  q <- p <- 1
  method <- match.arg(method)
  bscov <- match.arg(bscov)
  if(impute.na){
    data <- low.weight(y,Slist,impute.var)
    y <- data$effects
    Slist <- data$wscovar
    for (k in 1:K) {
      if (!matrixcalc::is.positive.definite(Slist[[k]])) {
        Slist[[k]] <- as.matrix(Matrix::nearPD(Slist[[k]],keepDiag = TRUE,maxit = 500)$mat)
      }
    }
    nay[nay] <-FALSE
  }
  if (any(sapply(Slist, function(x) any(is.na(x))))){
    stop("Error: Slist contains missing values.")
  }
  S <- t(sapply(Slist,vechMatrix))
  rep <- matrix(rep(1,K),nrow = K,ncol = 1)
  nalist <- lapply(1:K, function(j) nay[j,])
  ylist <- lapply(seq(K), function(i) c(t(y[i,]))[!nalist[[i]]])
  X <- matrix(1,nrow = K,ncol = 1)
  Xlist <- lapply(seq(K), function(i) (X[i, , drop = FALSE] %x% diag(N))[!nalist[[i]], , drop = FALSE])
  gls <- glsfit(Xlist, ylist, Slist, onlycoef = FALSE)
  Q <- drop(crossprod(gls$invtUy - gls$invtUX %*% gls$coef))
  df <- nall-N
  if (N > 1L) {
    coef <- matrix(gls$coef, ncol = N, byrow = TRUE)
    indS <- diag(xpndMatrix(seq(N * (N + 1)/2)))
    Q <- c(Q, sapply(seq(N), function(i) sum((y[, i] - X %*%coef[, i])^2/S[, indS[i]])))
    df <- c(df, colSums(!nay, na.rm = TRUE) - p)
  }
  I2 <- pmax((Q - df)/Q * 100, 0)
  if(is.null(Psi)){
    Psi <- diag(0.001, nrow = N, ncol = N)
    if (method == "reml"){
      psi_var <- estimateBSVar(Psi, Xlist, Zlist = NULL, ylist, Slist, nalist, rep, N, q, nall,rigls.maxiter)
      if(bscov == "unstructured"){
        if(impute.na){
          y.imp <- Y
          for (i in 1:N){
            y.imp[,i][is.na(Y[,i])] <- mean(stats::na.omit(Y[,i]))
          }
        }else{
          y.imp <- y
        }
        cormat <- estimateCorMat(y.imp,shrinkCor)
        psi <- getCovMat(sqrt(psi_var),cormat)
      }else if(bscov == "diag"){
        psi <- diag(psi_var, nrow = N , ncol = N)
      }
    }else if (method == "fixed"){
      psi <- diag(0, nrow = N , ncol = N)
    }
  }else{
    psi <- Psi
  }
  A <- matrix(0, ncol = N, nrow = N)
  B <- matrix(0, nrow = N)
  for(k in 1:K) {
    val1 <- solve(Slist[[k]] + psi)
    A <- A + val1
    val2 <- val1 %*% y[k,]
    B <- B + val2
  }
  estimate <- as.numeric(t(solve(A)%*%B))
  std.err <- sqrt(diag(solve(A)))
  zval <- estimate/std.err
  pvalue <- 2 * (1 - stats::pnorm(abs(zval)))
  return(list(estimate = estimate,
              std.err = std.err,
              pVal = pvalue,
              I2.stat = I2))
}

estimateCorMat <- function(Y,shrinkCor = TRUE){
  N <- ncol(Y)
  K <- nrow(Y)
  if (K <= 2){
    cormat <- 0
  }else {
    if (N > K){
      if(shrinkCor){
        cormat <- corpcor::cor.shrink(Y,verbose = FALSE)[1:N,1:N]
      }else{
        cormat <- stats::cor(Y)
      }
    }else{
      cormat <- stats::cor(Y)
    }
  }
  return(cormat)
}

low.weight <- function(Y,Slist,impute.var = 10^4){
  N <- ncol(Y)
  K <- nrow(Y)
  Y[is.na(Y)] <- 0
  S <- t(sapply(Slist,vechMatrix))
  n <- 1
  for (j in 1:ncol(S)) {
    if (j == (((2*N + 3)*n - (n^2) - (2*N))/2)){
      S[,j][is.na(S[,j])] <- max(S[,j],na.rm = TRUE)*impute.var
      n <- n + 1
    }else{
      S[,j][is.na(S[,j])] <- 0
    }
  }
  Slist <- lapply(seq(nrow(S)), function(i) xpndMatrix(S[i,]))
  return(list(effects = Y,
              wscovar = Slist))
}

# Following codes are adapted from library(mixmeta) to estimate between-study variances within MetaHD

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
  cpp_result <- cpp_XtVX(Qlist, Alist, ind1, ind2)
  XtVX <- xpndMatrix(cpp_result)
  XtVy <- sapply(seq(Qlist[[1L]]), function(k) sum(sapply(seq(Qlist), function(i) sum(diag(Alist[[i]][[k]] %*% Blist[[i]])))))
  par <- as.numeric(chol2inv(chol(XtVX)) %*% XtVy)
  Psi <- lapply(lapply(seq_along(q), function(j) par[seq(c(0,
                                                           cumsum(q * N * (q * N + 1)/2))[j] + 1, cumsum(q * N *
                                                                                                           (q * N + 1)/2)[j])]), xpndMatrix)
  Psi <- checkPD(Psi, set.negeigen = sqrt(.Machine$double.eps),force = TRUE, error = FALSE)
  dropList(Psi)
}

checkPD <- function (x, set.negeigen = sqrt(.Machine$double.eps), force = TRUE, error = FALSE, label = "x") {
  x <- getList(x)
  x <- lapply(x, function(mat) {
    eig <- eigen(mat)
    if (any(ind <- eig$values < 0) && error)
      stop(paste("Problems with positive-definiteness in '", label, "'. ", sep = ""))
    if (any(ind) && force) {
      eig$values[ind] <- set.negeigen
      mat <- eig$vectors %*% diag(eig$values, ncol(mat)) %*% t(eig$vectors)
    }
    return(mat)
  })
  dropList(x)
}

reml.newton <- function (Psi, Xlist, Zlist, ylist, Slist, nalist, rep, N, q, nall, const) {
  par <- log(diag(Psi))
  fn <- reml.loglik.fn
  gr <- NULL
  opt <- stats::optim(par = par, fn = fn, gr = gr, Xlist = Xlist,
               Zlist = Zlist, ylist = ylist, Slist = Slist, nalist = nalist,
               rep = rep, N = N, q = q, nall = nall, const = const, method = "BFGS",
               control = list(fnscale=-1, maxit=100, reltol=sqrt(.Machine$double.eps)), hessian = FALSE)
  d <- N*q
  Psi <- diag(exp(opt$par), d)
  return(Psi)
}

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

glsfit <- function (Xlist, ylist, Sigmalist, onlycoef = TRUE){
  Ulist <- lapply(Sigmalist, chol)
  Ulist <- mapply(function(U) {
    if (det(U) == 0) {
      diag_mat <- diag(1e-10, nrow = nrow(U))
      U <- U + diag_mat
    }
    return(U)
  }, Ulist, SIMPLIFY = FALSE)
  invUlist <- lapply(Ulist, function(U) backsolve(U, diag(ncol(U))))
  invtUXlist <- mapply(function(invU, X) crossprod(invU, X),
                       invUlist, Xlist, SIMPLIFY = FALSE)
  invtUylist <- mapply(function(invU, y) crossprod(invU, y),
                       invUlist, ylist, SIMPLIFY = FALSE)
  invtUX <- do.call("rbind", invtUXlist)
  invtUy <- do.call("rbind", invtUylist)
  coef <- as.numeric(qr.solve(invtUX, invtUy))
  if (onlycoef)
    return(coef)
  list(coef = coef, Ulist = Ulist, invUlist = invUlist, invtUXlist = invtUXlist,
       invtUX = invtUX, invtUy = invtUy)
}

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

getSigmaList <- function (Zlist, nalist, Psi, Slist){
  if (is.null(Psi))
    return(Slist)
  Psi <- getList(Psi)
  if (is.null(Zlist))
    return(mapply(function(S, na) S + Psi[[1L]][!na, !na, drop = FALSE], Slist, nalist, SIMPLIFY = FALSE))
  Psi <- getList(Psi)
  lapply(seq_along(Zlist), function(i) sumList(lapply(seq_along(Psi), function(j) blockDiagMat(lapply(Zlist[[i]][[j]], function(x) x %*% Psi[[j]] %*% t(x))))) + Slist[[i]])
}

blockDiagMat <- function (x){
  if (is.matrix(x))
    return(x)
  if (!all(sapply(x, is.matrix))) {
    warning("non-matrix components trasformed in matrices")
    x <- lapply(x, as.matrix)
  }
  if (length(x) == 1L)
    return(x[[1]])
  dim <- t(sapply(x, dim))
  end <- apply(dim, 2, cumsum)
  start <- apply(end, 2, function(x) c(1, x[-length(x)] + 1))
  matind <- array(seq(prod(colSums(dim))), colSums(dim))
  ind <- unlist(lapply(seq(nrow(dim)), function(i) matind[start[i, 1]:end[i, 1], start[i, 2]:end[i, 2]]))
  mat <- matrix(0, sum(dim[, 1]), sum(dim[, 2]))
  mat[ind] <- unlist(x)
  mat
}

sumList <- function (list){
  res <- 0
  for (i in seq(list)){
    res <- res + list[[i]]
  }
  res
}

getList <- function (object){
  if (is.list(object))
    object
  else list(object)
}

dropList <- function (object){
  if (is.list(object) && length(object) == 1L)
    object[[1L]]
  else object
}

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

vechMatrix <- function (mat, diag = TRUE){
  if (!is.matrix(mat)){
    mat <- as.matrix(mat)
  }
  if (diff(dim(mat)) != 0){
    stop("Error: Non-square matrix")
  }
  return(mat[lower.tri(mat, diag = diag)])
}

xpndMatrix <- function (vech){
  dim <- (-1 + sqrt(1 + 8 * length(vech)))/2
  if (dim%%1 != 0L){
    stop("Error: dimension of 'vech' not consistent")
  }
  mat <- matrix(nrow = dim, ncol = dim)
  mat[lower.tri(mat, diag = TRUE)] <- as.matrix(vech)
  mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
  return(mat)
}
