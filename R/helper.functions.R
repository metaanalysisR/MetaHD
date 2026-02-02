########################### MetaHD core function ##############################
.MetaHD_core <- function(y, Slist, Psi, method, bscov, est.wscor, shrinkCor, optim.algorithm, optim.maxiter, rigls.iter, initPsi, impute.na, impute.var, N, K, p, q) {
  nay <- is.na(y)
  nall <- sum(!nay)
  if (est.wscor) {
    WSVar <- dimCheck(Slist, K, N)
    wscormat <- estimateCorMat(y, shrinkCor, impute.na)
    Slist <- lapply(seq_len(K),function(k) getCovMat(sqrt(WSVar[k, ]), wscormat))
  }
  if (method == "multi" && !is.list(Slist) && !est.wscor) {
    stop("Require 'Slist' to be a K-dimensional list of N x N matrices representing within-study variances and covariances of the treatment effects. If within-study correlations are not available, set est.wscor = TRUE.")
  }
  if (impute.na) {
    imputed <- imputeNA(y, Slist, impute.var)
    y <- imputed$effects
    Slist <- lapply(imputed$wscovar, ensurePD)
    nay <- is.na(y)
  }
  if ((method == "REM" || method == "FEM") && !impute.na) {
    keep <- is.finite(y) & is.finite(Slist)
    if (!all(keep)) {
      warning("Dropping ", sum(!keep), " studies with missing estimates and variances")
    }
    y <- as.matrix(y[keep])
    Slist <- as.matrix(Slist[keep])
    K <- nrow(y)
  }
  if (any(is.na(y)) || any(sapply(Slist, function(x) any(is.na(x))))) {
    stop("Error: Y and/or Slist contains missing values.")
  }
  S <- t(sapply(Slist,vechMatrix))
  rep <- matrix(rep(1,K),nrow = K,ncol = 1)
  ylist <- lapply(seq(K), function(i) c(t(y[i,])))
  X <- matrix(1,nrow = K,ncol = 1)
  Xlist <- rep(list(diag(N)), K)
  if(is.null(Psi)){
    if(is.null(initPsi)){
      Psi <- diag(1, nrow = N, ncol = N)
    }else{
      Psi <- as.matrix(initPsi)
    }
    if (method != "FEM") {
      psi_var <- estimateBSVar(Psi, Xlist, Zlist = NULL, ylist, Slist, rep, N, K, q, nall, optim.algorithm, rigls.iter, optim.maxiter)
    }
    if (method == "multi"){
      if(bscov == "unstructured"){
        cormat <- estimateCorMat(y,shrinkCor,impute.na)
        psi <- getCovMat(sqrt(psi_var),cormat)
      }else if(bscov == "diag"){
        psi <- diag(psi_var, nrow = N , ncol = N)
      }else if(bscov == "none"){
        psi <- diag(0, nrow = N , ncol = N)
      }
    }else if (method == "REM"){
      psi <- diag(psi_var, nrow = N , ncol = N)
    }else if (method == "FEM"){
      psi <- diag(0, nrow = N , ncol = N)
    }
  }else{
    psi <- as.matrix(Psi)
  }
  A <- matrix(0, nrow = N, ncol = N)
  B <- matrix(0, nrow = N)
  for (k in seq_along(Slist)) {
    Sk_psi <- Slist[[k]] + psi
    R <- tryCatch(chol(Sk_psi), error = function(e) chol(ensurePD(Sk_psi)))
    Vinv <- chol2inv(R)
    A <- A + Vinv
    B <- B + Vinv %*% y[k, ]
  }
  Ainv <- chol2inv(chol(A))
  est <- as.numeric(Ainv %*% B)
  se  <- sqrt(diag(Ainv))
  p   <- 2 * (1 - pnorm(abs(est / se)))
  i2 <- i2Stat(X, y, Xlist, ylist, Slist, S, nay, N, nall, p)
  return(list(estimate = est,
              std.err = se,
              pVal = p,
              I2.stat = i2))
}

##################### Divide-and-conquer core function ########################

.MetaHD_DC_core <- function(vars_in_group, y, Slist, Psi, initPsi, method, bscov, est.wscor, shrinkCor, optim.algorithm, optim.maxiter, rigls.iter, impute.na, impute.var, K, p, q) {
  N_dc <- length(vars_in_group)
  y_dc <- y[, vars_in_group, drop = FALSE]
  Slist_dc <- if (is.list(Slist)) {
    lapply(Slist, function(s) s[vars_in_group, vars_in_group, drop = FALSE])
  } else {
    Slist[, vars_in_group, drop = FALSE]
  }
  Psi_dc <- if (is.null(Psi)) NULL else Psi[vars_in_group, vars_in_group, drop = FALSE]
  initPsi_dc <- if (is.null(initPsi)) NULL else initPsi[vars_in_group, vars_in_group, drop = FALSE]
  return(
    .MetaHD_core(y = y_dc, Slist = Slist_dc, Psi = Psi_dc, method = method, bscov = bscov, est.wscor = est.wscor, shrinkCor = shrinkCor, optim.algorithm = optim.algorithm, optim.maxiter = optim.maxiter, rigls.iter = rigls.iter, initPsi = initPsi_dc, impute.na = impute.na, impute.var = impute.var, N = N_dc, K = K, p = p, q = q)
  )
}

############################ estimate correlations ############################

estimateCorMat <- function(Y,shrinkCor = TRUE,impute.na = FALSE){
  N <- ncol(Y)
  K <- nrow(Y)
  if(impute.na){
    for (i in 1:N){
      Y[,i][is.na(Y[,i])] <- mean(na.omit(Y[,i]))
    }
  }else{
    Y <- Y
  }
  if (K <= 2){
    cormat <- 0
  }else {
    if (N > K){
      if(shrinkCor){
        cormat <- cor.shrink(Y,verbose = FALSE)[1:N,1:N]
      }else{
        cormat <- cor(Y)
      }
    }else{
      cormat <- cor(Y)
    }
  }
  return(cormat)
}

############################ impute missing values ############################

imputeNA <- function(Y,Slist,impute.var = 10^4){
  N <- ncol(Y)
  K <- nrow(Y)
  Y[is.na(Y)] <- 0
  if (is.list(Slist)){
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
  }else{
    for (i in 1:N) {
      Slist[,i][is.na(Slist[,i])] <- max(Slist[,i],na.rm = TRUE)*impute.var
    }
  }
  return(list(effects = Y,
              wscovar = Slist))
}

############################ estimate I2 statistic ############################

i2Stat <- function(X, y,  Xlist, ylist, Slist, S, nay, N, nall, p){
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
  return(I2)
}

########################### other helper functions ############################

dimCheck <- function(Slist, K, N) {
  if (is.list(Slist)) {
    stop("Require within-study variances to be in the form of a K x N matrix, where K is the number of studies and N is the number of outcomes")
  }
  WSVar <- as.matrix(Slist)
  if (nrow(WSVar) != K || ncol(WSVar) != N) {
    stop("Dimensions of within-study variance matrix and Y are not consistent")
  }
  return(WSVar)
}

ensurePD <- function(M, maxit = 1000) {
  if (!is.positive.definite(M)) {
    M <- as.matrix(nearPD(M, keepDiag = TRUE, maxit = maxit)$mat)
  }
  M
}

# Following codes are adapted from 'mixmeta' R package to estimate between-study variances within MetaHD function

estimateBSVar <- function (Psi, Xlist, Zlist, ylist, Slist, rep, N, K, q, nall, optim.algorithm, rigls.iter, optim.maxiter) {
  const <- -0.5 * (nall - ncol(Xlist[[1L]])) * log(2 * pi) + sum(log(diag(chol(sumList(lapply(Xlist, crossprod))))))
  if(optim.algorithm == "hybrid"){
    Qlist <- getQlist(Zlist, rep, N, K, q)
    niter <- 0
    converged <- FALSE
    reltol <- sqrt(.Machine$double.eps)
    while (!converged && niter < rigls.iter) {
      old <- unlist(Psi)
      Psi <- rigls.iter(Psi, Qlist, Xlist, Zlist, ylist, Slist, rep, N, q)
      niter <- niter + 1
      converged <- all(abs(unlist(Psi) - old) < reltol * abs(unlist(Psi) + reltol))
    }
    Psi <- reml.newton(Psi, Xlist, Zlist, ylist, Slist, rep, N, q, nall, const, optim.algorithm, optim.maxiter)
  }else if(optim.algorithm == "L-BFGS-B"){
    Psi <- reml.newton(Psi, Xlist, Zlist, ylist, Slist, rep, N, q, nall, const, optim.algorithm, optim.maxiter)
  }else if(optim.algorithm == "BOBYQA"){
    Psi <- reml.bobyqa(Psi, Xlist, Zlist, ylist, Slist, rep, N, q, nall, const, optim.maxiter)
  }
  return(diag(Psi))
}

rigls.iter <- function (Psi, Qlist, Xlist, Zlist, ylist, Slist, rep, N, q){
  Sigmalist <- lapply(Slist, function(x) x + Psi)
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

reml.newton <- function (Psi, Xlist, Zlist, ylist, Slist, rep, N, q, nall, const, optim.algorithm, optim.maxiter) {
  par <- log(diag(Psi))
  fn <- reml.loglik.fn
  gr <- NULL
  if(optim.algorithm == "hybrid"){
    opt <- optim(par = par, fn = fn, gr = gr, Xlist = Xlist,
                 Zlist = Zlist, ylist = ylist, Slist = Slist,
                 rep = rep, N = N, q = q, nall = nall, const = const, method = "BFGS",
                 control = list(fnscale=-1, maxit=optim.maxiter, reltol=sqrt(.Machine$double.eps)), hessian = FALSE)
  }else{
    opt <- optim(par = par, fn = fn, gr = gr, Xlist = Xlist,
                 Zlist = Zlist, ylist = ylist, Slist = Slist,
                 rep = rep, N = N, q = q, nall = nall, const = const, method = "L-BFGS-B",
                 control = list(fnscale=-1, maxit=optim.maxiter, factr = 1e+07, pgtol = 0), hessian = FALSE)
  }
  converged <- opt$convergence == 0
  if (!is.null(converged) && !converged) {
    print(opt$message)
  }
  d <- N*q
  Psi <- diag(exp(opt$par), d)
  return(Psi)
}

reml.bobyqa <- function(Psi, Xlist, Zlist, ylist, Slist, rep, N, q, nall, const, optim.maxiter) {
  par <- log(diag(Psi))
  fn <- function(par, Xlist, Zlist, ylist, Slist, rep, N, q, nall, const) {
    -reml.loglik.fn(par, Xlist, Zlist, ylist, Slist, rep, N, q, nall, const)
  }
  opt <- bobyqa(x0 = par, fn = fn, Xlist = Xlist, Zlist = Zlist, ylist = ylist,
                Slist = Slist, rep = rep, N = N, q = q, nall = nall, const = const,
                nl.info = FALSE,control = list(maxeval = optim.maxiter, ftol_rel = 1e-8, xtol_rel = 1e-8))
  converged <- opt$convergence > 0
  if (!is.null(converged) && !converged) {
    print(opt$message)
  }
  d <- N * q
  Psi <- diag(exp(opt$par), d)
  return(Psi)
}

reml.loglik.fn <- function (par, Xlist, Zlist, ylist, Slist, rep, N, q, nall, const){
  d <- N*q
  Psi <- diag(exp(par), d)
  Sigmalist <- lapply(Slist, function(x) x + Psi)
  gls <- tryCatch({
    glsfit(Xlist, ylist, Sigmalist, onlycoef = FALSE)
  }, error = function(e) {
    message("Skipping iteration: ", e$message)
    return(NULL)
  })
  if (is.null(gls)) {
    return(Inf)
  }
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

getQlist <- function (Zlist, rep, N, K, q){
  if (is.null(Zlist))
    Zlist <- lapply(1:K, function(id) list(list(diag(N))))
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
