## Necessary libraries
# library(corpcor)
# library(Matrix)
# library(matrixcalc)

#' The function performs a multivariate meta-analysis for combining summary estimates obtained from multiple metabolomic studies by using restricted maximum likelihood estimation. 
#' Assuming a meta-analysis is based on N outcomes/metabolites and K studies:
#' @param Y : treatment effect sizes of the outcomes. This should be in the form of a K x N matrix
#' @param Slist : K-dimensional list of N x N matrices representing within-study variances and covariances of the treatment effects
#' @param Psi : N x N matrix representing between-study variances and covariances of the treatment effects. (optional, if not specified this will be estimated internally by "MetaHD" using "estimateBSvar" and "estimateCorMat" functions in "MetaHD" package
#' @param shrinkCor : a logical value indicating whether a shrinkage estimator should be used to estimate between-study correlation matrix. Default is TRUE
#' @param method : estimation method: "fixed" for fixed-effects models,"reml" for random-effects models fitted through restricted maximum likelihood
#' @param bscov : a character vector defining the structure of the random-effects covariance matrix. Among available covariance structures, the user can select "unstructured" to obtain between-study covariance matrix with diagonal elements (variances) estimated using restricted maximul likelihood and off-diagonal elements (co-variances) reflecting the correlations estimated via shrinkage and "diag" (diagonal) for between-study variances as diagonal elements and zero co-variances
#' @param rigls.maxiter : maximum number of iterations of the restricted iterative generalized least square algorithm. Default is set to 5
#' @param impute.na : a logical value indicating whether missing values need to be imputed or not. Default is FALSE
#' @param impute.var : multiplier for replacing the missing variances in Slist.(a large value, default is 10^4

#' @return A list of objects containing:
#' estimate : a N-dimensional vector of the combined estimates
#' std.err : a N-dimensional vector of the associated standard errors
#' I2.stat : I2 statistic

MetaHD <- function(Y,Slist,Psi = NULL,shrinkCor = TRUE,method = c("reml","fixed"),bscov = c("unstructured","diag"),rigls.maxiter = 5,impute.na = FALSE,impute.var = 10^4){
  # ENSURING Y IS A MATRIX
  y <- Y
  if (!is.matrix(y)){
    y <- as.matrix(y)
  }
  # SET DIMENSIONS
  nay <- is.na(y)
  nall <- sum(!nay)
  N <- ncol(y) # NO.OF METABOLITES
  K <- nrow(y) # NO.OF STUDIES
  q <- p <- 1
  # USER INPUTS
  method <- match.arg(method)
  bscov <- match.arg(bscov)
  # IMPUTE MISSING
  if(impute.na){
    data <- low.weight(y,Slist,impute.var)
    y <- data$effects
    Slist <- data$wscovar
    # CHECK POSITIVE DEFINITENESS
    for (k in 1:K) {
      if (!is.positive.definite(Slist[[k]])) {
        Slist[[k]] <- as.matrix(nearPD(Slist[[k]],keepDiag = TRUE,maxit = 500)$mat)
      }
    }
    nay[nay] <-FALSE
  }
  # IF Slist CONTAINS MISSINGS GENERATE AN ERROR
  if (any(sapply(Slist, function(x) any(is.na(x))))){
    stop("Error: Slist contains missing values.")
  }
  # CONVERTS EACH MATRIX IN Slist INTO A VECTOR BY EXTRACTING LOWER TRIANGULAR ELEMENTS AND TRANSPOSES THE RESULT
  S <- t(sapply(Slist,vechMatrix))
  rep <- matrix(rep(1,K),nrow = K,ncol = 1)
  nalist <- lapply(1:K, function(j) nay[j,])
  # TRANSFORM y, VECTORIZING THEM BY ROW
  ylist <- lapply(seq(K), function(i) c(t(y[i,]))[!nalist[[i]]])
  X <- matrix(1,nrow = K,ncol = 1)
  # TRANSFORM X
  Xlist <- lapply(seq(K), function(i) (X[i, , drop = FALSE] %x% diag(N))[!nalist[[i]], , drop = FALSE])
  # COMPUTE I2 STATISTIC
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
    # SET INTIAL Psi IF NULL
    Psi <- diag(0.001, nrow = N, ncol = N)
    # ESTIMATE PSI 
    if (method == "reml"){
      psi_var <- estimateBSVar(Psi, Xlist, Zlist = NULL, ylist, Slist, nalist, rep, N, q, nall,rigls.maxiter)
      if(bscov == "unstructured"){
        # IMPUTE Y FOR ESTIMATING BETWEEN-STUDY CORRELATIONS
        if(impute.na){
          y.imp <- Y
          for (i in 1:N){
            y.imp[,i][is.na(Y[,i])] <- mean(na.omit(Y[,i]))
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
    # KNOWN PSI : USER INPUT
    psi <- Psi
  }
  # OBTAIN COMBINED ESTIMATES AND STANDARD ERRORS
  A <- matrix(0, ncol = N, nrow = N)
  B <- matrix(0, nrow = N)
  for(k in 1:K) {
    val1 <- solve(Slist[[k]] + psi)
    A <- A + val1
    val2 <- val1 %*% y[k,]
    B <- B + val2
  }
  return(list(estimate = as.numeric(t(solve(A)%*%B)),
              std.err = sqrt(diag(solve(A))),
              I2.stat = I2))
}

