#' A Multivariate Meta-Analysis Model for High-Dimensional Data
#'
#' @description
#' The MetaHD function performs a multivariate meta-analysis for high-dimensional data, combining summary estimates obtained from multiple studies by using restricted maximum likelihood estimation. In its default settings, the function fits the fastMetaHD model, which provides a memory-efficient and computationally faster implementation of the MetaHD methodology.
#' Assuming a meta-analysis is based on \eqn{N} outcomes and \eqn{K} studies:
#'
#' @usage MetaHD(
#'   Y,
#'   Slist,
#'   Psi = NULL,
#'   method = c("multi","REM","FEM"),
#'   bscov = c("unstructured","diag","none"),
#'   useDivideConquer = FALSE,
#'   parallel = FALSE,
#'   est.wscor = FALSE,
#'   shrinkCor = TRUE,
#'   impute.na = FALSE,
#'   optim.algorithm = c("BOBYQA","hybrid","L-BFGS-B"),
#'   optim.maxiter = 2000,
#'   rigls.iter = 1,
#'   initPsi = NULL,
#'   impute.var = 10^4
#' )
#'
#' @param Y treatment effect sizes of the outcomes. This should be in the form of a K x N matrix.
#' @param Slist A K-dimensional list of N × N matrices representing within-study variances and covariances of the treatment effects. If within-study correlations are not available, provide the associated variances of the treatment effects as a K × N matrix and set est.wscor = TRUE. For method = "REM" or method = "FEM", provide the associated variances of the treatment effects as a K × N matrix.
#' @param Psi N x N matrix representing between-study variances and covariances of the treatment effects. (optional, if not specified this will be estimated internally by "MetaHD" using "estimateBSvar" and "estimateCorMat" functions in "MetaHD" package).
#' @param method estimation method: "multi" for multivarite meta-analysis model fitted through restricted maximum likelihood estimation where the between-study covariance structure can be selected via 'bscov', "REM" for univariate random-effects model fitted through restricted maximum likelihood estimation and "FEM" for univariate fixed-effects model.
#' @param bscov a character vector defining the structure of the random-effects covariance matrix. Among available covariance structures, the user can select "unstructured" to obtain between-study covariance matrix with diagonal elements (variances) estimated using restricted maximum likelihood and off-diagonal elements (co-variances) reflecting the correlations estimated via shrinkage, "diag" (diagonal) for between-study variances as diagonal elements and zero co-variances, and "none" for zero between-study variances and co-variances.
#' @param useDivideConquer a logical value indicating whether to use the divide-and-conquer implementation of the fastMetaHD model. This option is used only when method = "multi". Default is FALSE.
#' @param parallel a logical value indicating whether to enable parallel computation for the divide-and-conquer approach. Default is \code{FALSE}. See also Details.
#' @param est.wscor a logical value indicating whether the within-study correlation matrix needs to be estimated or not. Default is \code{FALSE}.
#' @param shrinkCor a logical value indicating whether a shrinkage estimator should be used to estimate within- or between-study correlation matrix. \code{TRUE}.
#' @param impute.na a logical value indicating whether missing values need to be imputed or not. Default is \code{FALSE}.
#' @param optim.algorithm specifies the algorithm used to maximize the restricted log-likelihood function for estimating between-study variances. The default algorithm is "BOBYQA", which offers derivative-free, bound-constrained optimization by iteratively constructing a quadratic approximation of the objective function. The "hybrid" option performs up to rigls.iter iterations of the RIGLS algorithm, followed by quasi-Newton (BFGS algorithm) iterations until convergence. If rigls.iter is set to zero, only the quasi-Newton method (BFGS algorithm) is used for estimation. The "L-BFGS-B" algorithm is a limited-memory version of the BFGS quasi-Newton method, which supports box constraints, allowing each variable to have specified lower and/or upper bounds.
#' @param optim.maxiter maximum number of iterations in methods involving optimization procedures.
#' @param rigls.iter number of iterations of the restricted iterative generalized least square algorithm (RIGLS) when used in the initial phase of hybrid optimization procedure. Default is set to 1.
#' @param initPsi N x N diagonal matrix representing the starting values of the between-study variances to be used in the optimization procedures. If not specified, the starting values in Psi default to a diagonal matrix with variances set to 1.
#' @param impute.var multiplier for replacing the missing variances in Slist.(a large value, default is 10^4).
#'
#' @details
#' If \code{parallel = TRUE}, the divide-and-conquer approach may be evaluated in parallel. Parallel execution is implemented using the \code{future} R package.
#'
#' On Windows, users must set a future plan (e.g., \code{future::plan(future::multisession, workers = ncores)}) before calling \code{MetaHD()} in order to enable parallel computation.
#'
#' On Linux and macOS, users may alternatively use \code{future::plan(future::multicore, workers = ncores)}.
#'
#' If no future plan is set, or if \code{parallel = FALSE}, computations are performed sequentially.
#'
#' @return A list of objects containing :
#' \itemize{
#'   \item \code{estimate}: An \eqn{N}-dimensional vector of the combined estimates.
#'   \item \code{std.err}: An \eqn{N}-dimensional vector of the associated standard errors.
#'   \item \code{pVal}: An \eqn{N}-dimensional vector of the \eqn{p}-values.
#'   \item \code{I2.stat}: \eqn{I^2} statistics.
#' }
#'
#' @examples
#' # CREATE INPUT DATA
#' input_data <- MetaHDInput(realdata)
#' Y <- input_data$Y
#' Slist <- input_data$Slist
#'
#' N <- ncol(Y)
#' K <- nrow(Y)
#'
#' Smat <- matrix(0, nrow = K, ncol = N)
#' for (i in 1:K) {
#'  Smat[i, ] <- diag(Slist[[i]])
#' }
#'
#' # MULTIVARIATE RANDOM-EFFECTS META-ANALYSIS 
#' model <- MetaHD(Y = Y, Slist = Slist, method = "multi")
#' model$estimate
#' model$pVal
#'
#' # UNIVARIATE RANDOM-EFFECTS META-ANALYSIS
#' model <- MetaHD(Y = Y, Slist = Smat, method = "REM")
#' model$estimate
#' model$pVal
#'
#' # UNIVARIATE FIXED-EFFECTS META-ANALYSIS
#' model <- MetaHD(Y = Y, Slist = Smat, method = "FEM")
#' model$estimate
#' model$pVal
#'
#' @references
#' Liyanage JC, Prendergast L, Staudte R, De Livera AM (2024).
#' \emph{MetaHD: a multivariate meta-analysis model for metabolomics data}.
#' Bioinformatics, 40(7), btae470.
#' \doi{10.1093/bioinformatics/btae470}
#' @references
#' Powell MJ (2009).
#' \emph{The BOBYQA algorithm for bound constrained optimization without derivatives}.
#' Cambridge NA Report NA2009/06, University of Cambridge, 26, 26–46.
#' @references
#' Sera F, Armstrong B, Blangiardo M, et al. (2019).
#' \emph{An extended mixed-effects framework for meta-analysis}.
#' Statistics in Medicine, 38, 5429--5444.
#' @references
#' Schäfer J, Strimmer K (2005).
#' \emph{A shrinkage approach to large-scale covariance estimation and implications for functional genomics}.
#' Statistical Applications in Genetics and Molecular Biology, 4, 32.
#'
#' @export MetaHD
#'
#' @name MetaHD
#' @useDynLib MetaHD, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom matrixcalc is.positive.definite
#' @importFrom Matrix nearPD
#' @importFrom stats na.omit cor optim pnorm as.dist hclust
#' @importFrom utils capture.output
#' @importFrom dynamicTreeCut cutreeDynamic
#' @importFrom nloptr bobyqa
#' @importFrom corpcor cor.shrink
#' @importFrom future.apply future_lapply
NULL

sourceCpp("src/cpp_XtVX.cpp")

MetaHD <- function(Y,Slist,Psi = NULL,method = c("multi","REM","FEM"),bscov = c("unstructured", "diag", "none"),useDivideConquer = FALSE,parallel = FALSE,est.wscor = FALSE,shrinkCor = TRUE,impute.na = FALSE,optim.algorithm = c("BOBYQA","hybrid","L-BFGS-B"),optim.maxiter = 2000,rigls.iter = 1,initPsi = NULL,impute.var = 10^4){
  method <- match.arg(method)
  bscov <- match.arg(bscov)
  optim.algorithm <- match.arg(optim.algorithm)
  y <- Y
  if (!is.matrix(y)){
    y <- as.matrix(y)
  }
  N <- ncol(y)
  K <- nrow(y)
  q <- p <- 1
  if (is.null(rownames(y)))
    rownames(y) <- paste0("Study", seq_len(K))
  if (is.null(colnames(y)))
    colnames(y) <- paste0("Var", seq_len(N))
  var_names <- colnames(y)
  study_names <- rownames(y)
  if (is.list(Slist)) {
    Slist <- lapply(Slist, function(Sk) {
      if (is.null(rownames(Sk)))
        rownames(Sk) <- var_names
      if (is.null(colnames(Sk)))
        colnames(Sk) <- var_names
      Sk
    })
  } else {
    if (is.null(rownames(Slist)))
      rownames(Slist) <- study_names
    if (is.null(colnames(Slist)))
      colnames(Slist) <- var_names
  }
  if (!is.null(Psi) && (is.null(colnames(Psi)) || is.null(rownames(Psi)))) {
    colnames(Psi) <- rownames(Psi) <- var_names
  }
  if (!is.null(initPsi) && (is.null(colnames(initPsi)) || is.null(rownames(initPsi)))) {
    colnames(initPsi) <- rownames(initPsi) <- var_names
  }
  if (method == "multi") {
    if (useDivideConquer) {
      cormat <- estimateCorMat(y, shrinkCor, impute.na)
      distmat  <- as.dist(1 - abs(cormat))
      hc <- hclust(distmat, method = "complete")
      clusters <- {
        invisible(capture.output(
          res <- cutreeDynamic(hc,distM = as.matrix(distmat)),
          type = "output"
        ))
        res
      }
      names(clusters) <- hc$labels
      groups <- split(names(clusters), clusters)
      grouped_vars <- unlist(groups, use.names = FALSE)
      if (parallel) {
        res_list <- future.apply::future_lapply(groups, function(vars_in_group) {
          .MetaHD_DC_core(vars_in_group = vars_in_group,y = y,Slist = Slist,Psi = Psi,initPsi = initPsi,method = method,bscov = bscov,est.wscor = est.wscor,shrinkCor = shrinkCor,optim.algorithm = optim.algorithm,optim.maxiter = optim.maxiter,rigls.iter = rigls.iter,impute.na = impute.na,impute.var = impute.var,K = K,p = p,q = q)
        },future.seed = TRUE)
      } else {
        res_list <- lapply(groups, function(vars_in_group) {
          .MetaHD_DC_core(vars_in_group = vars_in_group,y = y,Slist = Slist,Psi = Psi,initPsi = initPsi,method = method,bscov = bscov,est.wscor = est.wscor,shrinkCor = shrinkCor,optim.algorithm = optim.algorithm,optim.maxiter = optim.maxiter,rigls.iter = rigls.iter,impute.na = impute.na,impute.var = impute.var,K = K,p = p,q = q)
        })
      }
      est_DC  <- unname(unlist(lapply(res_list, `[[`, "estimate")))
      se_DC   <- unname(unlist(lapply(res_list, `[[`, "std.err")))
      pval_DC <- unname(unlist(lapply(res_list, `[[`, "pVal")))
      i2_DC   <- unname(unlist(lapply(res_list, `[[`, "I2.stat")))
      idx <- match(var_names, grouped_vars)
      stopifnot(!anyNA(idx))
      return(list(
        estimate = est_DC[idx],
        std.err  = se_DC[idx],
        pVal     = pval_DC[idx],
        I2.stat  = i2_DC[idx]
      ))
    }else {
      if (N > 250) {
        warning("Large number of outcomes (N = ", N, "). Computation may be slow; consider using the divide-and-conquer approach (set useDivideConquer = TRUE).", call. = FALSE)
      }
      return(.MetaHD_core(y = y, Slist = Slist, Psi = Psi, method = method, bscov = bscov, est.wscor = est.wscor, shrinkCor = shrinkCor, optim.algorithm = optim.algorithm, optim.maxiter = optim.maxiter, rigls.iter = rigls.iter, initPsi = initPsi, impute.na = impute.na, impute.var = impute.var, N = N, K = K, p = p, q = q))
    }
  }else if (method == "REM" || method == "FEM") {
    WSVar <- dimCheck(Slist, K, N)
    estimate <- std.err <- pvalue <- I2 <- numeric(N)
    for (n in 1:N) {
      y_n <- y[, n]
      WSVar_n <- WSVar[, n]
      Psi_n <- if (is.null(Psi)) NULL else Psi[n, n, drop = FALSE]
      initPsi_n <- if (is.null(initPsi)) NULL else initPsi[n, n, drop = FALSE]
      res <- .MetaHD_core(y = y_n, Slist = WSVar_n, Psi = Psi_n, method = method, bscov = bscov, est.wscor = est.wscor, shrinkCor = shrinkCor, optim.algorithm = optim.algorithm, optim.maxiter = optim.maxiter, rigls.iter = rigls.iter, initPsi = initPsi_n, impute.na = impute.na, impute.var = impute.var, N = 1, K = K, p = p, q = q)
      estimate[n] <- res$estimate
      std.err[n] <- res$std.err
      pvalue[n] <- res$pVal
      I2[n] <- res$I2.stat
    }
    return(list(estimate = estimate,
                std.err = std.err,
                pVal = pvalue,
                I2.stat = I2))
  }
}

