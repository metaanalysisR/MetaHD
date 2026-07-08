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
#'   DCgroups = NULL,
#'   parallel = FALSE,
#'   dendro.method = c("dynamicTreeCut", "fixedHeight", "fixedK", "optimalK"),
#'   dendro.height = NULL,
#'   dendro.k = NULL,
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
#' @param Slist A K-dimensional list of N x N matrices representing within-study variances and covariances of the treatment effects. If within-study correlations are not available, provide the associated variances of the treatment effects as a K x N matrix and set est.wscor = TRUE. For method = "REM" or method = "FEM", provide the associated variances of the treatment effects as a K x N matrix.
#' @param Psi N x N matrix representing between-study variances and covariances of the treatment effects. (optional, if not specified this will be estimated internally by "MetaHD" using "estimateBSvar" and "estimateCorMat" functions in "MetaHD" package).
#' @param method estimation method: "multi" for multivarite meta-analysis model fitted through restricted maximum likelihood estimation where the between-study covariance structure can be selected via 'bscov', "REM" for univariate random-effects model fitted through restricted maximum likelihood estimation and "FEM" for univariate fixed-effects model.
#' @param bscov a character vector defining the structure of the random-effects covariance matrix. Among available covariance structures, the user can select "unstructured" to obtain between-study covariance matrix with diagonal elements (variances) estimated using restricted maximum likelihood and off-diagonal elements (co-variances) reflecting the correlations estimated via shrinkage, "diag" (diagonal) for between-study variances as diagonal elements and zero co-variances, and "none" for zero between-study variances and co-variances.
#' @param useDivideConquer a logical value indicating whether to use the divide-and-conquer implementation of the fastMetaHD model. This option is used only when method = "multi". Default is FALSE.
#' @param DCgroups A list of outcome groups for the divide-and-conquer approach. Each element should be a character vector containing the names of outcomes belonging to the same cluster. Outcome names must match column names in the input data (Y), and each outcome may appear in at most one group. If `NULL`, groups are determined automatically within the function using hierarchical clustering followed by a dynamic tree cut.
#' @param parallel a logical value indicating whether to enable parallel computation for the divide-and-conquer approach. Default is \code{FALSE}. See also Details.
#' @param dendro.method Character string specifying the method used to cut the dendrogram when \code{useDivideConquer = TRUE} and \code{DCgroups = NULL}.
#'   One of:
#'   \describe{
#'     \item{\code{"dynamicTreeCut"}}{Default. Automatically determines the
#'       number of clusters using the dynamic tree cut algorithm.}
#'     \item{\code{"fixedHeight"}}{Cuts the dendrogram at a fixed height
#'       specified via \code{dendro.height}.}
#'     \item{\code{"fixedK"}}{Cuts the dendrogram into a fixed number of
#'       clusters specified via \code{dendro.k}.}
#'     \item{\code{"optimalK"}}{Automatically selects the optimal number of
#'       clusters by maximising the average silhouette width across all
#'       candidate partitions.}
#'   }
#' @param dendro.height A numeric value specifying the height at which to cut the dendrogram. Only used when \code{dendro.method = "fixedHeight"}. Must be a positive number within the range of the dendrogram heights. Use \code{plot_dendrogram(Y)} to visualize the dendrogram and choose an appropriate height before running \code{MetaHD}. If \code{NULL} (default) and \code{dendro.method = "fixedHeight"}, an error is thrown.
#' @param dendro.k A positive integer specifying the desired number of clusters. Only used when \code{dendro.method = "fixedK"}. Must be between 2 and \code{ncol(Y) - 1}. If \code{NULL} (default) and \code{dendro.method = "fixedK"}, an error is thrown.
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
#'   \item \code{clustering_diagnostics}: A list of clustering diagnostic outputs
#'     returned when \code{useDivideConquer = TRUE} and \code{DCgroups = NULL}.
#'     Contains:
#'     \itemize{
#'       \item \code{avg_silhouette}: Overall average silhouette width across
#'         all outcomes.
#'       \item \code{cluster_avg}: Average silhouette width per cluster.
#'       \item \code{n_clusters}: Number of clusters identified.
#'       \item \code{cluster_sizes}: Number of outcomes in each cluster.
#'       \item \code{silhouette_object}: The full silhouette object of class
#'         \code{"silhouette"}, which can be passed directly to \code{plot()}
#'         for visualization.
#'       \item \code{silhouette_by_k}: A data frame with columns \code{k} and
#'         \code{avg_silhouette}, returned only when
#'         \code{dendro.method = "optimalK"}, showing the average silhouette
#'         width at each candidate partition used to select the optimal
#'         number of clusters.
#'     }
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
#' Cambridge NA Report NA2009/06, University of Cambridge, 26, 26--46.
#' @references
#' Sera F, Armstrong B, Blangiardo M, et al. (2019).
#' \emph{An extended mixed-effects framework for meta-analysis}.
#' Statistics in Medicine, 38, 5429--5444.
#' @references
#' Schaefer J, Strimmer K (2005).
#' \emph{A shrinkage approach to large-scale covariance estimation and implications for functional genomics}.
#' Statistical Applications in Genetics and Molecular Biology, 4, 32.
#' @references 
#' Langfelder, P., Zhang, B. and Horvath, S. (2008).
#' \emph{Defining clusters from a hierarchical cluster tree: the Dynamic Tree Cut package for R}.
#' Bioinformatics, 24(5), pp. 719--720.
#' @references
#' Rousseeuw, P.J. (1987). 
#' \emph{Silhouettes: a graphical aid to the interpretation and validation of cluster analysis}.
#' Journal of Computational and Applied Mathematics, 20, pp. 53--65.
#' 
#' @export MetaHD
#'
#' @name MetaHD
#' @useDynLib MetaHD, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom matrixcalc is.positive.definite
#' @importFrom Matrix nearPD
#' @importFrom stats na.omit cor optim pnorm as.dist hclust cutree
#' @importFrom utils capture.output
#' @importFrom dynamicTreeCut cutreeDynamic
#' @importFrom cluster silhouette
#' @importFrom nloptr bobyqa
#' @importFrom corpcor cor.shrink
#' @importFrom future.apply future_lapply
NULL

sourceCpp("src/cpp_XtVX.cpp")

MetaHD <- function(Y, Slist, Psi = NULL, method = c("multi","REM","FEM"),
                   bscov = c("unstructured", "diag", "none"),
                   useDivideConquer = FALSE, DCgroups = NULL, parallel = FALSE,
                   dendro.method = c("dynamicTreeCut", "fixedHeight", "fixedK", "optimalK"),
                   dendro.height = NULL, dendro.k = NULL,
                   est.wscor = FALSE, shrinkCor = TRUE, impute.na = FALSE,
                   optim.algorithm = c("BOBYQA","hybrid","L-BFGS-B"), 
                   optim.maxiter = 2000, rigls.iter = 1, initPsi = NULL, impute.var = 10^4){
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
  study_names <- rownames(y)
  var_names <- colnames(y)
  slist_colnames <- if (is.list(Slist)) colnames(Slist[[1]]) else colnames(Slist)
  if (is.null(study_names)) {
    study_names <- paste0("Study", seq_len(K))
    rownames(y) <- study_names
  }
  if (is.null(var_names) && is.null(slist_colnames)) {
    var_names <- paste0("Var", seq_len(N))
    colnames(y) <- var_names
  } else if (is.null(var_names) && !is.null(slist_colnames)) {
    colnames(y) <- slist_colnames
    var_names <- slist_colnames
  } else if (!is.null(var_names) && is.null(slist_colnames)) {
    if (is.list(Slist)) {
      Slist <- lapply(Slist, function(Sk) {
        colnames(Sk) <- var_names
        rownames(Sk) <- var_names
        Sk
      })
    } else {
      colnames(Slist) <- var_names
      rownames(Slist) <- study_names
    }
  } else {
    if (!identical(var_names, slist_colnames)) {
      stop(
        "Column names of Y and Slist do not match.\n",
        "Please ensure Y and Slist refer to the same outcomes in the same order."
      )
    }
  }
  if (is.list(Slist)) {
    Slist <- lapply(Slist, function(Sk) {
      if (is.null(rownames(Sk))) rownames(Sk) <- var_names
      if (is.null(colnames(Sk))) colnames(Sk) <- var_names
      Sk
    })
  } else {
    if (is.null(rownames(Slist))) rownames(Slist) <- study_names
    if (is.null(colnames(Slist))) colnames(Slist) <- var_names
  }
  if (!is.null(Psi) && (is.null(colnames(Psi)) || is.null(rownames(Psi)))) {
    colnames(Psi) <- rownames(Psi) <- var_names
  }
  if (!is.null(initPsi) && (is.null(colnames(initPsi)) || is.null(rownames(initPsi)))) {
    colnames(initPsi) <- rownames(initPsi) <- var_names
  }
  if (method == "multi") {
    if (useDivideConquer) {
      if (is.null(DCgroups)) {
        cormat <- estimateCorMat(y, shrinkCor, impute.na)
        distMat <- as.dist(1 - abs(cormat))
        hc <- hclust(distMat, method = "complete")
        cut_result <- cut_dendrogram(hc, distMat,
                                     method = dendro.method,
                                     height = dendro.height,
                                     k = dendro.k)
        clusters <- cut_result$clusters
        silhouette_by_k <- cut_result$silhouette_by_k
        names(clusters) <- hc$labels
        groups <- split(names(clusters), clusters)
        clust_diagnostics <- compute_clustering_diagnostics(clusters, distMat)
        clust_diagnostics$silhouette_by_k <- silhouette_by_k
        if (is.na(clust_diagnostics$avg_silhouette)) {
          message(
            "Clustering diagnostics: only 1 cluster found ",
            "(N = ", sum(clust_diagnostics$cluster_sizes), " outcomes). ",
            "Silhouette scores require at least 2 clusters. "
          )
        } else {
          message(
            "Clustering diagnostics: average silhouette width = ",
            round(clust_diagnostics$avg_silhouette, 3),
            " (", clust_diagnostics$n_clusters, " clusters, ",
            "N = ", sum(clust_diagnostics$cluster_sizes), " outcomes)."
          )
        }
      } else {
        if (!is.list(DCgroups)) {
          stop("`DCgroups` must be a list of character vectors.")
        }
        if (!all(vapply(DCgroups, is.character, logical(1)))) {
          stop("Each element of `DCgroups` must be a character vector of outcome names.")
        }
        if (any(lengths(DCgroups) == 0)) {
          stop("`DCgroups` contains empty groups.")
        }
        all_outcomes <- unlist(DCgroups, use.names = FALSE)
        if (any(duplicated(all_outcomes))) {
          stop("Each outcome must appear in at most one group in `DCgroups`.")
        }
        if (!all(all_outcomes %in% colnames(y))) {
          missing <- setdiff(all_outcomes, colnames(y))
          stop("The following outcomes in `DCgroups` are not found in `Y`: ",
               paste(missing, collapse = ", "))
        }
        groups <- DCgroups
        clust_diagnostics <- NULL   
      }
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
      est_DC <- unname(unlist(lapply(res_list, `[[`, "estimate")))
      se_DC <- unname(unlist(lapply(res_list, `[[`, "std.err")))
      pval_DC <- unname(unlist(lapply(res_list, `[[`, "pVal")))
      i2_DC <- unname(unlist(lapply(res_list, function(r) {
        v <- r$I2.stat
        if (length(v) == length(r$estimate) + 1L) v[-1L] else v
      })))
      idx <- match(var_names, grouped_vars)
      stopifnot(!anyNA(idx))
      out <- list(
        estimate = est_DC[idx],
        std.err = se_DC[idx],
        pVal = pval_DC[idx],
        I2.stat = i2_DC[idx],
        clustering_diagnostics = if (useDivideConquer && is.null(DCgroups)) clust_diagnostics else NULL
      )
      class(out) <- "MetaHD"
      return(out)
    }else {
      if (N > 250) {
        warning("Large number of outcomes (N = ", N, "). Computation may be slow; consider using the divide-and-conquer approach (set useDivideConquer = TRUE).", call. = FALSE)
      }
      out <- .MetaHD_core(y = y, Slist = Slist, Psi = Psi, method = method, bscov = bscov, est.wscor = est.wscor, shrinkCor = shrinkCor, optim.algorithm = optim.algorithm, optim.maxiter = optim.maxiter, rigls.iter = rigls.iter, initPsi = initPsi, impute.na = impute.na, impute.var = impute.var, N = N, K = K, p = p, q = q)
      class(out) <- "MetaHD"
      return(out)
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
    out <- list(estimate = estimate,
                std.err = std.err,
                pVal = pvalue,
                I2.stat = I2)
    class(out) <- "MetaHD"
    return(out)
  }
}

#' Summarise or print MetaHD results
#'
#' \code{summary()} returns the combined meta-analysis results from a
#' \code{\link{MetaHD}} fit as a tidy data frame, one row per outcome. It is
#' convenient for inspecting the results and for exporting them to a
#' \code{.csv} file with \code{\link[utils]{write.csv}}. \code{print()} shows
#' a concise overview of the fit.
#'
#' @param object,x A \code{"MetaHD"} object returned by \code{\link{MetaHD}}.
#' @param outcome_names Optional character vector of outcome labels, one per
#'   outcome. If \code{NULL} (default), the names of \code{x$estimate} are used
#'   when available, otherwise \code{Var1, Var2, ...}.
#' @param ... Additional arguments. For \code{print()} these are forwarded to
#'   \code{summary()} (e.g. \code{outcome_names}); \code{summary()} itself
#'   ignores them.
#'
#' @return \code{summary()} returns a data frame with one row per outcome and
#'   the columns \code{outcome}, \code{estimate}, \code{std.err}, \code{pVal}
#'   and \code{I2.stat}. \code{print()} displays a concise overview and
#'   invisibly returns \code{x}.
#'
#' @examples
#' Y <- simdata.1$Y
#' Slist <- simdata.1$Slist
#' model <- MetaHD(Y, Slist)
#'
#' model # concise overview
#' results <- summary(model, outcome_names = colnames(Y))
#' head(results)
#' 
#' \dontrun{
#' # Export the results table to CSV
#' write.csv(summary(model, outcome_names = colnames(Y)),
#'           "metahd_results.csv", row.names = FALSE)
#' }
#'
#' @name MetaHD-outputs
#' @seealso \code{\link{MetaHD}}
NULL

#' @rdname MetaHD-outputs
#' @method summary MetaHD
#' @export
summary.MetaHD <- function(object, outcome_names = NULL, ...) {
  n <- length(object$estimate)
  if (is.null(outcome_names)) {
    outcome_names <- if (!is.null(names(object$estimate))) {
      names(object$estimate)
    } else {
      paste0("Var", seq_len(n))
    }
  }
  if (length(outcome_names) != n) {
    stop(
      "`outcome_names` length (", length(outcome_names), ") ",
      "must match number of outcomes (", n, ").",
      call. = FALSE
    )
  }
  i2 <- object$I2.stat
  if (length(i2) == n + 1L) {
    i2 <- i2[-1L]
  } else if (length(i2) != n) {
    i2 <- rep(NA_real_, n)
  }
  data.frame(
    outcome = outcome_names,
    estimate = object$estimate,
    std.err = object$std.err,
    pVal = object$pVal,
    I2.stat = i2,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

#' @rdname MetaHD-outputs
#' @method print MetaHD
#' @export
print.MetaHD <- function(x, ...) {
  n <- length(x$estimate)
  cat("MetaHD multivariate meta-analysis\n")
  cat("Outcomes:", n, "\n\n")
  tab <- summary(x, ...)
  print(utils::head(tab, 10L), row.names = FALSE)
  if (n > 10L) {
    cat("... and ", n - 10L,
        " more outcome(s); use summary() for the full table.\n", sep = "")
  }
  invisible(x)
}

