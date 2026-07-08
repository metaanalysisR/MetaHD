#' Create a `MetaHD` result object
#'
#' Constructs an object containing combined results from different meta-analysis
#' methods, for downstream visualisation via Venn diagrams (Tav, 2025), UpSet
#' plots (Conway et al., 2017), and ROC curves (Robin et al., 2011). The
#' \code{sets} argument accepts numeric vectors of p-values (thresholded by
#' \code{alpha}) or effect sizes (ranked and selected by \code{top_n}). Note
#' that ROC curves require numeric p-values and a \code{truth} vector.
#'
#' @param sets A named list where each element is a numeric vector of p-values
#'   or effect sizes for one method, one value per outcome. All vectors must
#'   have the same length. Selection of outcomes requires either \code{alpha}
#'   (for p-values; available for Venn, UpSet, and ROC plots) or \code{top_n}
#'   (for effect sizes; available for Venn and UpSet plots).
#'   Names of the list are used as method labels in the plots.
#' @param outcome_names An optional character vector of outcome names matching
#'   the length of the numeric vectors in \code{sets}. If \code{NULL} (default),
#'   outcomes are indexed by position.
#' @param truth An optional logical or numeric (0/1) vector of length N
#'   indicating true signals. Required for ROC curve visualisation via
#'   \code{plot(res, type = "ROC")}.
#' @param alpha A numeric value specifying the significance threshold used
#'   when \code{sets} contains p-value vectors. Default is \code{0.05}.
#'   Ignored when \code{top_n} is specified.
#' @param top_n Optional positive integer. If specified, the \code{top_n}
#'   outcomes with the largest absolute values are selected from each numeric
#'   vector in \code{sets}. Useful when \code{sets} contains effect sizes
#'   rather than p-values. Overrides \code{alpha} when specified.
#'   Note: ROC curves are not available when \code{top_n} is used.
#'
#' @return An object of class \code{"MetaHDResult"} containing:
#' \describe{
#'   \item{\code{sig_df}}{A data frame of binary indicators
#'     (1 = significant/selected, 0 = otherwise).}
#'   \item{\code{sets}}{The original input list.}
#'   \item{\code{truth}}{The provided ground truth vector (if any).}
#'   \item{\code{alpha}}{The significance threshold used (if applicable).}
#'   \item{\code{top_n}}{The top-N value used (if applicable).}
#' }
#'
#' @details
#' This object is designed to be used with the S3 method
#' \code{\link{plot.MetaHDResult}} for visualisation of \pkg{MetaHD} results
#' via Venn diagrams (Tav, 2025), UpSet plots (Conway et al., 2017), and 
#' ROC curves (Robin et al., 2011).
#'
#' The function supports two input modes:
#' \enumerate{
#'   \item \strong{P-value thresholding}: Supply numeric p-value vectors
#'     with \code{alpha} (and optionally \code{outcome_names}). Available for
#'     Venn, UpSet, and ROC plots.
#'   \item \strong{Top-N selection}: Supply numeric vectors (e.g. effect
#'     sizes) with \code{top_n} (and optionally \code{outcome_names}).
#'     Available for Venn and UpSet plots only.
#' }
#'
#' @seealso
#' \code{\link{plot.MetaHDResult}} for Venn diagrams, UpSet plots, and ROC
#' curves based on this object.
#' \code{\link{plot_effect_heatmap}} for comparing pooled effect size
#' estimates across methods.
#' \code{\link{plot_correlation_heatmap}} for visualising correlations
#' among outcomes.
#'
#' @references
#' Tav, C. (2025).
#' \emph{gVenn: Proportional Venn and UpSet Diagrams for Gene Sets and
#' Genomic Regions}.
#' \doi{10.18129/B9.bioc.gVenn}.
#' R package version 1.1.1,
#' \url{https://bioconductor.org/packages/gVenn}
#'
#' @references
#' Conway, J.R., Lex, A., and Gehlenborg, N. (2017).
#' \emph{UpSetR: an R package for the visualization of intersecting sets
#' and their properties}.
#' Bioinformatics, 33(18), 2938--2940.
#' \doi{10.1093/bioinformatics/btx364}
#'
#' @references
#' Robin, X., Turck, N., Hainard, A., Tiberti, N., Lisacek, F.,
#' Sanchez, J., and Mueller, M. (2011).
#' \emph{pROC: an open-source package for R and S+ to analyze and
#' compare ROC curves}.
#' BMC Bioinformatics, 12, 77.
#' \doi{10.1186/1471-2105-12-77}
#'
#' @examples
#' set.seed(123)
#' N <- 100
#' truth <- rbinom(N, 1, 0.2)   # 20% of features are true signals
#'
#' # Example 1: p-value threshold (default alpha = 0.05)
#' # Supports venn, upset, and ROC plots
#' res <- MetaHDResult(
#'   sets = list(
#'     method_A = runif(N)^ifelse(truth, 5, 1),
#'     method_B = runif(N)^ifelse(truth, 3, 1),
#'     method_C = runif(N)^ifelse(truth, 2, 1)
#'   ),
#'   truth = truth
#' )
#'
#' # Example 2: top-N selection by effect size
#' # Supports venn and upset plots only
#' res2 <- MetaHDResult(
#'   sets = list(
#'     method_A = rnorm(N),
#'     method_B = rnorm(N),
#'     method_C = rnorm(N)
#'   ),
#'   outcome_names = paste0("Var", seq_len(N)),
#'   top_n = 10
#' )
#'
#' @export
MetaHDResult <- function(sets,
                         outcome_names = NULL,
                         truth = NULL,
                         alpha = 0.05,
                         top_n = NULL) {
  if (!is.list(sets) || length(sets) == 0L) {
    stop("`sets` must be a non-empty named list.", call. = FALSE)
  }
  if (is.null(names(sets)) || any(!nzchar(names(sets)))) {
    stop("`sets` must be a named list (one name per method).", call. = FALSE)
  }
  if (anyDuplicated(names(sets))) {
    stop("`sets` names must be unique.", call. = FALSE)
  }
  if (!all(vapply(sets, is.numeric, logical(1)))) {
    stop(
      "`sets` must contain numeric vectors (p-values or effect sizes).",
      call. = FALSE
    )
  }
  lens <- vapply(sets, length, integer(1))
  if (length(unique(lens)) != 1L) {
    stop("All numeric vectors in `sets` must have the same length.",
         call. = FALSE)
  }
  if (is.null(outcome_names)) {
    outcome_names <- as.character(seq_len(lens[[1L]]))
  }
  if (length(outcome_names) != lens[[1L]]) {
    stop(
      "`outcome_names` length (", length(outcome_names), ") ",
      "must match length of numeric vectors (", lens[[1L]], ").",
      call. = FALSE
    )
  }
  if (!is.null(top_n)) {
    if (!is.numeric(top_n) || length(top_n) != 1L || top_n < 1) {
      stop("`top_n` must be a single positive integer.", call. = FALSE)
    }
    top_n <- as.integer(top_n)
    sig_df <- data.frame(
      lapply(sets, function(x) {
        idx <- order(abs(x), decreasing = TRUE)[seq_len(min(top_n, length(x)))]
        result <- integer(length(x))
        result[idx] <- 1L
        result
      }),
      check.names = FALSE
    )
    rownames(sig_df) <- outcome_names
  } else {
    if (!is.numeric(alpha) || length(alpha) != 1L ||
        alpha <= 0 || alpha >= 1) {
      stop("`alpha` must be a numeric value in (0, 1).", call. = FALSE)
    }
    sig_df <- data.frame(
      lapply(sets, function(p) as.integer(p < alpha)),
      check.names = FALSE
    )
    rownames(sig_df) <- outcome_names
  }
  if (!is.null(truth)) {
    if (!is.numeric(truth) && !is.logical(truth)) {
      stop("`truth` must be a logical or numeric (0/1) vector.", call. = FALSE)
    }
    if (length(truth) != nrow(sig_df)) {
      stop("`truth` must have the same length as the outcome vectors.",
           call. = FALSE)
    }
    if (is.numeric(truth) && !all(truth %in% c(0, 1))) {
      stop("`truth` must contain only 0 and 1 if numeric.", call. = FALSE)
    }
    truth <- as.integer(as.logical(truth))
  }
  result <- list(
    sig_df = sig_df,
    sets = sets,
    truth = truth,
    alpha = alpha,
    top_n = top_n
  )
  class(result) <- "MetaHDResult"
  return(result)
}
