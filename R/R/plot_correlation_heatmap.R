#' Plot a heatmap of correlations among outcomes
#'
#' Visualises correlations among outcomes from a single covariance or
#' correlation matrix. A covariance matrix is converted to correlations
#' internally; a correlation matrix can be supplied directly via
#' \code{is.corr = TRUE} (for example, the within-study or between-study
#' correlations among outcomes). This is useful for understanding the
#' correlation structure of the outcomes and assessing whether the multivariate
#' meta-analysis approach is likely to offer advantages over univariate
#' methods. This function uses the \pkg{ComplexHeatmap} package (Gu, 2022) for
#' rendering.
#'
#' @param mat An N x N covariance or correlation matrix, where N is the
#'   number of outcomes, such as one of the within-study covariance matrices
#'   returned by \code{MetaHDInput()} or a correlation matrix computed by the
#'   user. By default it is treated as a covariance matrix and converted to
#'   correlations internally; set \code{is.corr = TRUE} to supply a correlation
#'   matrix directly. 
#' @param is.corr Logical. If \code{FALSE} (default), \code{mat} is treated as
#'   a covariance matrix and converted to correlations with
#'   \code{\link[stats]{cov2cor}}. If \code{TRUE}, it is used directly as a
#'   correlation matrix.
#' @param outcome_names An optional character vector of outcome names for
#'   axis labels. If \code{NULL} (default), row and column names of
#'   \code{mat} are used.
#' @param title Character string for the heatmap title. If \code{NULL}
#'   (default), a default title is used. Use \code{title = ""} to omit the
#'   title.
#' @param col_low Character string specifying the colour for correlation
#'   of -1. Default is \code{"blue"}.
#' @param col_mid Character string specifying the colour for correlation
#'   of 0. Default is \code{"white"}.
#' @param col_high Character string specifying the colour for correlation
#'   of 1. Default is \code{"red"}.
#' @param cluster_rows Logical. Whether to cluster rows. Default is
#'   \code{TRUE}.
#' @param cluster_columns Logical. Whether to cluster columns. Default is
#'   \code{TRUE}.
#' @param row_km Integer >= 1. Number of groups into which outcomes are
#'   partitioned by k-means clustering, drawn as separate heatmap slices,
#'   showing groups of correlated outcomes (Gu, 2022). Because the correlation 
#'   matrix is symmetric, the same partition is applied to both rows and columns 
#'   so the diagonal blocks stay aligned. Default \code{1} (no partitioning).
#' @param show_row_names,show_column_names Logical or \code{NULL}. Whether to
#'   display row / column (outcome) labels. If \code{NULL} (default), labels
#'   are shown when there are at most 50 outcomes and hidden otherwise (with a
#'   message). Set \code{TRUE} or \code{FALSE} to override.
#' @param names_fontsize Numeric. Font size for row and column labels.
#'   Default is \code{7}.
#' @param legend_title Character string for the colour legend title.
#'   Default is \code{"Correlation"}.
#' @param file Optional file path to save the plot. If \code{NULL} (default),
#'   the plot is drawn to the active graphics device. Supported extensions:
#'   \code{.pdf}, \code{.svg}, \code{.png}, \code{.jpeg}/\code{.jpg},
#'   \code{.tiff}/\code{.tif}, \code{.bmp}.
#' @param width,height Plot dimensions. \emph{Used only when \code{file} is
#'   set.} Default 8 x 6.
#' @param units Units for \code{width} and \code{height}: \code{"in"}
#'   (default), \code{"cm"}, or \code{"mm"}. \emph{Used only when
#'   \code{file} is set.}
#' @param dpi Resolution for raster formats (PNG, JPEG, TIFF, BMP). Ignored
#'   for vector formats (PDF, SVG). Default 300. \emph{Used only when
#'   \code{file} is set.}
#' @param ... Further arguments passed to \code{ComplexHeatmap::Heatmap()}.
#'
#' @details
#' Unless \code{is.corr = TRUE}, \code{mat} is treated as a covariance matrix
#' and converted to a correlation matrix internally using \code{cov2cor()}
#' before plotting. The colour scale runs from \code{col_low} at -1 through
#' \code{col_mid} at 0 to \code{col_high} at 1.
#'
#' @section Large numbers of outcomes:
#' With many outcomes the row and column labels become unreadable and individual
#' cells shrink to sub-pixel size. To keep the plot readable, axis labels are
#' hidden automatically beyond 50 outcomes (override via \code{show_row_names} /
#' \code{show_column_names}). Very large heatmaps are rasterised automatically by
#' \pkg{ComplexHeatmap} to keep file size and rendering manageable.
#'
#' @return Invisibly returns the \code{Heatmap} object.
#'
#' @seealso
#' \code{\link{plot_effect_heatmap}} for comparing pooled effect sizes
#' across methods.
#' \code{\link{plot.MetaHDResult}} for Venn diagrams, UpSet plots, and ROC
#' curves.
#'
#' @references
#' Gu, Z. (2022).
#' \emph{Complex heatmap visualization}.
#' iMeta, 1, e43.
#' \doi{10.1002/imt2.43}
#'
#' @examples
#' Y <- simdata.1$Y
#' Slist <- simdata.1$Slist
#'
#' # Example 1: a within-study covariance matrix (converted to correlations)
#' plot_correlation_heatmap(Slist[[1]])
#'
#' # Example 2: partition into k-means blocks
#' set.seed(123)  # reproducible k-means blocks
#' plot_correlation_heatmap(Slist[[1]], row_km = 4)
#'
#' # Example 3: supply a correlation matrix directly
#' cormat <- estimateCorMat(Y)
#' plot_correlation_heatmap(cormat, is.corr = TRUE)
#'
#' \dontrun{
#' # Example 4: save to PNG
#' plot_correlation_heatmap(
#'   Slist[[1]],
#'   file = "correlation_heatmap.png",
#'   width = 10,
#'   height = 8,
#'   dpi = 600
#' )
#' }
#'
#' @importFrom grDevices dev.off
#' @importFrom ComplexHeatmap Heatmap draw
#' @importFrom circlize colorRamp2
#' @importFrom stats cov2cor
#' @export
plot_correlation_heatmap <- function(mat,
                                     is.corr = FALSE,
                                     outcome_names = NULL,
                                     title = NULL,
                                     col_low = "blue",
                                     col_mid = "white",
                                     col_high = "red",
                                     cluster_rows = TRUE,
                                     cluster_columns = TRUE,
                                     row_km = 1,
                                     show_row_names = NULL,
                                     show_column_names = NULL,
                                     names_fontsize = 7,
                                     legend_title = "Correlation",
                                     file = NULL,
                                     width = 8,
                                     height = 7,
                                     units = c("in", "cm", "mm"),
                                     dpi = 300,
                                     ...) {
  if (!is.matrix(mat)) {
    stop("`mat` must be a covariance or correlation matrix.", call. = FALSE)
  }
  if (nrow(mat) != ncol(mat)) {
    stop("`mat` must be a square matrix.", call. = FALSE)
  }
  units <- match.arg(units)
  if (!is.null(file)) {
    .open_plot_device(file, width = width, height = height,
                      units = units, dpi = dpi)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  cormat <- if (isTRUE(is.corr)) mat else cov2cor(mat)
  if (!is.null(outcome_names)) {
    if (length(outcome_names) != nrow(cormat)) {
      stop(
        "`outcome_names` length (", length(outcome_names), ") ",
        "must match number of outcomes (", nrow(cormat), ").",
        call. = FALSE
      )
    }
    rownames(cormat) <- colnames(cormat) <- outcome_names
  }
  if (is.null(title)) title <- "Heatmap of Correlations"
  if (!is.numeric(row_km) || length(row_km) != 1L || row_km < 1) {
    stop("`row_km` must be a single integer >= 1.", call. = FALSE)
  }
  row_km <- as.integer(row_km)
  n_out  <- nrow(cormat)
  split <- NULL
  if (row_km > 1L) {
    km <- stats::kmeans(cormat, centers = min(row_km, n_out))
    split <- factor(km$cluster)
  }
  hide_msg <- FALSE
  if (is.null(show_row_names)) {
    show_row_names <- n_out <= 50L
    if (!show_row_names) hide_msg <- TRUE
  }
  if (is.null(show_column_names)) {
    show_column_names <- n_out <= 50L
    if (!show_column_names) hide_msg <- TRUE
  }
  if (hide_msg) {
    message(
      "Hiding axis labels for ", n_out, " outcomes (> 50). Set ",
      "`show_row_names`/`show_column_names = TRUE` to force labels."
    )
  }
  col_fun <- circlize::colorRamp2(
    c(-1, 0, 1),
    c(col_low, col_mid, col_high)
  )
  ht <- ComplexHeatmap::Heatmap(
    cormat,
    name = legend_title,
    col = col_fun,
    column_title = title,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    row_split = split,
    column_split = split,
    show_row_names = show_row_names,
    show_column_names = show_column_names,
    row_names_gp = grid::gpar(fontsize = names_fontsize),
    column_names_gp = grid::gpar(fontsize = names_fontsize),
    heatmap_legend_param = list(title = legend_title),
    ...
  )
  ComplexHeatmap::draw(ht)
  invisible(ht)
}
