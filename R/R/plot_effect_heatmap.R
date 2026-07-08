#' Plot a heatmap of pooled effect sizes across meta-analysis methods
#'
#' Produces a heatmap with outcomes on the y-axis and meta-analysis methods
#' on the x-axis, allowing visual comparison of pooled effect sizes across
#' methods. Rows (outcomes) are clustered by default to reveal patterns of
#' agreement and disagreement across methods. This function uses the
#' \pkg{ComplexHeatmap} package (Gu, 2022) for rendering.
#'
#' @param estimates A named list of numeric vectors, where each element
#'   corresponds to a meta-analysis method and contains the pooled effect
#'   size estimates for all outcomes. All vectors must have the same length.
#'   Names of the list are used as column labels (methods). 
#' @param outcome_names An optional character vector of outcome names to use
#'   as row labels. If \code{NULL} (default), outcomes are labelled
#'   \code{Var1, Var2, ...}.
#' @param title Character string for the heatmap title. Default is
#'   \code{"Pooled Effect Sizes Across Methods"}. Use \code{title = ""} to
#'   omit the title.
#' @param col_low Character string specifying the colour for the lowest
#'   values. Default is \code{"blue"}.
#' @param col_mid Character string specifying the colour for the midpoint
#'   (zero). Default is \code{"white"}.
#' @param col_high Character string specifying the colour for the highest
#'   values. Default is \code{"red"}.
#' @param cluster_rows Logical. Whether to cluster rows (outcomes). Default
#'   is \code{TRUE}.
#' @param cluster_columns Logical. Whether to cluster columns (methods).
#'   Default is \code{FALSE}.
#' @param row_km Integer >= 1. Number of groups into which outcomes (rows)
#'   are partitioned by k-means clustering, each drawn as a separate heatmap
#'   slice, showing subgroups of outcomes with similar effect patterns
#'   (Gu, 2022). Default \code{1} (no partitioning).
#' @param show_row_names Logical or \code{NULL}. Whether to display row
#'   (outcome) labels. If \code{NULL} (default), labels are shown when there
#'   are at most 50 outcomes (and no \code{label_top_n} is given) and hidden
#'   otherwise (with a message). Set \code{TRUE} or \code{FALSE} to override.
#' @param label_top_n Optional positive integer. If set, the \code{label_top_n}
#'   outcomes with the largest absolute pooled effect (across methods) are
#'   labelled with connector lines via a mark annotation (Gu, 2022). Useful for
#'   highlighting key outcomes when row labels are hidden. When labels are
#'   marked this way, the full row labels are hidden by default.
#' @param row_names_fontsize Numeric. Font size for row (outcome) labels.
#'   Default is \code{8}.
#' @param column_names_fontsize Numeric. Font size for column (method)
#'   labels. Default is \code{10}.
#' @param legend_title Character string for the colour legend title. Default
#'   is \code{"Pooled\nEffect Size"}.
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
#' The colour scale is centred at zero by default, with \code{col_low}
#' indicating negative effect sizes and \code{col_high} indicating positive
#' effect sizes. This allows easy visual identification of outcomes with
#' consistent direction across methods as well as those where methods disagree.
#'
#' @section Large numbers of outcomes:
#' With many outcomes per-row labels become unreadable and individual rows are
#' squeezed to sub-pixel height. To keep the plot readable, row labels are hidden
#' automatically beyond 50 outcomes (override via \code{show_row_names}). Alternatively,
#' label only the key outcomes with \code{label_top_n}, which marks the outcomes 
#' with the largest absolute effect with connector lines while leaving the rest unlabelled. 
#' Very large heatmaps are rasterised automatically by \pkg{ComplexHeatmap} 
#' to keep file size and rendering manageable.
#'
#' @return Invisibly returns the \code{Heatmap} object.
#'
#' @seealso
#' \code{\link{plot_correlation_heatmap}} for visualising correlations
#' among outcomes.
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
#' K <- nrow(Y)
#' N <- ncol(Y)
#' Smat <- matrix(0, nrow = K, ncol = N,
#'                dimnames = list(rownames(Y), colnames(Y)))
#' for (i in 1:K) Smat[i, ] <- diag(Slist[[i]])
#'
#' model_multi <- MetaHD(Y = Y, Slist = Slist, method = "multi")
#' model_rem <- MetaHD(Y = Y, Slist = Smat,  method = "REM")
#' model_fem <- MetaHD(Y = Y, Slist = Smat,  method = "FEM")
#'
#' # Example 1: basic heatmap with custom colours
#' plot_effect_heatmap(
#'   estimates = list(
#'     fastMetaHD = model_multi$estimate,
#'     REM = model_rem$estimate,
#'     FEM = model_fem$estimate
#'   ),
#'   outcome_names = colnames(Y),
#'   col_low = "purple",
#'   col_high = "orange"
#' )
#'
#' # Example 2: partition rows into k-means slices
#' set.seed(123)  # reproducible k-means slices
#' plot_effect_heatmap(
#'   estimates = list(
#'     fastMetaHD = model_multi$estimate,
#'     REM = model_rem$estimate,
#'     FEM = model_fem$estimate
#'   ),
#'   outcome_names = colnames(Y),
#'   row_km = 4
#' )
#'
#' # Example 3: label the outcomes with the largest absolute effect with connector lines.
#' # useful when there are many outcomes. 
#' plot_effect_heatmap(
#'   estimates = list(
#'     fastMetaHD = model_multi$estimate,
#'     REM = model_rem$estimate,
#'     FEM = model_fem$estimate
#'   ),
#'   outcome_names = colnames(Y),
#'   label_top_n = 15
#' )
#'
#' \dontrun{
#' # Example 4: save to PNG
#' plot_effect_heatmap(
#'   estimates = list(
#'     fastMetaHD = model_multi$estimate,
#'     REM = model_rem$estimate,
#'     FEM = model_fem$estimate
#'   ),
#'   outcome_names = colnames(Y),
#'   file = "effect_heatmap.png",
#'   width = 10,
#'   height = 8,
#'   dpi = 600
#' )
#' }
#'
#' @importFrom grDevices dev.off
#' @importFrom ComplexHeatmap Heatmap draw
#' @importFrom circlize colorRamp2
#' @export
plot_effect_heatmap <- function(estimates,
                                outcome_names = NULL,
                                title = "Pooled Effect Sizes Across Methods",
                                col_low = "blue",
                                col_mid = "white",
                                col_high = "red",
                                cluster_rows = TRUE,
                                cluster_columns = FALSE,
                                row_km = 1,
                                show_row_names = NULL,
                                label_top_n = NULL,
                                row_names_fontsize = 8,
                                column_names_fontsize = 10,
                                legend_title = "Pooled\nEffect Size",
                                file = NULL,
                                width = 8,
                                height = 6,
                                units = c("in", "cm", "mm"),
                                dpi = 300,
                                ...) {
  if (!is.list(estimates) || is.null(names(estimates))) {
    stop("`estimates` must be a named list of numeric vectors.", call. = FALSE)
  }
  if (!all(vapply(estimates, is.numeric, logical(1)))) {
    stop("All elements in `estimates` must be numeric vectors.", call. = FALSE)
  }
  if (length(unique(lengths(estimates))) != 1L) {
    stop("All vectors in `estimates` must have the same length.", call. = FALSE)
  }
  units <- match.arg(units)
  if (!is.null(file)) {
    .open_plot_device(file, width = width, height = height,
                      units = units, dpi = dpi)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  effect_mat <- do.call(cbind, estimates)
  if (!is.null(outcome_names)) {
    if (length(outcome_names) != nrow(effect_mat)) {
      stop(
        "`outcome_names` length (", length(outcome_names), ") ",
        "must match number of outcomes (", nrow(effect_mat), ").",
        call. = FALSE
      )
    }
    rownames(effect_mat) <- outcome_names
  } else {
    rownames(effect_mat) <- paste0("Var", seq_len(nrow(effect_mat)))
  }
  if (!is.numeric(row_km) || length(row_km) != 1L || row_km < 1) {
    stop("`row_km` must be a single integer >= 1.", call. = FALSE)
  }
  row_km <- as.integer(row_km)
  n_out <- nrow(effect_mat)
  mark_idx <- integer(0)
  if (!is.null(label_top_n)) {
    if (!is.numeric(label_top_n) || length(label_top_n) != 1L ||
        label_top_n < 1) {
      stop("`label_top_n` must be a single positive integer.", call. = FALSE)
    }
    score <- apply(abs(effect_mat), 1, max)
    k  <- min(as.integer(label_top_n), n_out)
    mark_idx <- order(score, decreasing = TRUE)[seq_len(k)]
  }
  right_annotation <- NULL
  if (length(mark_idx) > 0L) {
    mark_idx <- sort(unique(mark_idx))
    right_annotation <- ComplexHeatmap::rowAnnotation(
      outcome = ComplexHeatmap::anno_mark(
        at = mark_idx,
        labels = rownames(effect_mat)[mark_idx],
        labels_gp = grid::gpar(fontsize = row_names_fontsize)
      )
    )
  }
  if (is.null(show_row_names)) {
    show_row_names <- (n_out <= 50L) && length(mark_idx) == 0L
    if (!show_row_names && length(mark_idx) == 0L) {
      message(
        "Hiding row labels for ", n_out, " outcomes (> 50). ",
        "Set `show_row_names = TRUE` to force labels."
      )
    }
  }
  col_fun <- circlize::colorRamp2(
    c(min(effect_mat), 0, max(effect_mat)),
    c(col_low, col_mid, col_high)
  )
  ht_args <- list(
    effect_mat,
    name = legend_title,
    col = col_fun,
    column_title = title,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    row_km = row_km,
    show_row_names = show_row_names,
    row_names_gp = grid::gpar(fontsize = row_names_fontsize),
    column_names_gp = grid::gpar(fontsize = column_names_fontsize),
    heatmap_legend_param = list(title = legend_title)
  )
  if (!is.null(right_annotation)) {
    ht_args$right_annotation <- right_annotation
  }
  ht <- do.call(ComplexHeatmap::Heatmap, c(ht_args, list(...)))
  ComplexHeatmap::draw(ht)
  invisible(ht)
}
