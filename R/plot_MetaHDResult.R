#' Plot `MetaHD` results
#'
#' Visualizes a \code{MetaHDResult} object using a Venn diagram (Tav, 2025)
#' for a compact view of overlapping outcomes across methods, an UpSet plot 
#' (Conway et al., 2017) to show overlap of significant outcomes across methods, 
#' or ROC curves (Robin et al., 2011) to assess performance when ground truth is 
#' available. Optionally saves the plot to a file in any of several formats.
#'
#' @param x An object of class \code{"MetaHDResult"}.
#' @param type Character string specifying the type of plot:
#'   \code{"upset"} (default) for an UpSet plot (Conway et al., 2017);
#'   \code{"venn"} for a Venn diagram (Tav, 2025); or \code{"ROC"} for ROC
#'   curves (Robin et al., 2011). See Details for more information.
#' @param roc.colors Optional vector of colors for ROC curves. Must match the
#'   number of methods. If \code{NULL} (default), colors are generated
#'   automatically via \code{grDevices::hcl.colors()}. Used only when
#'   \code{type = "ROC"}.
#' @param venn.colors Optional vector of colors for the Venn diagram sets.
#'   Must match the number of methods. If \code{NULL} (default), colors are
#'   generated automatically via \code{grDevices::hcl.colors()}. Used only
#'   when \code{type = "venn"}.
#' @param queries For \code{type = "upset"}: a list of pre-built UpSetR query
#'   lists for highlighting, or \code{NULL}. For details and examples
#'   on creating UpSetR queries see R package \pkg{UpSetR} (Conway et al.,
#'   2017).
#' @param highlight For \code{type = "upset"}: a convenience argument for
#'   highlighting intersections. A list of character vectors; each vector
#'   names the methods (or \code{"Truth"}) whose intersection should be
#'   highlighted. Ignored if \code{queries} is supplied.
#' @param highlight.colors Optional colors for \code{highlight}. If
#'   \code{NULL} (default), distinct colors are generated automatically
#'   via \code{grDevices::hcl.colors()}.
#' @param show.truth Logical. If \code{TRUE} (default) and \code{x$truth} is
#'   available, includes \code{"Truth"} as an additional set in the UpSet
#'   plot. Not applicable for \code{type = "ROC"} or \code{type = "venn"}.
#' @param file Optional file path to save the plot. If \code{NULL} (default),
#'   the plot is drawn to the active graphics device. See the
#'   \strong{Saving the plot} section below.
#' @param width,height Plot dimensions. \emph{Used only when \code{file} is
#'   set.} Default 8 x 6.
#' @param units Units for \code{width} and \code{height}: \code{"in"}
#'   (default), \code{"cm"}, or \code{"mm"}. \emph{Used only when \code{file}
#'   is set.}
#' @param dpi Resolution for raster formats (PNG, JPEG, TIFF, BMP). Ignored
#'   for vector formats (PDF, SVG). Default 300. \emph{Used only when
#'   \code{file} is set.}
#' @param ... Further arguments passed to \code{gVenn::plotVenn()} for 
#'   \code{type = "venn"}, \code{UpSetR::upset()} for \code{type = "upset"}, 
#'   or \code{graphics::plot()} for \code{type = "ROC"}.
#'
#' @details
#' \strong{Venn diagram} (Tav, 2025): Displays proportional overlap of significant 
#' or selected outcomes across methods, implemented via the \pkg{gVenn} package.
#' Recommended for up to three methods; for more methods use \code{type = "upset"}. 
#' Requires R >= 4.5 and \pkg{gVenn} to be installed:
#' \preformatted{
#'   BiocManager::install("gVenn")
#' }
#' Available for both input modes supported by \code{\link{MetaHDResult}}.
#'
#' \strong{UpSet plot} (Conway et al., 2017): Visualizes the intersection of
#' significant or selected outcomes across meta-analysis methods, implemented
#' via the \pkg{UpSetR} package. Available for both input modes supported
#' by \code{\link{MetaHDResult}}.
#'
#' \strong{ROC curves} (Robin et al., 2011): Requires a ground truth vector
#' (provided when creating the object) and displays ROC curves along with AUC
#' values for each method, implemented via the \pkg{pROC} package. AUC values
#' are computed as the area under each ROC curve and displayed in the plot
#' legend. \strong{ROC curves require numeric p-values} --- they are not
#' available when \code{top_n} was used to create the \code{MetaHDResult}
#' object.
#'
#' @section Saving the plot:
#' To save the plot directly, pass a file path to \code{file}. The graphics
#' device is auto-selected from the file extension, the plot is written to
#' that file, and the plot is \emph{not} also drawn to the screen. Supported
#' extensions:
#'
#' \itemize{
#'   \item Vector formats: \code{.pdf}, \code{.svg}
#'   \item Raster formats: \code{.png}, \code{.jpeg}/\code{.jpg},
#'         \code{.tiff}/\code{.tif}, \code{.bmp}
#' }
#'
#' Use \code{width}, \code{height}, \code{units}, and \code{dpi} to control
#' the saved output. \code{dpi} is ignored for vector formats. The path may
#' be absolute or relative to \code{getwd()}; missing parent directories are
#' created automatically.
#'
#' @return Invisibly returns: a named list of character vectors containing the
#'   significant or top-ranked outcome names per method for \code{type = "venn"};
#'   and a named numeric vector of AUCs for \code{type = "ROC"}.
#'
#' @seealso
#' \code{\link{MetaHDResult}} for creating objects to be plotted.
#' \code{\link{plot_effect_heatmap}} for comparing pooled effect size
#' estimates across methods.
#' \code{\link{plot_correlation_heatmap}} for visualising correlations
#' among outcomes.
#'
#' @examples
#' set.seed(123)
#' N <- 100
#' truth <- rbinom(N, 1, 0.2)   # 20% of features are true signals
#' res <- MetaHDResult(
#'   sets = list(
#'     method_A = runif(N)^ifelse(truth, 5, 1),
#'     method_B = runif(N)^ifelse(truth, 3, 1),
#'     method_C = runif(N)^ifelse(truth, 2, 1)
#'   ),
#'   truth = truth
#' )
#'
#' # Example 1: simple UpSet plot 
#' plot(res, type = "upset")
#'
#' # Example 2: highlights via the `highlight` argument
#' plot(res, type = "upset",
#'      highlight = list(c("Truth", "method_A", "method_B", "method_C")),
#'      highlight.colors = "darkgreen")
#'
#' # Example 3: passing in-built queries
#' plot(res, type = "upset",
#'      queries = list(
#'        list(
#'          query = UpSetR::intersects,
#'          params = list(c("method_A", "method_B", "method_C")),
#'          color = "dodgerblue3",
#'          active = TRUE,
#'          query.name = "Identified by all methods"
#'        ),
#'        list(
#'          query = function(row) {
#'            row["method_A"] == 1 && sum(row) < length(row)
#'          },
#'          color = "orange",
#'          active = TRUE,
#'          query.name = "Others identified by method_A"
#'        )
#'      ),
#'      show.truth = FALSE)
#'
#' # Example 4: ROC curves with AUC values 
#' # Requires numeric p-values in sets
#' aucs <- plot(res, type = "ROC")
#' print(aucs)
#'
#' # Example 5: Venn diagram (Tav, 2025)
#' # Requires R >= 4.5 and gVenn installed
#' if (requireNamespace("gVenn", quietly = TRUE)) {
#'   plot(res, type = "venn")
#' }
#'
#' \dontrun{
#' # Example 6: save plots to file
#' plot(res, type = "venn",  file = "venn.pdf",  width = 8,  height = 7)
#' plot(res, type = "upset", file = "upset.pdf", width = 10, height = 6)
#' plot(res, type = "ROC",   file = "roc.png",   width = 8,  height = 6, dpi = 600)
#' }
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
#' @method plot MetaHDResult
#' @export
plot.MetaHDResult <- function(x,
                              type = c("upset", "venn", "ROC"),
                              roc.colors = NULL,
                              venn.colors = NULL,
                              queries = NULL,
                              highlight = NULL,
                              highlight.colors = NULL,
                              show.truth = TRUE,
                              file = NULL,
                              width = 8,
                              height = 6,
                              units = c("in", "cm", "mm"),
                              dpi = 300,
                              ...) {
  type <- match.arg(type)
  units <- match.arg(units)
  if (!is.null(file)) {
    .open_plot_device(file, width = width, height = height,
                      units = units, dpi = dpi)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  if (type == "upset") {
    plot_df <- x$sig_df
    if (isTRUE(show.truth) && !is.null(x$truth)) {
      plot_df <- cbind(Truth = x$truth, plot_df)
    }
    set_names <- colnames(plot_df)
    if (is.null(queries) && !is.null(highlight)) {
      queries <- .build_intersect_queries(
        highlight = highlight,
        available_sets = set_names,
        highlight_colors = highlight.colors
      )
    }
    print(withCallingHandlers(
      UpSetR::upset(
        plot_df,
        sets = set_names,
        keep.order = TRUE,
        mainbar.y.label = "Shared Significant Outcomes",
        sets.x.label = "Significant Outcomes\nPer Method",
        order.by = "freq",
        queries = queries,
        query.legend = if (!is.null(queries)) "top" else "none",
        ...
      ),
      warning = function(w) {
        msg <- conditionMessage(w)
        if (grepl("deprecat|aes_string|linewidth|empty aesthetic|Ignoring unknown",
                  msg, ignore.case = TRUE)) {
          invokeRestart("muffleWarning")
        }
      }
    ))
    return(invisible(NULL))
  }
  if (type == "ROC") {
    if (is.null(x$truth)) {
      stop("True signals must be provided to compute ROC curves.",
           call. = FALSE)
    }
    if (!is.null(x$top_n)) {
      stop(
        "ROC curves require numeric p-values. ",
        "They are not available when `top_n` was used to create the ",
        "MetaHDResult object.",
        call. = FALSE
      )
    }
    pvals_list <- lapply(x$sets, function(p) pmax(p, 1e-10))
    logp_list <- lapply(pvals_list, function(p) -log10(p))
    roc_list <- lapply(logp_list, function(lp) {
      pROC::roc(x$truth, lp, quiet = TRUE)
    })
    if (is.null(roc.colors)) {
      cols <- grDevices::hcl.colors(length(roc_list), "Dark 3")
    } else if (length(roc.colors) != length(roc_list)) {
      stop("Length of `roc.colors` must match number of methods.",
           call. = FALSE)
    } else {
      cols <- roc.colors
    }
    plot(0, 0, type = "n",
         xlim = c(0, 1), ylim = c(0, 1),
         xlab = "False Positive Rate",
         ylab = "True Positive Rate",
         main = "ROC Curves for Different Meta-Analysis Methods",
         cex.lab = 0.8, cex.axis = 0.8, cex.main = 0.9,
         ...)
    graphics::abline(0, 1, lty = 2, col = "gray")
    for (i in seq_along(roc_list)) {
      graphics::lines(1 - roc_list[[i]]$specificities,
                      roc_list[[i]]$sensitivities,
                      col = cols[i], lwd = 2)
    }
    auc_vals <- vapply(roc_list, function(r) as.numeric(pROC::auc(r)),
                              numeric(1))
    names(auc_vals) <- names(x$sets)
    graphics::legend("bottomright",
                     legend = paste0(names(auc_vals),
                                     " (AUC = ", round(auc_vals, 3), ")"),
                     col = cols, lwd = 2, cex = 0.8, bty = "n")

    return(invisible(auc_vals))
  }
  if (type == "venn") {
    if (!requireNamespace("gVenn", quietly = TRUE)) {
      stop(
        "Package 'gVenn' is required for Venn diagrams.\n",
        "It requires R >= 4.5. Install it with:\n",
        "  BiocManager::install('gVenn')",
        call. = FALSE
      )
    }
    if (length(x$sets) > 3) {
      warning(
        "Venn diagrams with more than 3 sets can be hard to read. ",
        "Consider using type = 'upset' instead.",
        call. = FALSE
      )
    }
    sig_sets <- lapply(names(x$sets), function(nm) {
      rownames(x$sig_df)[x$sig_df[[nm]] == 1L]
    })
    names(sig_sets) <- names(x$sets)
    if (is.null(venn.colors)) {
      venn.colors <- grDevices::hcl.colors(length(x$sets), "Dark 3")
    } else if (length(venn.colors) != length(x$sets)) {
      stop(
        "`venn.colors` length (", length(venn.colors), ") ",
        "must match number of methods (", length(x$sets), ").",
        call. = FALSE
      )
    }
    overlaps <- gVenn::computeOverlaps(sig_sets)
    print(gVenn::plotVenn(overlaps,
                          fills = venn.colors,
                          ...))
    return(invisible(sig_sets))
  }
}
