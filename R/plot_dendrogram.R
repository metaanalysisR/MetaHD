#' Plot dendrogram 
#'
#' Plots the dendrogram of outcomes to help users choose an appropriate
#' cutting height when using \code{dendro.method = "fixedHeight"} in
#' \code{MetaHD}.
#'
#' @param Y A K x N matrix of treatment effect sizes of the outcomes (K studies, N outcomes).
#' @param shrinkCor Logical. Whether to use shrinkage estimation for the
#'   correlation matrix. Default is \code{TRUE}.
#' @param impute.na Logical. Whether to impute missing values. Default is
#'   \code{FALSE}.
#' @param h Optional numeric value. If provided, draws a horizontal red
#'   dashed line at this height to preview the cut.
#' @param ... Further arguments passed to \code{plot()}.
#'
#' @return Invisibly returns the \code{hclust} object.
#' @export
#' 
#' @importFrom graphics abline text

plot_dendrogram <- function(Y, shrinkCor = TRUE, impute.na = FALSE, 
                            h = NULL, ...) {
  cormat <- estimateCorMat(Y, shrinkCor, impute.na)
  distMat <- as.dist(1 - abs(cormat))
  hc <- hclust(distMat, method = "complete")
  plot(hc,
       main = "Dendrogram",
       xlab = "Outcomes",
       ylab = "Height",
       sub = "Use this plot to choose a cutting height for dendro.method = 'fixedHeight'",
       cex = 0.7,
       ...)
  if (!is.null(h)) {
    graphics::abline(h = h,
                     col = "red",
                     lty = 2,
                     lwd = 1.5)
    graphics::text(x = 1,
                   y = h + 0.01,
                   labels = paste0("h = ", h),
                   col = "red",
                   cex = 0.8,
                   adj = 0)
  }
  message("Dendrogram height range: [",
          round(min(hc$height), 3), ", ",
          round(max(hc$height), 3), "]",
          "\nUse the 'h' argument to preview a cut at a specific height.")
  invisible(hc)
}