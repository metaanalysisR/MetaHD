#' Simulated Study-Level Complete Data
#'
#' A simulated dataset for demonstrating meta-analysis workflows in
#' \pkg{MetaHD}. It contains the observed treatment effects, within-study
#' covariance matrices, the associated variances of the treatment effects,
#' individual two-sided \emph{p}-values, true effect sizes, and sample
#' sizes for each study. 
#'
#' @format A list with six elements:
#' \describe{
#'   \item{Y}{A numeric matrix of dimension 10 x 200 containing the
#'     observed effect sizes. Each row corresponds to a study and each
#'     column corresponds to one of 200 outcomes.}
#'   \item{Slist}{A list of length 10 containing the within-study
#'     covariance matrices, one per study.}
#'   \item{wsvar}{A numeric matrix of dimension 10 x 200 containing the
#'     within-study variances for each study and outcome.}
#'   \item{pvals}{A numeric matrix of dimension 10 x 200 containing the
#'     two-sided p-values for each study and outcome.}
#'   \item{true.theta}{A numeric vector of length 200 containing the
#'     true effect sizes.}
#'   \item{sample.size}{An integer vector of length 10 giving the sample
#'     size of each study.}
#' }
#' 
#' @seealso \code{\link{simData.missing}} for the version with missing values.
"simData.complete"


#' Simulated Study-Level Data with Missing Values
#'
#' A simulated dataset for demonstrating meta-analysis workflows in
#' \pkg{MetaHD} when some outcomes are unobserved in some studies. The
#' structure matches \code{\link{simData.complete}}, with missing values
#' generated under a missing-at-random (MAR) mechanism using a logistic
#' regression model that determines the probability of missingness based
#' on the values of other outcomes. Missing values were introduced on ten
#' randomly selected outcomes; for each selected outcome, approximately
#' 50\% of the values were set to missing. 
#'
#' @format A list with six elements:
#' \describe{
#'   \item{Y}{A numeric matrix of dimension 10 x 200 containing the
#'     observed effect sizes, with \code{NA} values where outcomes are
#'     missing. Each row corresponds to a study and each column to one
#'     of 200 outcomes.}
#'   \item{Slist}{A list of length 10 containing the within-study
#'     covariance matrices, one per study, with \code{NA} values where
#'     outcomes are missing.}
#'   \item{wsvar}{A numeric matrix of dimension 10 x 200 containing the
#'     within-study variances for each study and outcome, with \code{NA}
#'     values where outcomes are missing.}
#'   \item{pvals}{A numeric matrix of dimension 10 x 200 containing the
#'     two-sided p-values for each study and outcome, with \code{NA}
#'     values where outcomes are missing.}
#'   \item{true.theta}{A numeric vector of length 200 containing the
#'     true effect sizes.}
#'   \item{sample.size}{An integer vector of length 10 giving the sample
#'     size of each study.}
#' }
#' 
#' @seealso \code{\link{simData.complete}} for the version without missing values.
"simData.missing"