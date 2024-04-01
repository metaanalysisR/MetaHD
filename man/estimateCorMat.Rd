\name{estimateCorMat}
\alias{estimateCorMat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimate Between-Study Correlation Matrix
}
\description{
This function estimate between-study correlation matrix using the available treatment effect estimates. This is used internally in \code{"MetaHD"}.
}
\usage{
estimateCorMat(Y,shrinkCor=TRUE)
}
\arguments{
Assuming a meta-analysis is based on N outcomes/metabolites and K studies:
  \item{Y}{treatment effect sizes of the outcomes. This should be in the form of a K \eqn{\times} N matrix.}
  \item{shrinkCor}{a logical value indicating whether a shrinkage estimator should be used to estimate between-study correlation matrix. Default is TRUE.}
}
\details{
This function estimate between-study correlations using the available treatment effect estimates. Specially, if N > K and if a shrinkage estimator is needs to be used, correlations are estimated using \code{"cor.shrink"} function in \code{"corpcor"} package. If not, standard \code{"cov"} function in R is used.
}
\value{
A N \eqn{\times} N matrix representing between-study correlation matrix.
}
\references{
Sch\eqn{\"a}fer, J., and K. Strimmer. (2005). A shrinkage approach to large-scale covariance estimation and implications for functional genomics. \emph{Statistical Applications in Genetics and Molecular Biology} 2005; 4:32.
}
\author{
Jayamini C. Liyanage and Alysha M. De Livera
}
