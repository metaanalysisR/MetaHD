\name{MetaHD}
\alias{MetaHD}
\title{A Multivariate Meta-Analysis Model for Metabolomics Data
}
\description{
The function \code{MetaHD} performs a multivariate meta-analysis for combining summary estimates obtained from multiple metabolomic studies by using restricted maximum likelihood estimation.
}
\usage{
MetaHD(Y,Slist,Psi = NULL,shrinkCor = TRUE,method = c("reml","fixed"),bscov = c("unstructured","diag"),rigls.maxiter = 5,impute.na = FALSE,impute.var = 10^4)
}
\arguments{
Assuming a meta-analysis is based on N outcomes/metabolites and K studies:
  \item{Y}{treatment effect sizes of the outcomes. This should be in the form of a K \eqn{\times} N matrix.}
  \item{Slist}{K-dimensional list of N \eqn{\times} N matrices representing within-study variances and covariances of the treatment effects.}
  \item{Psi}{N \eqn{\times} N matrix representing between-study variances and covariances of the treatment effects.(optional, if not specified this will be estimated internally by \code{"MetaHD"} using \code{"estimateBSvar"} and \code{"estimateCorMat"} functions in \code{"MetaHD"} package.)}
  \item{shrinkCor}{a logical value indicating whether a shrinkage estimator should be used to estimate between-study correlation matrix. Default is TRUE.}
  \item{method}{estimation method: \code{"fixed"} for fixed-effects models, \code{"reml"} for random-effects models fitted through restricted maximum likelihood}
  \item{bscov}{a character vector defining the structure of the random-effects covariance matrix. Among available covariance structures, the user can select \code{"unstructured"} to obtain between-study covariance matrix with diagonal elements (variances) estimated using restricted maximul likelihood and off-diagonal elements (co-variances) reflecting the correlations estimated via shrinkage and \code{"diag"}(diagonal)}
  \item{rigls.maxiter}{maximum number of iterations of the restricted iterative generalized least square algorithm. Default is set to 5}
  \item{impute.na}{a logical value indicating whether missing values need to be imputed or not. Default is FALSE.}
  \item{impute.var}{multiplier for replacing the missing variances in Slist.(a large value, default is \eqn{{10}^{4}})}
}
\details{
}
\value{
The \code{MetaHD} function returns a list of objects:
  \item{estimate }{a N-dimensional vector of the combined estimates}
  \item{std.err }{a N-dimensional vector of the associated standard errors}
  \item{I2.stat}{\eqn{I^2} statistic}
}
\references{
Sera F, Gasparrini A. (2019). An extended mixed-effects framework for meta-analysis.\emph{Statistics in Medicine}. 2019;38(29):5429-5444.

Sch\eqn{\"a}fer, J., and K. Strimmer. (2005). A shrinkage approach to large-scale covariance estimation and implications for functional genomics. \emph{Statistical Applications in Genetics and Molecular Biology} 2005; 4:32.

}
\author{
Jayamini C. Liyanage and Alysha M. De Livera
}
