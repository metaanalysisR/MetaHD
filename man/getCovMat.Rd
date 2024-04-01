\name{getCovMat}
\alias{getCovMat}
\title{
Create (Co)variance Matrices
}
\description{
This function create (co)variance matrices of a set of outcomes when the corresponding standard deviation and correlation values are given.
}
\usage{
getCovMat(sd, cor = NULL)
}

\arguments{
  \item{sd}{
  a k \eqn{\times} n matrix of standard deviations for n outcomes in k studies, or a vector for n outcomes
}
  \item{cor}{
  either a vector of length 1, k or n(n-1)/2, or alternatively a n \eqn{\times} n or k \eqn{\times} n(n-1)/2 matrix. See Details.
}
}
\details{
Depending on the number of outcomes n and studies k, the argument cor is interpreted as:

if a vector of length 1 (a scalar), the same correlation for all the n outcomes for all the k studies;

if a vector of length k, the same correlation for all the n outcomes for each of the k studies;

if a vector of length n(n-1)/2, the lower triangular elements (without diagonal, taken by column) of the correlation matrix of the n outcomes, the same for all the k studies;

if a n \eqn{\times} n matrix, the correlation matrix of the n outcomes, the same for all the k studies;

if a k \eqn{\times} n(n-1)/2 matrix, each row represents the lower triangular elements (without diagonal, taken by column) of the correlation matrix of the n outcomes for each of the k studies.
}
\value{
(co)variance marix
}

\author{
Jayamini C. Liyanage
}
\seealso{
See \code{\link{inputcov}} in \code{\link{mixmeta}} package
}
