\name{simdata.1}
\alias{simdata.1}
\docType{data}
\title{
Simulated Dataset 1 : With Complete Data
}
\description{
This dataset consists of a list of two data frames containing treatment effect-sizes and within-study covariance matrices
}
\usage{simdata.1}
\format{
   A list of data frames as follows:
  \describe{
    \item{\code{Y}}{treatment effect sizes of the metabolites in the form of a 12 x 30 matrix, where 12 is the number of studies and 30 is the number of metabolites.}
    \item{\code{Slist}}{12-dimensional list of 30 x 30 matrices representing within-study variances and covariances of the treatment effects}
  }
}

\examples{
Y <- simdata.1$Y
Slist <- simdata.1$Slist

head(Y)
head(Slist[[1]])
head(Slist[[12]])
}
\keyword{datasets}
