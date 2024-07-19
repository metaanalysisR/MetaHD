\name{simdata.2}
\alias{simdata.2}
\docType{data}
\title{
Simulated Dataset 2 : With Data Missing-At-Random
}
\description{
This dataset consists of a list of two data frames containing treatment effect-sizes and within-study covariance matrices with missing values
}
\usage{simdata.2}
\format{
   A list of data frames as follows:
  \describe{
    \item{\code{Y}}{treatment effect sizes of the metabolites in the form of a 12 x 30 matrix, where 12 is the number of studies and 30 is the number of metabolites.}
    \item{\code{Slist}}{12-dimensional list of 30 x 30 matrices representing within-study variances and covariances of the treatment effects}
  }
}

\examples{
Y <- simdata.2$Y
Slist <- simdata.2$Slist

head(Y)
head(Slist[[1]])
head(Slist[[12]])
}
\keyword{datasets}
