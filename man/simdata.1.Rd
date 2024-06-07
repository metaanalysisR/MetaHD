\name{simdata.1}
\alias{simdata.1}
\docType{data}
\title{
Simulated Dataset 1 : With Complete Data
}
\description{
This dataset consists of a list of three data frames containing treatment effect estimates, within-study variances and upper trainangular elements of the within-study correlation matrices.
}
\usage{simdata.1}
\format{
   A list of data frames as follows:
  \describe{
    \item{\code{effects}}{A dataframe with 12 rows and 30 columns, containing simulated treatment effect estimates of 30 metabolites for 12 studies.}
    \item{\code{wsvar}}{A dataframe with 12 rows and 30 columns, containing simulated within-study variances of 30 metabolites for 12 studies.}
    \item{\code{wscor}}{A dataframe with 12 rows and 435 columns, containing the upper traingular elements of simulated within-study correlation matrices for 12 studies.}
  }
}

\examples{
head(simdata.1$effects)
head(simdata.1$wsvar)
head(simdata.1$wscor)
}
\keyword{datasets}
