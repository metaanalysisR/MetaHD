\name{simdata.2}
\alias{simdata.2}
\docType{data}
\title{
Simulated Dataset 2 : With Data Missing-At-Random
}
\description{
This dataset consists of a list of three data frames containing treatment effect estimates and within-study variances with missing values and upper trainangular elements of the within-study correlation matrices.
}
\usage{simdata.2}
\format{
  A list of data frames as follows:
  \describe{
    \item{\code{effects}}{A dataframe with 12 rows and 30 columns, containing simulated treatment effect estimates of 30 metabolites for 12 studies with missing values.}
    \item{\code{wsvar}}{A dataframe with 12 rows and 30 columns, containing simulated within-study variances of 30 metabolites for 12 studies with missing values.}
    \item{\code{wscor}}{A dataframe with 12 rows and 435 columns, containing the upper traingular elements of simulated within-study correlation matrices for 12 studies.}
  }
}

\examples{
head(simdata.2$effects)
head(simdata.2$wsvar)
head(simdata.2$wscor)
}
\keyword{datasets}
