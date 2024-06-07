\name{realdata}
\alias{realdata}
\docType{data}
\title{
An Individual-Level Metabolomics Dataset
}
\description{
This dataset consists of a list of three data frames containing individual-level data, treatment effect estimates, and their associated variances.
}
\usage{realdata}
\format{
  A list of data frames as follows:
  \describe{
    \item{\code{all}}{A dataframe with 12 rows and 14 columns containing individual-level data of 14 metabolites in columns and a total of 12 samples from two groups that have been run separately in two different instruments, leading to two separate studies in rows.}
    \item{\code{effects}}{A dataframe with 2 rows and 14 columns, containing treatment effect estimates of the 14 metabolites for the 2 studies.}
    \item{\code{var}}{A dataframe with 2 rows and 14 columns, containing associated variances of the treatment effects.}
  }
}

\examples{
head(realdata$all)
head(realdata$effects)
head(realdata$var)
}
\keyword{datasets}
