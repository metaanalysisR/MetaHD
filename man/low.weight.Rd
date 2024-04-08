\name{low.weight}
\alias{low.weight}
\title{Replace Missing Values
}
\description{
This function augment data by replacing missing values in both outcomes and their associated covariances. It replaces missing outcomes and covariances with 0, while missing variances are replaced with a value obtained by multiplying the largest observed variance by impute.var. Typically, impute.var is set to a significantly high value, by default \eqn{{10}^{4}}, so the missing outcome is allocated a lower weight in the meta analysis. This functionality is employed internally within the \code{"MetaHD"} function to impute missing outcomes, variances and covariances for subsequent use in model estimation.
}
\usage{
low.weight(Y, Slist, impute.var=10^4)
}
\arguments{
  \item{Y}{treatment effect sizes of the outcomes. This should be in the form of a K \eqn{\times} N matrix.}
  \item{Slist}{K-dimensional list of N \eqn{\times} N matrices representing within-study variances and covariances of the treatment effects.}
  \item{impute.var}{multiplier for replacing the missing variances in Slist.(a large value, default is \eqn{{10}^{4}})}
}
\details{
}
\value{
The \code{low.weight} function returns a list of objects:
  \item{effects}{a K \eqn{\times} N matrix corresponding to the augmented outcomes}
  \item{wscovar}{a K-dimensional list of N \eqn{\times} N matrices corresponding to the augmented within-study variances and covariances of the treatment effects}
}
\references{
}
\author{
Jayamini C. Liyanage and Alysha M. De Livera
}
\note{
}

\seealso{
See \code{\link{inputna}} in \code{\link{mixmeta}} package
}
\examples{

}
