\name{xpndMatrix}
\alias{xpndMatrix}
\title{Vectorization and Expansion of Symmetric Matrices
}
\description{
The vechMatrix function converts a symmetric matrix into a vector by extracting its lower triangular elements column-wise. Conversely, the xpndMatrix function performs the inverse operation, restoring the original matrix from the vector.
}
\usage{
vechMatrix(mat, diag = TRUE)
xpndMatrix(vech)
}

\arguments{
  \item{mat}{a square matrix.}
  \item{vech}{a vector.}
  \item{diag}{a logical switch indicating if the diagonal entries must be included.}
}
\value{
A vector containing lower triangular elements of the square matrix column-wise for vechMatrix, a symmetric matrix for xnpdMatrix.
}
\author{
Jayamini C. Liyanage
}
\seealso{
See functions \code{\link{vechMat}} and \code{\link{xpndMat}} in \code{\link{mixmeta}} package
}
\examples{
# GENERATE A POSITIVE-DEFINITE MATRIX, VECTORIZE IT AND THEN RE-EXPAND
(M <- crossprod(matrix(rnorm(9),3)))
(v <- vechMatrix(M))
xpndMatrix(v)
}
