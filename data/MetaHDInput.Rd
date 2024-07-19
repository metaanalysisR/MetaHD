\name{MetaHDInput}
\alias{MetaHDInput}
\title{
Creating Input Data for MetaHD When Individual-Level Data are Available
}
\description{
The MetaHDInput function creates input data Y (treatment effects) and Slist (within-study covariance matrices) for MetaHD when individual-level data are available. Assuming that the individual-level data are in the following format, with 'study' in column 1, 'group' in column 2 and outcomes in rest of the columns, with samples in rows.
}
\usage{
MetaHDInput(data)
}
\arguments{
  \item{data}{ a dataframe consisting of individual-level data in the format, where 'study' in column 1, 'group' in column 2 and outcomes in rest of the columns and samples in rows.}
}
\value{
A list of objects containing :
 \item{Y}{treatment effect sizes of the outcomes in the form of a K x N matrix, where K is the number of studies and N is the number of outcomes.}
 \item{Slist}{K-dimensional list of N x N matrices representing within-study variances and covariances of the treatment effects}
}
\examples{
input_data <- MetaHDInput(realdata)

Y <- input_data$Y
Slist <- input_data$Slist

## MULTIVARIATE RANDOM-EFFECTS META-ANALYSIS, ESTIMATED WITH REML
model <- MetaHD(Y, Slist, method = "reml", bscov = "unstructured")
model$estimate
model$pVal
}
