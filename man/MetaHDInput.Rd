\name{MetaHDInput}
\alias{MetaHDInput}
\title{
Creating Input Data for MetaHD When Individual-Level Data are Available
}
\description{
The MetaHDInput function creates input data Y (treatment effects) and Slist (within-study covariance matrices) for MetaHD when individual-level data are available. Assuming that the individual-level data are in the following format, with 'study' in column 1, 'group' in column 2 and outcomes in rest of the columns, with samples in rows.
}
\usage{
MetaHDInput(data,isMissing=FALSE,isLgTransformed=FALSE)
}
\arguments{
  \item{data}{ a dataframe consisting of individual-level data in the format, where 'study' in column 1, 'group' in column 2 and outcomes in rest of the columns and samples in rows.}
  \item{isMissing}{ a logical value indicating whether data are missing or not. Default is FALSE.}
  \item{isLgTransformed}{ a logical value indicating whether data are log transformed or not. Default is FALSE.}
}
\value{
A list of objects containing :
 \item{Y}{treatment effect sizes of the outcomes in the form of a K x N matrix, where K is the number of studies and N is the number of outcomes.}
 \item{Slist}{K-dimensional list of N x N matrices representing within-study variances and covariances of the treatment effects}
}
\examples{
input_data <- MetaHDInput(realdata$all)

Y <- input_data$Y
Slist <- input_data$Slist

## MULTIVARIATE RANDOM-EFFECTS META-ANALYSIS, ESTIMATED WITH REML
model <- MetaHD(Y, Slist, method = "reml", bscov = "unstructured")
model$estimate
model$pVal
}
