# LOAD LIBRARIES
library(corpcor) #cor.shrink()

estimateCorMat <- function(Y,shrinkCor=TRUE){
  N <- ncol(Y) # NO.OF METABOLITES/OUTCOMES
  K <- nrow(Y) # NO.OF STUDIES
  # IF THE NO.OF STUDIES ARE LESS THAN OR EQUAL TO 2, SET BETWEEN-STUDY CORRELATIONS TO BE ZERO
  if (K <= 2){
    cormat <- 0
  }else {
    if (N > K){
      if(shrinkCor){
        # ESTIMATE CORRELATIONS VIA SHRINKAGE
        cormat <- cor.shrink(Y,verbose = FALSE)[1:N,1:N]
      }else{
        cormat <- cor(Y)
      }
    }else{
      cormat <- cor(Y)
    }
  }
  return(cormat)
}
