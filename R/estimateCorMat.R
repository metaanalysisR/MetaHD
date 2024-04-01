estimateCorMat <- function(Y,shrinkCor=TRUE){
  N <- ncol(Y)
  K <- nrow(Y)
  if (K <= 2){
    cormat <- 0
  }else {
    if (N > K){
      if(shrinkCor){
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
