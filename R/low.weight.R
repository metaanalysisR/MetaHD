low.weight <- function(Y,Slist,impute.var=10^4){
  N <- ncol(Y)
  K <- nrow(Y)
  Y[is.na(Y)] <- 0
  S <- t(sapply(Slist,vechMatrix))
  n <- 1
  for (j in 1:ncol(S)) {
    if (j == (((2*N + 3)*n - (n^2) - (2*N))/2)){
      S[,j][is.na(S[,j])] <- max(S[,j],na.rm = TRUE)*impute.var
      n <- n + 1
    }else{
      S[,j][is.na(S[,j])] <- 0
    }
  }
  Slist <- lapply(seq(nrow(S)), function(i) xpndMatrix(S[i,]))
  return(list(effects = Y,
              wscovar = Slist))
}
