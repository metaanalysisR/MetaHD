#' HANDLE MISSING VALUES

low.weight <- function(Y,Slist,impute.var=10^4){
  N <- ncol(Y)
  K <- nrow(Y)
  # REPLACES MISSING VALUES IN THE MATRIX Y (OUTCOMES) WITH ZEROS.
  Y[is.na(Y)] <- 0
  # APPLY THE FUNCTION VechMatrix TO EACH MATRIX IN THE Slist AND TRANSPOSES THE RESULT
  S <- t(sapply(Slist,vechMatrix))
  n <- 1
  for (j in 1:ncol(S)) {
    # CHECK FOR THE COLUMNS THAT CORRESPOND TO THE ASSOCIATED VARIANCES OF THE OUTCOMES
    if (j == (((2*N + 3)*n - (n^2) - (2*N))/2)){
      # REPLACES THE MISSING VARIANCES WITH THE MAXIMUM OBSERVED VALUE IN THAT COLUMN, MULTIPLIED BY THE VALUE OF impute.var
      S[,j][is.na(S[,j])] <- max(S[,j],na.rm = TRUE)*impute.var
      n <- n + 1
    }else{
      # REPLACES THE MISSING CO-VARIANCES WITH ZERO
      S[,j][is.na(S[,j])] <- 0
    }
  }
  # APPLY THE XPNDMATRIX FUNCTION TO EACH ROW OF THE MATRIX S AND EXPAND TO A LIST
  Slist <- lapply(seq(nrow(S)), function(i) xpndMatrix(S[i,]))
  return(list(effects = Y,
              wscovar = Slist))
}
