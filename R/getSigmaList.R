getSigmaList <- function (Zlist, nalist, Psi, Slist){
  if (is.null(Psi))
    return(Slist)
  Psi <- getList(Psi)
  if (is.null(Zlist))
    return(mapply(function(S, na) S + Psi[[1L]][!na, !na, drop = FALSE], Slist, nalist, SIMPLIFY = FALSE))
  Psi <- getList(Psi)
  lapply(seq_along(Zlist), function(i) sumList(lapply(seq_along(Psi),
                                                      function(j) blockDiagMat(lapply(Zlist[[i]][[j]], function(x) x %*%
                                                                                    Psi[[j]] %*% t(x))))) + Slist[[i]])
}
