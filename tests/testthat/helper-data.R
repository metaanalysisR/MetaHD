# Small test dataset
make_test_data <- function(K = 10, N = 100, seed = 123) {
  set.seed(seed)
  
  # Y: K x N matrix of effect sizes
  Y <- matrix(rnorm(K * N), nrow = K, ncol = N)
  rownames(Y) <- paste0("Study", seq_len(K))
  colnames(Y) <- paste0("Var",   seq_len(N))
  
  # Slist: K-dimensional list of N x N matrices
  Slist <- lapply(seq_len(K), function(i) {
    S <- diag(runif(N, 0.1, 1))
    rownames(S) <- colnames(S) <- colnames(Y)
    S
  })
  
  list(Y = Y, Slist = Slist, K = K, N = N)
}