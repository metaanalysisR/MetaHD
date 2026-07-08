# ---- MetaHDResult construction ----

test_that("MetaHDResult builds binary indicators in p-value mode", {
  set.seed(1)
  N <- 30
  res <- MetaHDResult(sets = list(A = runif(N), B = runif(N)), alpha = 0.05)

  expect_s3_class(res, "MetaHDResult")
  expect_equal(nrow(res$sig_df), N)
  expect_equal(ncol(res$sig_df), 2L)
  expect_true(all(unlist(res$sig_df) %in% c(0L, 1L)))
})

test_that("MetaHDResult top_n selects exactly top_n outcomes per method", {
  set.seed(1)
  N <- 30
  res <- MetaHDResult(sets = list(A = rnorm(N), B = rnorm(N)), top_n = 5)

  expect_equal(unname(colSums(res$sig_df)), c(5, 5))
  expect_equal(res$top_n, 5L)
})

test_that("MetaHDResult uses default positional outcome names", {
  set.seed(1)
  res <- MetaHDResult(sets = list(A = runif(10)))
  expect_equal(rownames(res$sig_df), as.character(seq_len(10)))
})

test_that("MetaHDResult applies provided outcome names", {
  set.seed(1)
  nm <- paste0("Var", seq_len(10))
  res <- MetaHDResult(sets = list(A = runif(10)), outcome_names = nm)
  expect_equal(rownames(res$sig_df), nm)
})

# ---- Input validation ----

test_that("MetaHDResult rejects character (pre-defined) sets", {
  expect_error(
    MetaHDResult(sets = list(A = c("v1", "v2"), B = c("v2", "v3"))),
    "must contain numeric"
  )
})

test_that("MetaHDResult errors on unequal vector lengths", {
  expect_error(
    MetaHDResult(sets = list(A = runif(10), B = runif(8))),
    "same length"
  )
})

test_that("MetaHDResult errors on outcome_names length mismatch", {
  expect_error(
    MetaHDResult(sets = list(A = runif(10)), outcome_names = letters[1:3]),
    "must match length"
  )
})

test_that("MetaHDResult validates alpha and top_n", {
  expect_error(MetaHDResult(sets = list(A = runif(10)), alpha = 1.5), "alpha")
  expect_error(MetaHDResult(sets = list(A = rnorm(10)), top_n = 0), "top_n")
})

test_that("MetaHDResult requires a named list", {
  expect_error(MetaHDResult(sets = list(runif(10))), "named list")
})

test_that("MetaHDResult validates truth length", {
  expect_error(
    MetaHDResult(sets = list(A = runif(10)), truth = c(0, 1, 1)),
    "same length"
  )
})
