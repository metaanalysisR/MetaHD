# Tests for the heatmap visualisations (ComplexHeatmap-based)

# ---- plot_effect_heatmap ----

test_that("plot_effect_heatmap returns a Heatmap object invisibly", {
  set.seed(1)
  est <- list(A = rnorm(20), B = rnorm(20), C = rnorm(20))
  path <- tempfile(fileext = ".pdf")

  out <- withVisible(plot_effect_heatmap(est, file = path))
  expect_false(out$visible)
  expect_s4_class(out$value, "Heatmap")
  expect_true(file.exists(path))
  unlink(path)
})

test_that("plot_effect_heatmap supports row_km partitioning", {
  set.seed(1)
  est <- list(A = rnorm(60), B = rnorm(60))
  path <- tempfile(fileext = ".pdf")

  expect_no_error(plot_effect_heatmap(est, row_km = 3, file = path))
  unlink(path)
})

test_that("plot_effect_heatmap validates row_km", {
  est <- list(A = rnorm(10), B = rnorm(10))
  expect_error(plot_effect_heatmap(est, row_km = 0), "row_km")
})

test_that("plot_effect_heatmap errors on unequal vector lengths", {
  est <- list(A = rnorm(10), B = rnorm(8))
  expect_error(plot_effect_heatmap(est), "same length")
})

test_that("plot_effect_heatmap supports label_top_n", {
  set.seed(1)
  est <- list(A = rnorm(60), B = rnorm(60))
  path <- tempfile(fileext = ".pdf")

  expect_no_error(plot_effect_heatmap(est, label_top_n = 5, file = path))
  unlink(path)
})

test_that("plot_effect_heatmap validates label_top_n", {
  est <- list(A = rnorm(10), B = rnorm(10))
  expect_error(plot_effect_heatmap(est, label_top_n = 0), "label_top_n")
})

# ---- plot_correlation_heatmap ----

test_that("plot_correlation_heatmap returns a Heatmap object invisibly", {
  d <- make_test_data(K = 5, N = 20)
  path <- tempfile(fileext = ".pdf")

  out <- withVisible(plot_correlation_heatmap(d$Slist[[1]], file = path))
  expect_false(out$visible)
  expect_s4_class(out$value, "Heatmap")
  expect_true(file.exists(path))
  unlink(path)
})

test_that("plot_correlation_heatmap supports symmetric row_km", {
  set.seed(1)
  d <- make_test_data(K = 5, N = 20)
  path <- tempfile(fileext = ".pdf")

  expect_no_error(
    plot_correlation_heatmap(d$Slist[[1]], row_km = 3, file = path)
  )
  unlink(path)
})

test_that("plot_correlation_heatmap validates its inputs", {
  d <- make_test_data(K = 5, N = 20)
  
  expect_error(
    plot_correlation_heatmap(d$Slist),
    "covariance or correlation matrix"
  )
  # non-square matrix
  expect_error(
    plot_correlation_heatmap(d$Slist[[1]][, 1:5]),
    "square"
  )
  expect_error(
    plot_correlation_heatmap(d$Slist[[1]], row_km = 0),
    "row_km"
  )
})

test_that("plot_correlation_heatmap accepts a correlation matrix directly", {
  d <- make_test_data(K = 5, N = 20)
  cormat <- stats::cov2cor(d$Slist[[1]])
  path <- tempfile(fileext = ".pdf")

  expect_no_error(
    plot_correlation_heatmap(cormat, is.corr = TRUE, file = path)
  )
  unlink(path)
})

