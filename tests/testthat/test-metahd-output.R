# Methods for viewing/exporting MetaHD results (summary, print)

fit_small_model <- function() {
  d <- make_test_data(K = 6, N = 20)
  model <- MetaHD(Y = d$Y, Slist = d$Slist, method = "multi")
  list(model = model, d = d)
}

test_that("MetaHD returns an object of class MetaHD", {
  fit <- fit_small_model()
  expect_s3_class(fit$model, "MetaHD")
  # classing must not disturb list-based access
  expect_true(is.list(fit$model))
  expect_equal(length(fit$model$estimate), fit$d$N)
})

test_that("summary.MetaHD returns a tidy data frame", {
  fit <- fit_small_model()
  tab <- summary(fit$model, outcome_names = colnames(fit$d$Y))

  expect_s3_class(tab, "data.frame")
  expect_equal(nrow(tab), fit$d$N)
  expect_named(tab, c("outcome", "estimate", "std.err", "pVal", "I2.stat"))
  expect_equal(tab$outcome, colnames(fit$d$Y))
})

test_that("summary.MetaHD falls back to positional outcome names", {
  fit <- fit_small_model()
  tab <- summary(fit$model)
  expect_equal(tab$outcome, paste0("Var", seq_len(fit$d$N)))
})

test_that("summary.MetaHD errors on wrong outcome_names length", {
  fit <- fit_small_model()
  expect_error(
    summary(fit$model, outcome_names = c("a", "b")),
    "must match number of outcomes"
  )
})

test_that("print.MetaHD returns x invisibly and prints an overview", {
  fit <- fit_small_model()

  out <- withVisible(print(fit$model))
  expect_false(out$visible)
  expect_s3_class(out$value, "MetaHD")

  expect_output(print(fit$model), "MetaHD multivariate meta-analysis")
})
