setup_res <- function() {
  set.seed(123)
  N <- 100
  truth <- rbinom(N, 1, 0.2)
  MetaHDResult(
    sets = list(
      method_A = runif(N)^ifelse(truth, 5, 1),
      method_B = runif(N)^ifelse(truth, 3, 1),
      method_C = runif(N)^ifelse(truth, 2, 1)
    ),
    truth = truth
  )
}

setup_res_no_truth <- function() {
  set.seed(123)
  N <- 100
  MetaHDResult(
    sets = list(
      method_A = runif(N)^5,
      method_B = runif(N)^3,
      method_C = runif(N)^2
    )
  )
}

# ---- Input validation ----

test_that("plot errors for unknown type", {
  res <- setup_res()
  expect_error(plot(res, type = "unknown"), "should be one of")
})

test_that("ROC plot errors when truth is NULL", {
  res <- setup_res_no_truth()
  expect_error(
    plot(res, type = "ROC"),
    "True signals must be provided"
  )
})

test_that("ROC plot errors when roc.colors length does not match methods", {
  res <- setup_res()
  expect_error(
    plot(res, type = "ROC", roc.colors = c("red", "blue")),
    "Length of `roc.colors` must match number of methods"
  )
})

# ---- Example 1: simple upset call ----

test_that("simple upset plot runs without error", {
  res <- setup_res()
  expect_no_error(plot(res, type = "upset"))
})

test_that("upset returns NULL invisibly", {
  res <- setup_res()
  result <- plot(res, type = "upset")
  expect_null(result)
})

test_that("upset with show.truth = FALSE runs without error", {
  res <- setup_res()
  expect_no_error(plot(res, type = "upset", show.truth = FALSE))
})

test_that("upset without truth runs without error", {
  res <- setup_res_no_truth()
  expect_no_error(plot(res, type = "upset"))
})

# ---- Example 2: highlight argument ----

test_that("upset with highlight runs without error", {
  res <- setup_res()
  expect_no_error(
    plot(res, type = "upset",
         highlight = list(
           c("Truth", "method_A", "method_B", "method_C")
         ),
         highlight.colors = "darkgreen")
  )
})

# ---- Example 3: queries argument ----

test_that("upset with intersection query runs without error", {
  res <- setup_res()
  expect_no_error(
    plot(res, type = "upset",
         queries = list(
           list(
             query = UpSetR::intersects,
             params = list(c("method_A", "method_B", "method_C")),
             color = "dodgerblue3",
             active = TRUE,
             query.name = "Identified by all methods"
           ),
           list(
             query = function(row) {
               row["method_A"] == 1 && sum(row) < length(row)
             },
             color = "orange",
             active = TRUE,
             query.name = "Others identified by method_A"
           )
         ),
         show.truth = FALSE)
  )
})

# ---- Example 4: ROC plot and AUCs ----

test_that("ROC plot runs without error", {
  res <- setup_res()
  expect_no_error(plot(res, type = "ROC"))
})

test_that("ROC plot returns named numeric AUC vector", {
  res <- setup_res()
  aucs <- plot(res, type = "ROC")
  
  expect_true(is.numeric(aucs))
  expect_equal(length(aucs), 3)
  expect_equal(names(aucs), c("method_A", "method_B", "method_C"))
  expect_true(all(aucs >= 0 & aucs <= 1))
})

test_that("AUCs returned invisibly", {
  res <- setup_res()
  result <- withVisible(plot(res, type = "ROC"))
  expect_false(result$visible)
  expect_true(is.numeric(result$value))
})

test_that("ROC with custom colors runs without error", {
  res <- setup_res()
  expect_no_error(
    plot(res, type = "ROC",
         roc.colors = c("red", "blue", "green"))
  )
})

# ---- Venn diagram (requires gVenn, R >= 4.5) ----

test_that("venn plot runs and returns named list invisibly", {
  skip_if_not_installed("gVenn")
  skip_if(getRversion() < "4.5")
  res <- setup_res()
  path <- tempfile(fileext = ".pdf")
  out <- withVisible(plot(res, type = "venn", file = path))
  expect_false(out$visible)
  expect_type(out$value, "list")
  expect_named(out$value, c("method_A", "method_B", "method_C"))
  expect_true(file.exists(path))
  unlink(path)
})

# ---- Example 5: save to file ----

test_that("upset saves to PDF without error", {
  res <- setup_res()
  path <- tempfile(fileext = ".pdf")
  
  expect_no_error(
    plot(res, type = "upset", file = path, width = 10, height = 6)
  )
  expect_true(file.exists(path))
  expect_gt(file.size(path), 0)
  unlink(path)
})

test_that("ROC saves to PNG without error", {
  res <- setup_res()
  path <- tempfile(fileext = ".png")
  
  expect_no_error(
    plot(res, type = "ROC", file = path, width = 8, height = 6, dpi = 600)
  )
  expect_true(file.exists(path))
  expect_gt(file.size(path), 0)
  unlink(path)
})

test_that("ROC saves to SVG without error", {
  res <- setup_res()
  path <- tempfile(fileext = ".svg")
  
  expect_no_error(
    plot(res, type = "ROC", file = path)
  )
  expect_true(file.exists(path))
  expect_gt(file.size(path), 0)
  unlink(path)
})

# ---- plot_dendrogram ----

test_that("plot_dendrogram runs without error", {
  d <- make_test_data()
  expect_no_error(plot_dendrogram(d$Y))
})

test_that("plot_dendrogram returns hclust object invisibly", {
  d <- make_test_data()
  result <- withVisible(plot_dendrogram(d$Y))
  expect_false(result$visible)
  expect_s3_class(result$value, "hclust")
})







