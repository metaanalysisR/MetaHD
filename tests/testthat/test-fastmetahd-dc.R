test_that("fastMetaHD DC with dynamicTreeCut returns diagnostics", {
  d <- make_test_data()
  model <- MetaHD(Y = d$Y,
                  Slist = d$Slist,
                  method = "multi",
                  useDivideConquer = TRUE,
                  dendro.method = "dynamicTreeCut")
  
  expect_false(is.null(model$clustering_diagnostics))
  expect_true(model$clustering_diagnostics$n_clusters >= 1)
  
  # Only check silhouette range if more than 1 cluster found
  if (model$clustering_diagnostics$n_clusters > 1) {
    expect_true(model$clustering_diagnostics$avg_silhouette >= -1)
    expect_true(model$clustering_diagnostics$avg_silhouette <= 1)
  } else {
    expect_true(is.na(model$clustering_diagnostics$avg_silhouette))
  }
})

test_that("fastMetaHD DC with optimalK returns silhouette_by_k", {
  d <- make_test_data()
  model <- MetaHD(Y = d$Y,
                  Slist = d$Slist,
                  method = "multi",
                  useDivideConquer = TRUE,
                  dendro.method = "optimalK")
  
  expect_false(is.null(model$clustering_diagnostics$silhouette_by_k))
  expect_true(nrow(model$clustering_diagnostics$silhouette_by_k) >= 1)
})

test_that("fastMetaHD DC with user DCgroups returns NULL diagnostics", {
  d <- make_test_data()
  
  outcome_names <- colnames(d$Y)
  half <- floor(length(outcome_names) / 2)
  
  groups <- list(
    g1 = outcome_names[seq_len(half)],
    g2 = outcome_names[seq(half + 1, length(outcome_names))]
  )
  
  model <- MetaHD(Y = d$Y,
                  Slist = d$Slist,
                  method = "multi",
                  useDivideConquer = TRUE,
                  DCgroups = groups)
  
  expect_null(model$clustering_diagnostics)
})

test_that("fastMetaHD DC results consistent across dendro methods", {
  d <- make_test_data()
  
  model_dtc <- MetaHD(Y = d$Y,
                      Slist = d$Slist,
                      method = "multi",
                      useDivideConquer = TRUE,
                      dendro.method = "dynamicTreeCut")
  
  model_k3 <- MetaHD(Y = d$Y,
                     Slist = d$Slist,
                     method = "multi",
                     useDivideConquer = TRUE,
                     dendro.method = "fixedK",
                     dendro.k = 3)
  
  # Both should return estimates of length N
  expect_equal(length(model_dtc$estimate), d$N)
  expect_equal(length(model_k3$estimate),  d$N)
})

test_that("divide-and-conquer returns per-outcome I2.stat aligned with outcomes", {
  d <- make_test_data(K = 6, N = 20)

  full <- MetaHD(Y = d$Y, Slist = d$Slist, method = "multi")

  # A single DC group spanning all outcomes (in column order) reproduces the
  # full multivariate fit, so the per-outcome statistics must match exactly.
  one_group <- list(all = colnames(d$Y))
  dc <- MetaHD(Y = d$Y, Slist = d$Slist, method = "multi",
               useDivideConquer = TRUE, DCgroups = one_group)

  # I2.stat must be one value per outcome 
  expect_equal(length(dc$I2.stat), d$N)
  expect_equal(unname(dc$estimate), unname(full$estimate))

  # full$I2.stat carries a leading overall multivariate I^2; the remaining
  # per-outcome values must align with the divide-and-conquer output.
  # Compare values only -- the two paths differ in whether names are attached.
  expect_equal(unname(dc$I2.stat), unname(full$I2.stat[-1]))
})