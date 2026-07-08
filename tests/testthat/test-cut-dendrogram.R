test_that("dynamicTreeCut returns valid cluster labels", {
  d <- make_test_data()
  cormat <- estimateCorMat(d$Y)
  distMat <- as.dist(1 - abs(cormat))
  hc <- hclust(distMat, method = "complete")
  
  result <- cut_dendrogram(hc, distMat, method = "dynamicTreeCut")
  
  expect_true(is.numeric(result$clusters) || is.integer(result$clusters))
  expect_equal(length(result$clusters), d$N)
  expect_true(length(unique(result$clusters)) >= 1)
})

test_that("fixedK returns exactly k clusters", {
  d <- make_test_data()
  cormat <- estimateCorMat(d$Y)
  distMat <- as.dist(1 - abs(cormat))
  hc <- hclust(distMat, method = "complete")
  
  for (k in c(2, 3, 5)) {
    result <- cut_dendrogram(hc, distMat, method = "fixedK", k = k)
    expect_equal(length(unique(result$clusters)), k,
                 info = paste("Failed for k =", k))
  }
})

test_that("fixedK errors when k not provided", {
  d <- make_test_data()
  cormat <- estimateCorMat(d$Y)
  distMat <- as.dist(1 - abs(cormat))
  hc <- hclust(distMat, method = "complete")
  
  expect_error(
    cut_dendrogram(hc, distMat, method = "fixedK"),
    "Please provide 'k'"
  )
})

test_that("fixedHeight returns valid clusters", {
  d <- make_test_data()
  cormat <- estimateCorMat(d$Y)
  distMat <- as.dist(1 - abs(cormat))
  hc <- hclust(distMat, method = "complete")
  
  # Use median height of dendrogram
  h <- median(hc$height)
  result <- cut_dendrogram(hc, distMat, method = "fixedHeight", height = h)
  
  expect_true(is.numeric(result$clusters) || is.integer(result$clusters))
  expect_equal(length(result$clusters), d$N)
})

test_that("fixedHeight errors when height not provided", {
  d <- make_test_data()
  cormat <- estimateCorMat(d$Y)
  distMat <- as.dist(1 - abs(cormat))
  hc <- hclust(distMat, method = "complete")
  
  expect_error(
    cut_dendrogram(hc, distMat, method = "fixedHeight"),
    "Please provide 'height'"
  )
})

test_that("optimalK returns valid clusters and silhouette_by_k", {
  d <- make_test_data()
  cormat <- estimateCorMat(d$Y)
  distMat <- as.dist(1 - abs(cormat))
  hc <- hclust(distMat, method = "complete")
  
  result <- cut_dendrogram(hc, distMat, method = "optimalK")
  
  expect_true(is.numeric(result$clusters) || is.integer(result$clusters))
  expect_equal(length(result$clusters), d$N)
  expect_true(length(unique(result$clusters)) >= 2)
  expect_true(!is.null(result$silhouette_by_k))
  expect_true(all(c("k", "avg_silhouette") %in% names(result$silhouette_by_k)))
})