test_that("compute_clustering_diagnostics returns correct structure", {
  d <- make_test_data()
  cormat <- estimateCorMat(d$Y)
  distMat <- as.dist(1 - abs(cormat))
  hc <- hclust(distMat, method = "complete")
  clusters <- cutree(hc, k = 3)
  
  diag <- compute_clustering_diagnostics(clusters, distMat)
  
  # Check all expected elements present
  expect_true(all(c("avg_silhouette",
                    "cluster_avg",
                    "n_clusters",
                    "cluster_sizes",
                    "silhouette_object") %in% names(diag)))
  
  # Check avg silhouette is in valid range
  expect_true(diag$avg_silhouette >= -1 && diag$avg_silhouette <= 1)
  
  # Check n_clusters matches
  expect_equal(diag$n_clusters, 3)
  
  # Check cluster sizes sum to N
  expect_equal(sum(diag$cluster_sizes), d$N)
})

test_that("silhouette object can be plotted", {
  d <- make_test_data()
  cormat <- estimateCorMat(d$Y)
  distMat <- as.dist(1 - abs(cormat))
  hc <- hclust(distMat, method = "complete")
  clusters <- cutree(hc, k = 3)
  
  diag <- compute_clustering_diagnostics(clusters, distMat)
  
  expect_s3_class(diag$silhouette_object, "silhouette")
  expect_no_error(plot(diag$silhouette_object))
})