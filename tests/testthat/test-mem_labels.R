test_that("mem_labels()", {
  mem_mat <- matrix(
    c(5, -3, 0.5, 2, 8, -1),
    nrow = 2,
    dimnames = list(c("clusterA", "clusterB"), c("CD4", "CD8", "CD3"))
  )
  mem_result <- list(
    MEM_matrix = list(mem_mat),
    MAGpop = list(mem_mat),
    MAGref = list(mem_mat),
    IQRpop = list(mem_mat),
    IQRref = list(mem_mat),
    `File Order` = "test"
  )

  labels <- mem_labels(mem_result, min_score = 0)
  expect_type(labels, "character")
  expect_named(labels, c("clusterA", "clusterB"))
  expect_match(labels[["clusterA"]], "clusterA")
  expect_match(labels[["clusterA"]], "CD4")
  expect_match(labels[["clusterA"]], "CD3")

  # min_score filters low-scoring markers
  labels <- mem_labels(mem_result, min_score = 3)
  expect_no_match(labels[["clusterA"]], "CD8")
  expect_match(labels[["clusterA"]], "CD4")

  # max_markers caps markers per direction
  labels <- mem_labels(mem_result, min_score = 0, max_markers = 1)
  expect_match(labels[["clusterA"]], "CD3")
  expect_no_match(labels[["clusterA"]], "CD4")
  expect_match(labels[["clusterB"]], "CD4")
  expect_no_match(labels[["clusterB"]], "CD3")

  # errors on invalid input
  expect_error(mem_labels(list(a = 1)))
})
