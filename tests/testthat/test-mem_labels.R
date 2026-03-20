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

  labels <- mem_labels(mem_result, min_label_score = 0)
  expect_type(labels, "character")
  expect_named(labels, c("clusterA", "clusterB"))
  expect_match(labels[["clusterA"]], "clusterA")
  expect_match(labels[["clusterA"]], "CD4")
  expect_match(labels[["clusterA"]], "CD3")

  # min_score filters low-scoring markers
  labels <- mem_labels(mem_result, min_label_score = 3)
  expect_no_match(labels[["clusterA"]], "CD8")
  expect_match(labels[["clusterA"]], "CD4")

  # max_label_markers caps markers per direction
  labels <- mem_labels(mem_result, min_label_score = 0, max_label_markers = 1)
  expect_match(labels[["clusterA"]], "CD3")
  expect_no_match(labels[["clusterA"]], "CD4")
  expect_match(labels[["clusterB"]], "CD4")
  expect_no_match(labels[["clusterB"]], "CD3")

  # show_label_scores = FALSE omits numeric values
  labels <- mem_labels(mem_result, min_label_score = 0, show_label_scores = FALSE)
  expect_match(labels[["clusterA"]], "CD4")
  expect_no_match(labels[["clusterA"]], "\\+5")

  # errors on invalid input
  expect_error(mem_labels(list(a = 1)))
})
