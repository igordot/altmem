test_that("mem_heatmap()", {
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

  pdf(nullfile())
  on.exit(dev.off(), add = TRUE)

  ht <- mem_heatmap(mem_result)
  expect_s4_class(ht, "Heatmap")

  # min_heatmap_score drops low-scoring markers
  ht_filtered <- mem_heatmap(mem_result, min_heatmap_score = 3)
  expect_s4_class(ht_filtered, "Heatmap")
  expect_identical(colnames(ht_filtered@matrix), c("CD4", "CD3"))

  expect_error(mem_heatmap(list(a = 1)))
})
