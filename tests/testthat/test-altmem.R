test_that("altmem() with HDCytoData Samusik_01", {
  # CyTOF dataset containing 39 markers with 24 gated cell population labels
  se <- suppressMessages(HDCytoData::Samusik_01_SE())

  # Remove unassigned cells and downsample large populations
  # table(rowData(se)$population_id)
  se <- se[rowData(se)$population_id != "unassigned", ]
  set.seed(123)
  idx <- unlist(tapply(seq_len(nrow(se)), rowData(se)$population_id, \(i) {
    if (length(i) > 500) sample(i, 500) else i
  }))
  se <- se[idx, ]
  # Subset to cell type markers only
  se <- se[, colData(se)$marker_class == "type"]
  # Transform data using asinh with cofactor 5
  assay(se) <- asinh(assay(se) / 5)

  # Generare a data frame version of the input
  exprs_df <- as.data.frame(assay(se, "exprs"))
  exprs_df$population_id <- rowData(se)$population_id

  # Run MEM
  res <- altmem(se, cluster_col = "population_id")

  # Returns a list with expected MEM slots
  expect_type(res, "list")
  expect_named(
    res,
    c("MEM_matrix", "MAGpop", "MAGref", "IQRpop", "IQRref", "File Order"),
    ignore.order = TRUE
  )

  # Row names of MEM_matrix match original cluster labels
  labels <- sort(unique(as.character(rowData(se)$population_id)))
  expect_equal(length(labels), 24)
  expect_equal(sort(rownames(res$MEM_matrix[[1]])), labels)

  # MEM_matrix dimensions: clusters x markers
  expect_equal(nrow(res$MEM_matrix[[1]]), 24)
  expect_gt(ncol(res$MEM_matrix[[1]]), 0)

  # Run MEM using the data frame input
  res_df <- altmem(exprs_df, cluster_col = "population_id")

  # Returns the same list list
  expect_type(res_df, "list")
  expect_equal(names(res_df), names(res))
  expect_equal(res_df$MEM_matrix[[1]], res$MEM_matrix[[1]])
})
