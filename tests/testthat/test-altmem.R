suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(HDCytoData))

test_that("altmem works with Samusik_01_SE from HDCytoData", {
  # CyTOF dataset containing 39 markers with 24 gated cell population labels
  se <- suppressMessages(Samusik_01_SE())
  # Remove unassigned cells
  se <- se[rowData(se)$population_id != "unassigned", ]
  # Subset object to speed up testing
  set.seed(123)
  idx <- sample(nrow(se), 20000)
  se <- se[idx, ]
  # Subset to cell type markers only
  se <- se[, colData(se)$marker_class == "type"]
  # Transform data using asinh with cofactor 5
  assay(se) <- asinh(assay(se) / 5)

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
})
