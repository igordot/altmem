#' Marker Enrichment Modeling
#'
#' A wrapper around [cytoMEM::MEM()] that accepts a
#' [SummarizedExperiment][SummarizedExperiment::SummarizedExperiment],
#' [SingleCellExperiment][SingleCellExperiment::SingleCellExperiment], or
#' data frame and supports non-numeric cluster labels.
#'
#' The expression matrix is extracted from the specified assay and transposed
#' to cells-by-markers format. Cluster labels are converted to integer IDs
#' for [cytoMEM::MEM()], then mapped back to the original labels in the
#' output.
#'
#' @param x A [SummarizedExperiment][SummarizedExperiment::SummarizedExperiment],
#'   [SingleCellExperiment][SingleCellExperiment::SingleCellExperiment], or
#'   data frame with a cluster column and numeric marker columns.
#' @param cluster_col Name of the column containing cluster labels. For data
#'   frames, this is ignored; the column must be named `"cluster"`. For
#'   SummarizedExperiment, checked in `rowData(x)`, then `colData(x)`, then
#'   `CATALYST::cluster_codes(x)` (if CATALYST is installed). Labels can be
#'   any type (character, factor, numeric).
#' @param assay_name Name of the assay to use for the expression matrix.
#'   Only used for SummarizedExperiment input. Defaults to `"exprs"`.
#' @param markers Character vector of marker names to include, or `"all"` to
#'   use all markers.
#' @param transform,choose.markers,choose.ref,zero.ref,IQR.thresh Passed to
#'   [cytoMEM::MEM()].
#'
#' @returns The return value of [cytoMEM::MEM()], with row names of the
#'   MEM matrix restored to the original cluster labels.
#'
#' @importFrom cytoMEM MEM
#' @importFrom stats setNames
#' @importFrom SummarizedExperiment assay colData rowData
#' @export
altmem <- function(x,
                   cluster_col,
                   assay_name = "exprs",
                   markers = "all",
                   transform = FALSE,
                   choose.markers = FALSE,
                   choose.ref = FALSE,
                   zero.ref = FALSE,
                   IQR.thresh = NULL) {
  if (is.data.frame(x)) {
    # Rename the cluster column to "cluster" for data frames
    names(x)[names(x) == cluster_col] <- "cluster"
    clusters_original <- x[["cluster"]]
    exprs_mat <- as.matrix(x[, setdiff(names(x), "cluster")])
  } else if (inherits(x, "SummarizedExperiment")) {
    # Find cluster labels and determine orientation
    # Bioconductor convention: features in rows, cells in columns
    # HDCytoData follows cytometry convention: cells in rows, features in columns
    if (cluster_col %in% names(rowData(x)) && nrow(x) > 100) {
      clusters_original <- rowData(x)[[cluster_col]]
      cells_in_rows <- TRUE
    } else if (cluster_col %in% names(colData(x)) && ncol(x) > 100) {
      clusters_original <- colData(x)[[cluster_col]]
      cells_in_rows <- FALSE
    } else if (requireNamespace("CATALYST", quietly = TRUE) &&
      cluster_col %in% names(CATALYST::cluster_codes(x)) && ncol(x) > 100) {
      clusters_original <- CATALYST::cluster_ids(x, cluster_col)
      cells_in_rows <- FALSE
    } else {
      stop("Column '", cluster_col, "' not found in rowData or colData of x.")
    }

    exprs_mat <- assay(x, assay_name)
    if (!cells_in_rows) {
      exprs_mat <- t(exprs_mat)
    }
  } else {
    stop("'x' must be a SummarizedExperiment, SingleCellExperiment, or data frame.")
  }

  # Map original labels to integer IDs
  labels <- sort(unique(as.character(clusters_original)))
  label_to_int <- setNames(seq_along(labels), labels)
  clusters_int <- unname(label_to_int[as.character(clusters_original)])

  # Subset markers (MEM() expects markers to be column indices, not actual names)
  if (length(markers) > 1 && is.character(markers)) {
    exprs_mat <- exprs_mat[, markers]
    markers <- "all"
  }
  mem_df <- as.data.frame(exprs_mat)
  mem_df$cluster <- clusters_int

  mem_res <- cytoMEM::MEM(
    exp_data = mem_df,
    transform = transform,
    choose.markers = choose.markers,
    markers = markers,
    choose.ref = choose.ref,
    zero.ref = zero.ref,
    rename.markers = FALSE,
    file.is.clust = FALSE,
    add.fileID = FALSE,
    IQR.thresh = IQR.thresh
  )

  # Validate MEM output and restore original cluster labels as row names
  if (!is_mem_output(mem_res)) {
    stop("Invalid MEM() output.")
  }

  int_to_label <- setNames(labels, as.character(seq_along(labels)))
  expected_slots <- c("MEM_matrix", "MAGpop", "MAGref", "IQRpop", "IQRref")
  for (slot in expected_slots) {
    rn <- rownames(mem_res[[slot]][[1]])
    rownames(mem_res[[slot]][[1]]) <- unname(int_to_label[rn])
  }

  mem_res
}
