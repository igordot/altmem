#' MEM score heatmap
#'
#' Visualizes MEM enrichment scores as a heatmap using
#' [ComplexHeatmap][ComplexHeatmap::Heatmap]. Rows are clusters and columns
#' are markers. Positive scores (enriched) are shown in red and negative
#' scores (depleted) are shown in blue.
#'
#' @inheritParams mem_labels
#' @param show_mem_labels Whether to annotate rows with MEM labels. When `TRUE`,
#'   labels are generated via [mem_labels()] using `min_label_score`,
#'   `max_label_markers`, and `show_label_scores`. Defaults to `FALSE`.
#' @param show_label_scores Whether to include MEM scores alongside marker names
#'   in the MEM labels. Defaults to `FALSE`.
#' @param min_heatmap_score Minimum absolute MEM score across any cluster for
#'   a marker to be shown as a column in the heatmap. Increase this to remove
#'   low variance markers. Defaults to `0` (show all markers).
#' @param filename Optional file path to save the heatmap. The format is
#'   inferred from the extension (e.g., `".png"`, `".pdf"`, `".svg"`).
#'   Defaults to `NULL` (no file saved).
#' @param width,height Plot dimensions in inches. Only used when `filename` is
#'   provided.
#' @param cluster_rows Whether to cluster rows. Defaults to `TRUE`.
#' @param cluster_columns Whether to cluster columns. Defaults to `TRUE`.
#' @param ... Additional arguments passed to [ComplexHeatmap::Heatmap()].
#'
#' @returns A [ComplexHeatmap::Heatmap-class] object, invisibly.
#'
#' @importFrom circlize colorRamp2
#' @importFrom ComplexHeatmap Heatmap draw max_text_width
#' @importFrom RColorBrewer brewer.pal
#' @export
mem_heatmap <- function(x,
                        show_mem_labels = FALSE,
                        min_label_score = 1,
                        max_label_markers = 5,
                        show_label_scores = FALSE,
                        min_heatmap_score = 0,
                        filename = NULL,
                        width = 15,
                        height = 8,
                        cluster_rows = TRUE,
                        cluster_columns = TRUE,
                        ...) {
  if (!is_mem_output(x)) {
    stop("Invalid MEM output.")
  }

  mem_mat <- x[["MEM_matrix"]][[1]]

  # Drop markers where no cluster exceeds the threshold
  if (min_heatmap_score > 0) {
    max_per_marker <- apply(abs(mem_mat), 2, max)
    keep <- max_per_marker >= min_heatmap_score
    mem_mat <- mem_mat[, keep, drop = FALSE]
  }

  # Diverging color scale symmetric around zero
  max_abs <- max(abs(mem_mat))
  n <- 11
  pal <- rev(grDevices::colorRampPalette(
    RColorBrewer::brewer.pal(n, "RdBu")
  )(n))
  breaks <- seq(-max_abs, max_abs, length.out = n)
  col_fn <- circlize::colorRamp2(breaks, pal)

  # Replace cluster names with MEM-enriched labels
  if (show_mem_labels) {
    lab <- mem_labels(x, min_label_score = min_label_score, max_label_markers = max_label_markers, show_label_scores = show_label_scores)
    rownames(mem_mat) <- lab[rownames(mem_mat)]
  }

  # Clusters × markers heatmap
  ht <- ComplexHeatmap::Heatmap(
    mem_mat,
    name = "MEM\nScore",
    col = col_fn,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    row_names_max_width = ComplexHeatmap::max_text_width(rownames(mem_mat)),
    row_title = NULL,
    column_title = "Marker Enrichment Modeling (MEM)",
    ...
  )

  # Save to file if requested
  if (!is.null(filename)) {
    ext <- tolower(tools::file_ext(filename))
    device <- switch(ext,
      png = \(...) grDevices::png(..., units = "in", res = 300),
      pdf = grDevices::pdf,
      svg = grDevices::svg,
      stop("Unsupported file extension: ", ext)
    )
    device(filename, width = width, height = height)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  # Extra right padding to prevent legend clipping
  ComplexHeatmap::draw(ht, padding = grid::unit(c(2, 2, 2, 10), "mm"))
  invisible(ht)
}
