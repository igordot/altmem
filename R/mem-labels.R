#' MEM-enriched cluster labels
#'
#' Generates descriptive labels for each cluster that include measured features
#' specifically enriched on the population based on MEM scores, following the
#' cytoMEM format. Markers are ranked by absolute MEM score and separated into
#' positively and negatively enriched groups.
#'
#' @param x The return value of [altmem()] or [cytoMEM::MEM()].
#' @param min_label_score Minimum absolute MEM score for a marker to be
#'   included in the label. Defaults to `1`.
#' @param max_label_markers Maximum number of markers to include per direction
#'   (positive and negative) in each label. Defaults to `5`.
#' @param show_label_scores Whether to include MEM scores alongside marker names.
#'   Defaults to `TRUE`.
#'
#' @returns A named character vector.
#'
#' @importFrom stats setNames
#' @importFrom utils head
#' @export
mem_labels <- function(x, min_label_score = 1, max_label_markers = 5, show_label_scores = TRUE) {
  if (!is_mem_output(x)) {
    stop("Invalid MEM output.")
  }

  mem_mat <- x[["MEM_matrix"]][[1]]
  cluster_names <- rownames(mem_mat)
  labels <- character(nrow(mem_mat))

  for (i in seq_len(nrow(mem_mat))) {
    marker_values <- mem_mat[i, ]

    pos <- marker_values[marker_values >= min_label_score]
    neg <- marker_values[marker_values < 0 & abs(marker_values) >= min_label_score]

    pos <- pos[order(pos, decreasing = TRUE)]
    neg <- neg[order(abs(neg), decreasing = TRUE)]

    pos <- head(pos, max_label_markers)
    neg <- head(neg, max_label_markers)

    if (length(pos) > 0) {
      if (show_label_scores) {
        pos_str <- paste(paste0(names(pos), "+", abs(round(pos, 1))), collapse = " ")
      } else {
        pos_str <- paste(names(pos), collapse = " ")
      }
    } else {
      pos_str <- "None"
    }

    if (length(neg) > 0) {
      if (show_label_scores) {
        neg_str <- paste(paste0(names(neg), "-", abs(round(neg, 1))), collapse = " ")
      } else {
        neg_str <- paste(names(neg), collapse = " ")
      }
    } else {
      neg_str <- "None"
    }

    labels[[i]] <- paste(cluster_names[[i]], ":", "\u25b2", pos_str, "\u25bc", neg_str)
  }

  stats::setNames(labels, cluster_names)
}
