#' Create MEM-enriched cluster labels
#'
#' Generates descriptive labels for each cluster that include measured features
#' specifically enriched on the population based on MEM scores, following the
#' cytoMEM format. Markers are ranked by absolute MEM score and separated into
#' positively and negatively enriched groups.
#'
#' @param mem_result The return value of [altmem()] or [cytoMEM::MEM()].
#' @param min_score Minimum absolute MEM score for a marker to be
#'   included in the label. Defaults to `1`.
#' @param max_markers Maximum number of markers to include per direction
#'   (positive and negative) in each label. Defaults to `5`.
#'
#' @returns A named character vector.
#'
#' @importFrom stats setNames
#' @export
mem_labels <- function(mem_result, min_score = 1, max_markers = 5) {
  if (!is_mem_output(mem_result)) {
    stop("Invalid MEM output.")
  }

  mem_mat <- mem_result[["MEM_matrix"]][[1]]
  cluster_names <- rownames(mem_mat)
  labels <- character(nrow(mem_mat))

  for (i in seq_len(nrow(mem_mat))) {
    scores <- mem_mat[i, ]

    pos <- scores[scores >= min_score]
    neg <- scores[scores < 0 & abs(scores) >= min_score]

    pos <- pos[order(pos, decreasing = TRUE)]
    neg <- neg[order(abs(neg), decreasing = TRUE)]

    pos <- pos[seq_len(min(length(pos), max_markers))]
    neg <- neg[seq_len(min(length(neg), max_markers))]

    if (length(pos) > 0) {
      pos_str <- paste(
        paste0(names(pos), "+", abs(round(pos, 1))),
        collapse = " "
      )
    } else {
      pos_str <- "None"
    }

    if (length(neg) > 0) {
      neg_str <- paste(
        paste0(names(neg), "-", abs(round(neg, 1))),
        collapse = " "
      )
    } else {
      neg_str <- "None"
    }

    labels[[i]] <- paste(cluster_names[[i]], ":", "\u25b2", pos_str, "\u25bc", neg_str)
  }

  stats::setNames(labels, cluster_names)
}
