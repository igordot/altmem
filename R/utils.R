is_mem_output <- function(x) {
  if (!is.list(x)) {
    return(FALSE)
  }

  expected_slots <- c("MEM_matrix", "MAGpop", "MAGref", "IQRpop", "IQRref")
  missing <- setdiff(c(expected_slots, "File Order"), names(x))
  if (length(missing) > 0) {
    return(FALSE)
  }

  n_clusters <- NULL
  n_markers <- NULL

  for (slot in expected_slots) {
    val <- x[[slot]]
    if (!is.list(val) || length(val) != 1) {
      return(FALSE)
    }
    mat <- val[[1]]
    if (!is.matrix(mat)) {
      return(FALSE)
    }
    if (nrow(mat) == 0 || ncol(mat) == 0) {
      return(FALSE)
    }
    if (is.null(n_clusters)) {
      n_clusters <- nrow(mat)
      n_markers <- ncol(mat)
    } else {
      if (nrow(mat) != n_clusters) {
        return(FALSE)
      }
      if (ncol(mat) != n_markers) {
        return(FALSE)
      }
    }
    if (is.null(rownames(mat))) {
      return(FALSE)
    }
  }

  ref_rn <- rownames(x[[expected_slots[[1]]]][[1]])
  for (slot in expected_slots[-1]) {
    if (!identical(rownames(x[[slot]][[1]]), ref_rn)) {
      return(FALSE)
    }
  }

  TRUE
}
