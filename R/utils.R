is_mem_output <- function(x) {
  if (!is.list(x)) {
    message("Invalid MEM output: not a list")
    return(FALSE)
  }

  expected_slots <- c("MEM_matrix", "MAGpop", "MAGref", "IQRpop", "IQRref")
  missing <- setdiff(c(expected_slots, "File Order"), names(x))
  if (length(missing) > 0) {
    message("Invalid MEM output: missing slots: ", toString(missing))
    return(FALSE)
  }

  n_clusters <- NULL
  n_markers <- NULL

  for (slot in expected_slots) {
    val <- x[[slot]]
    if (!is.list(val) || length(val) != 1) {
      message("Invalid MEM output: ", slot, " is not a length-1 list")
      return(FALSE)
    }
    mat <- val[[1]]
    if (!is.matrix(mat)) {
      message("Invalid MEM output: ", slot, "[[1]] is not a matrix")
      return(FALSE)
    }
    if (nrow(mat) == 0 || ncol(mat) == 0) {
      message("Invalid MEM output: ", slot, "[[1]] has zero rows or columns")
      return(FALSE)
    }
    if (is.null(n_clusters)) {
      n_clusters <- nrow(mat)
      n_markers <- ncol(mat)
    } else {
      if (nrow(mat) != n_clusters) {
        message("Invalid MEM output: ", slot, " has ", nrow(mat), " rows, expected ", n_clusters)
        return(FALSE)
      }
      if (ncol(mat) != n_markers) {
        message("Invalid MEM output: ", slot, " has ", ncol(mat), " columns, expected ", n_markers)
        return(FALSE)
      }
    }
    if (is.null(rownames(mat))) {
      message("Invalid MEM output: ", slot, "[[1]] has no row names")
      return(FALSE)
    }
  }

  ref_rn <- rownames(x[[expected_slots[[1]]]][[1]])
  for (slot in expected_slots[-1]) {
    if (!identical(rownames(x[[slot]][[1]]), ref_rn)) {
      message("Invalid MEM output: row names of ", slot, " differ from MEM_matrix")
      return(FALSE)
    }
  }

  TRUE
}
