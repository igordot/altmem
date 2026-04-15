# Marker Enrichment Modeling

A wrapper around
[`cytoMEM::MEM()`](https://rdrr.io/pkg/cytoMEM/man/MEM.html) that
accepts a matrix, data frame,
[SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html),
or
[SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html),
and supports non-numeric cluster labels.

## Usage

``` r
altmem(
  x,
  cluster_col,
  assay_name = "exprs",
  markers = "all",
  transform = FALSE,
  cofactor = 1,
  choose.markers = FALSE,
  choose.ref = FALSE,
  zero.ref = FALSE,
  IQR.thresh = NULL
)
```

## Arguments

- x:

  A matrix, data frame,
  [SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html),
  or
  [SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html).

- cluster_col:

  Name of the column containing cluster labels.

- assay_name:

  Name of the assay to use for the expression matrix. Only used for
  SummarizedExperiment input. Defaults to `"exprs"`.

- markers:

  Character vector of marker names to include, or `"all"` to use all
  markers.

- transform, cofactor, choose.markers, choose.ref, zero.ref, IQR.thresh:

  Passed to
  [`cytoMEM::MEM()`](https://rdrr.io/pkg/cytoMEM/man/MEM.html).

## Value

A list of matrices (the return value of
[`cytoMEM::MEM()`](https://rdrr.io/pkg/cytoMEM/man/MEM.html)).
