# Heatmap of MEM scores

Visualizes MEM enrichment scores as a heatmap.

## Usage

``` r
mem_heatmap(
  x,
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
  ...
)
```

## Arguments

- x:

  The return value of
  [`altmem()`](https://igordot.github.io/altmem/reference/altmem.md) or
  [`cytoMEM::MEM()`](https://rdrr.io/pkg/cytoMEM/man/MEM.html).

- show_mem_labels:

  Whether to annotate rows with MEM labels. When `TRUE`, labels are
  generated via
  [`mem_labels()`](https://igordot.github.io/altmem/reference/mem_labels.md)
  using `min_label_score`, `max_label_markers`, and `show_label_scores`.
  Defaults to `FALSE`.

- min_label_score:

  Minimum absolute MEM score for a marker to be included in the label.
  Defaults to `1`.

- max_label_markers:

  Maximum number of markers to include per direction (positive and
  negative) in each label. Defaults to `5`.

- show_label_scores:

  Whether to include MEM scores alongside marker names in the MEM
  labels. Defaults to `FALSE`.

- min_heatmap_score:

  Minimum absolute MEM score across any cluster for a marker to be shown
  as a column in the heatmap. Increase this to remove low variance
  markers. Defaults to `0` (show all markers).

- filename:

  Optional file path to save the heatmap. The format is inferred from
  the extension (e.g., `".png"`, `".pdf"`, `".svg"`). Defaults to `NULL`
  (no file saved).

- width, height:

  Plot dimensions in inches. Only used when `filename` is provided.

- cluster_rows:

  Whether to cluster rows. Defaults to `TRUE`.

- cluster_columns:

  Whether to cluster columns. Defaults to `TRUE`.

- ...:

  Additional arguments passed to
  [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html).

## Value

A
[ComplexHeatmap::Heatmap](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap-class.html)
object, invisibly.
