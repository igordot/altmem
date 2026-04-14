# Create MEM-enriched cluster labels

Generates descriptive labels for each cluster that include measured
features specifically enriched on the population based on MEM scores,
following the cytoMEM format. Markers are ranked by absolute MEM score
and separated into positively and negatively enriched groups.

## Usage

``` r
mem_labels(
  x,
  min_label_score = 1,
  max_label_markers = 5,
  show_label_scores = TRUE
)
```

## Arguments

- x:

  The return value of
  [`altmem()`](https://igordot.github.io/altmem/reference/altmem.md) or
  [`cytoMEM::MEM()`](https://rdrr.io/pkg/cytoMEM/man/MEM.html).

- min_label_score:

  Minimum absolute MEM score for a marker to be included in the label.
  Defaults to `1`.

- max_label_markers:

  Maximum number of markers to include per direction (positive and
  negative) in each label. Defaults to `5`.

- show_label_scores:

  Whether to include MEM scores alongside marker names. Defaults to
  `TRUE`.

## Value

A named character vector.
