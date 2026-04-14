# altmem

<!-- badges: start -->
[![R-CMD-check](https://github.com/igordot/altmem/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/igordot/altmem/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/igordot/altmem/graph/badge.svg)](https://app.codecov.io/gh/igordot/altmem)
<!-- badges: end -->

altmem provides an alternative interface to [cytoMEM](https://doi.org/doi:10.18129/B9.bioc.cytoMEM) marker enrichment modeling (MEM) R package.
The enrichment scores are calculated using cytoMEM, but the inputs and outputs are adjusted for easier integration with other tools.

Modifications from the cytoMEM implementation:

- Accepts common data structures (data frame, SummarizedExperiment, SingleCellExperiment)
- Supports non-numeric cluster labels
- Generates heatmaps with MEM label annotations and marker filtering
- Returns MEM labels as a named vector with configurable score and marker cutoffs

You can install altmem from GitHub using:

```r
BiocManager::install("igordot/altmem")
```
