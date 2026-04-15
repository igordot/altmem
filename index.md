# altmem: marker enrichment modeling

altmem provides an alternative interface to
[cytoMEM](https://doi.org/doi:10.18129/B9.bioc.cytoMEM) marker
enrichment modeling (MEM). The enrichment scores are calculated using
cytoMEM, but the inputs and outputs are adjusted for more flexibility.

Modifications from the original cytoMEM implementation:

- Accepts multiple common data structures (data frame,
  SummarizedExperiment, SingleCellExperiment)
- Supports non-numeric cluster labels
- Generates heatmaps with MEM label annotations and marker filtering
- Returns MEM labels as a named vector with configurable score and
  marker cutoffs

Check the [documentation website](https://igordot.github.io/altmem/) for
more information.
