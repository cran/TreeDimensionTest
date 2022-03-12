---
title: The 'TreeDimensionTest' R package
bibliography: inst/REFERENCES.bib
---
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

### Overview

The package offers tests for trajectory presence and heterogeneity on multivariate data. Two statistical methods (Tenha and Song, 2022) are implemented. The tree dimension test quantifies the statistical evidence for trajectory presence. The subset specificity measure summarizes pattern heterogeneity using the minimum subtree cover. There is no user tunable parameters for either method. Examples are included to illustrate how to use the methods on single-cell data for studying gene and pathway expression dynamics and pathway expression specificity.
 
### To download and install the package

```{r}
install.packages("DimensionTreeTest")
```

### Examples

See the package vignette.

### Citing the package

Tenha L, Song M (2022). “Inference of trajectory presence by
tree dimension and subset specificity by subtree cover.” _PLOS
Computational Biology_, *18*(2), e1009829. doi:
10.1371/journal.pcbi.1009829 (URL:
<https://doi.org/10.1371/journal.pcbi.1009829>).
